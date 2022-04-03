#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 14:56:35 2020
Updated Wed Apr 23 13:20:00 2020
Updated on Sat Apr  2 18:57:50 2022

@author: cghiaus
"""
import numpy as np
import pandas as pd
import psychro as psy
import matplotlib.pyplot as plt

# global variables
# UA = 935.83                 # bldg conductance
# θIsp, wIsp = 18, 6.22e-3    # indoor conditions

θOd = -1                    # outdoor design conditions
mid = 2.18                  # infiltration design

# constants
c = 1e3                     # air specific heat J/kg K
l = 2496e3                  # latent heat J/kg


# *****************************************
# RECYCLED AIR
# *****************************************
def ModelRecAir(m, α, β, θS, θIsp, φIsp, θO, φO, Qsa, Qla, mi, UA):
    """
    Model:
        Heating and adiabatic humidification
        Recycled air
        CAV Constant Air Volume:
            mass flow rate calculated for design conditions
            maintained constant in all situations

    INPUTS:
        m       mass flow of supply dry air, kg/s
        α       mixing ratio of outdoor air, -
        β       by-pass factor of the adiabatic humidifier, -
        θS      supply air, °C
        θIsp    indoor air setpoint, °C
        φIsp    indoor relative humidity set point, -
        θO      outdoor temperature for design, °C
        φO      outdoor relative humidity for design, -
        Qsa     aux. sensible heat, W
        Qla     aux. latente heat, W
        mi      infiltration massflow rate, kg/s
        UA      global conductivity bldg, W/K

    System:
        MX1:    Mixing box
        HC1:    Heating Coil
        AH:     Adiabatic Humidifier
        MX2:    Mixing in humidifier model
        HC2:    Reheating coil
        TZ:     Thermal Zone
        BL:     Building
        Kw:     Controller - humidity
        Kθ:     Controller - temperature
        o:      outdoor conditions
        0..5    unknown points (temperature, humidity ratio)

        <----|<-------------------------------------------|
             |                                            |
             |              |-------|                     |
        -o->MX1--0->HC1--1->|       MX2--3->HC2--4->TZ--5-|
                    /       |       |        /      ||    |
                    |       |->AH-2-|        |      BL    |
                    |                        |            |
                    |                        |<-----Kθ----|<-t5
                    |<------------------------------Kw----|<-w5


    Returns
    -------
    x       vector 16 elem.:
            θ0, w0, t1, w1, t2, w2, t3, w3, t4, w4, t5, w5,...
                QHC1, QHC2, QsTZ, QlTZ

    """
    Kθ, Kw = 1e10, 1e10             # controller gain
    wO = psy.w(θO, φO)            # hum. out
    wIsp = psy.w(θIsp, φIsp)      # indoor mumidity ratio

    # Model
    θs0, Δ_θs = θS, 2             # initial guess saturation temp.

    A = np.zeros((16, 16))          # coefficents of unknowns
    b = np.zeros(16)                # vector of inputs
    while Δ_θs > 0.01:
        # MX1
        A[0, 0], A[0, 10], b[0] = m * c, -(1 - α) * m * c, α * m * c * θO
        A[1, 1], A[1, 11], b[1] = m * l, -(1 - α) * m * l, α * m * l * wO
        # HC1
        A[2, 0], A[2, 2], A[2, 12], b[2] = m * c, -m * c, 1, 0
        A[3, 1], A[3, 3], b[3] = m * l, -m * l, 0
        # AH
        A[4, 2], A[4, 3], A[4, 4], A[4, 5], b[4] = c, l, -c, -l, 0
        A[5, 4], A[5, 5] = psy.wsp(θs0), -1
        b[5] = psy.wsp(θs0) * θs0 - psy.w(θs0, 1)
        # MX2
        A[6, 2], A[6, 4], A[6, 6], b[6] = β * m * c, (1 - β) * m * c, -m * c, 0
        A[7, 3], A[7, 5], A[7, 7], b[7] = β * m * l, (1 - β) * m * l, -m * l, 0
        # HC2
        A[8, 6], A[8, 8], A[8, 13], b[8] = m * c, -m * c, 1, 0
        A[9, 7], A[9, 9], b[9] = m * l, -m * l, 0
        # TZ
        A[10, 8], A[10, 10], A[10, 14], b[10] = m * c, -m * c, 1, 0
        A[11, 9], A[11, 11], A[11, 15], b[11] = m * l, -m * l, 1, 0
        # BL
        A[12, 10], A[12, 14], b[12] = (UA + mi * c), 1, (UA + mi * c
                                                         ) * θO + Qsa
        A[13, 11], A[13, 15], b[13] = mi * l, 1, mi * l * wO + Qla
        # Kθ & Kw
        A[14, 10], A[14, 12], b[14] = Kθ, 1, Kθ * θIsp
        A[15, 11], A[15, 13], b[15] = Kw, 1, Kw * wIsp

        x = np.linalg.solve(A, b)
        Δ_θs = abs(θs0 - x[4])
        θs0 = x[4]
    return x


def RecAirCAV(α=1, β=0.1,
              θS=30, θIsp=18, φIsp=0.49, θO=-1, φO=1,
              Qsa=0, Qla=0, mi=2.18, UA=935.83):
    """
    Model:
        Heating and adiabatic humidification
        Recycled air
        CAV Constant Air Volume:
            mass flow rate calculated for design conditions
            maintained constant in all situations

    INPUTS:
        α   mixing ratio of outdoor air, -
        β    by-pass factor of the adiabatic humidifier, -
        θS      supply air, °C
        θIsp    indoor air setpoint, °C
        φIsp  indoor relative humidity set point, -
        θO      outdoor temperature for design, °C
        φO    outdoor relative humidity for design, -
        Qsa     aux. sensible heat, W
        Qla     aux. latente heat, W
        mi      infiltration massflow rate, kg/s
        UA      global conductivity bldg, W/K

    System:
        MX1:    Mixing box
        HC1:    Heating Coil
        AH:     Adiabatic Humidifier
        MX2:    Mixing in humidifier model
        HC2:    Reheating coil
        TZ:     Thermal Zone
        BL:     Building
        Kw:     Controller - humidity
        Kθ:     Controller - temperature
        o:      outdoor conditions
        0..5    6 unknown points (temperature, humidity ratio)

        <----|<-------------------------------------------|
             |                                            |
             |              |-------|                     |
        -o->MX1--0->HC1--1->|       MX2--3->HC2--4->TZ--5-|
                    |       |       |        |      ||    |
                    |       |->AH-2-|        |      BL    |
                    |                        |            |
                    |                        |<-----Kθ----|<-t5
                    |<------------------------------Kw----|<-w5

    16 Unknowns
        0..5: 2*6 points (temperature, humidity ratio)
        QsHC1, QsHC2, QsTZ, QlTZ
    Returns
    -------
    None
    """
    plt.close('all')
    wO = psy.w(θO, φO)            # hum. out

    # Mass flow rate for design conditions
    # Supplay air mass flow rate
    # QsZ = UA*(θO - θIsp) + mi*c*(θO - θIsp)
    # m = - QsZ/(c*(θS - θIsp)
    # where
    # θO, wO = -1, 3.5e-3           # outdoor
    # θS = 30                       # supply air
    # mid = 2.18                     # infiltration
    QsZ = UA * (θOd - θIsp) + mid * c * (θOd - θIsp)
    m = - QsZ / (c * (θS - θIsp))
    print(f'm = {m: 5.3f} kg/s constant for design conditions:')
    print(f'    [θSd = {θS: 3.1f} °C, mi = 2.18 kg/S, θO = -1°C, φ0 = 100%]')

    # Model
    x = ModelRecAir(m, α, β,
                    θS, θIsp, φIsp, θO, φO, Qsa, Qla, mi, UA)

    θ = np.append(θO, x[0:12:2])
    w = np.append(wO, x[1:12:2])

    # Adjancy matrix
    # Points calc.  o   0   1   2   3   4   5       Elements
    # Points pplot  0   1   2   3   4   5   6       Elements
    A = np.array([[-1, +1, +0, +0, +0, +0, -1],     # MX1
                  [+0, -1, +1, +0, +0, +0, +0],     # HC1
                  [+0, +0, -1, +1, +0, +0, +0],     # AH
                  [+0, +0, -1, -1, +1, +0, +0],     # MX2
                  [+0, +0, +0, +0, -1, +1, +0],     # HC2
                  [+0, +0, +0, +0, +0, -1, +1]])    # TZ

    psy.chartA(θ, w, A)

    θ = pd.Series(θ)
    w = 1000 * pd.Series(w)
    P = pd.concat([θ, w], axis=1)       # points
    P.columns = ['θ [°C]', 'w [g/kg]']

    output = P.to_string(formatters={
        't [°C]': '{:,.2f}'.format,
        'w [g/kg]': '{:,.2f}'.format
    })
    print()
    print(output)

    Q = pd.Series(x[12:], index=['QsHC1', 'QsHC2', 'QsTZ', 'QlTZ'])
    # Q.columns = ['kW']
    pd.options.display.float_format = '{:,.2f}'.format
    print()
    print(Q.to_frame().T / 1000, 'kW')

    return None


def RecAirVAV(α=1, β=0.1,
              θSsp=30, θIsp=18, φIsp=0.49, θO=-1, φO=1,
              Qsa=0, Qla=0, mi=2.18, UA=935.83):
    """
    Created on Fri Apr 10 13:57:22 2020
    Heating & Adiabatic humidification & Re-heating
    Recirculated air
    VAV Variable Air Volume:
        mass flow rate calculated to have const. supply temp.

    INPUTS:
        α   mixing ratio of outdoor air, -
        β    by-pass factor of the adiabatic humidifier, -
        θS      supply air, °C
        θIsp    indoor air setpoint, °C
        φIsp  indoor relative humidity set point, -
        θO      outdoor temperature for design, °C
        φO    outdoor relative humidity for design, -
        Qsa     aux. sensible heat, W
        Qla     aux. latente heat, W
        mi      infiltration massflow rate, kg/s
        UA      global conductivity bldg, W/K

    System:
        MX1:    Mixing box
        HC1:    Heating Coil
        AH:     Adiabatic Humidifier
        MX2:    Mixing in humidifier model
        HC2:    Reheating coil
        TZ:     Thermal Zone
        BL:     Building
        Kw:     Controller - humidity
        Kθ:     Controller - temperature
        o:      outdoor conditions
        0..5    unknown points (temperature, humidity ratio)

        <----|<-------------------------------------------------|
             |                                                  |
             |              |-------|                           |
        -o->MX1--0->HC1--1->|       MX2--3->HC2--F-----4->TZ--5-|
                    |       |->AH-2-|        |   |     |  ||    |
                    |                        |   |-Kθ4-|  BL    |
                    |                        |                  |
                    |                        |<-----Kθ----------|<-t5
                    |<------------------------------Kw----------|<-w5

    16 Unknowns
        0..5: 2*6 points (temperature, humidity ratio)
        QsHC1, QsHC2, QsTZ, QlTZ
    """
    from scipy.optimize import least_squares

    def Saturation(m):
        """
        Used in VAV to find the mass flow which solves θS = θSsp
        Parameters
        ----------
            m : mass flow rate of dry air

        Returns
        -------
            θS - θSsp: difference between supply temp. and its set point

        """
        x = ModelRecAir(m, α, β,
                        θSsp, θIsp, φIsp, θO, φO, Qsa, Qla, mi, UA)
        θS = x[8]
        return (θS - θSsp)

    plt.close('all')
    wO = psy.w(θO, φO)            # hum. out

    # Mass flow rate
    res = least_squares(Saturation, 5, bounds=(0, 10))
    if res.cost < 1e-10:
        m = float(res.x)
    else:
        print('RecAirVAV: No solution for m')

    print(f'm = {m: 5.3f} kg/s')
    x = ModelRecAir(m, α, β,
                    θSsp, θIsp, φIsp, θO, φO, Qsa, Qla, mi, UA)

    # ΔθS, m = 2, 1                   # initial temp; diff; flow rate
    # while ΔθS > 0.01:
    #     m = m + 0.01
    #     # Model
    #     x = ModelRecAir(m, θSsp, mi, θO, φO, α, β)
    #     θS = x[8]
    #     ΔθS = -(θSsp - θS)

    θ = np.append(θO, x[0:12:2])
    w = np.append(wO, x[1:12:2])

    # Adjancy matrix
    # Points calc.  o   0   1   2   3   4   5       Elements
    # Points pplot  0   1   2   3   4   5   6       Elements
    A = np.array([[-1, +1, +0, +0, +0, +0, -1],     # MX1
                  [+0, -1, +1, +0, +0, +0, +0],     # HC1
                  [+0, +0, -1, +1, +0, +0, +0],     # AH
                  [+0, +0, -1, -1, +1, +0, +0],     # MX2
                  [+0, +0, +0, +0, -1, +1, +0],     # HC2
                  [+0, +0, +0, +0, +0, -1, +1]])    # TZ

    psy.chartA(θ, w, A)

    θ = pd.Series(θ)
    w = 1000 * pd.Series(w)
    P = pd.concat([θ, w], axis=1)       # points
    P.columns = ['θ [°C]', 'w [g/kg]']

    output = P.to_string(formatters={
        'θ [°C]': '{:,.2f}'.format,
        'w [g/kg]': '{:,.2f}'.format
    })
    print()
    print(output)

    Q = pd.Series(x[12:], index=['QsHC1', 'QsHC2', 'QsTZ', 'QlTZ'])
    # Q.columns = ['kW']
    pd.options.display.float_format = '{:,.2f}'.format
    print()
    print(Q.to_frame().T / 1000, 'kW')

    return None


# RecAirCAV()
# RecAirVAV()
