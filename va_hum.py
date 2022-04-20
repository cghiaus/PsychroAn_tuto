"""
Created on Tue Apr 21 11:36:00 2020
Updated on Sat Apr  2 17:33:24 2022

Heating and vapor humidification
Two HVAC systems: all out air and recirculated air
NOTE for CAV: specify the design conditions in this module

Building characterised by:
    UA      global conductivity
    mi      air infiltration mass flow
    Qsa     auxiliary sensible loads
    Qla     auxiliary latent loads

Indoor conditions
    θI      indoor temperature
    φI    indoor relative humidity
Outdoor conditions
    θO      outdoor temperature
    φO    outdoor relative humidity
Constant air volume (CAV)
    m       mass flow rate constant (cf. design conditions)
    θS      supply air temperature variable
Variable air volume (VAV)
    θS      supply air temperature controlled
    m       mass flow rate varable

Inputs given in Jupyter Notebook

"""
import numpy as np
import pandas as pd
import psychro as psy
import matplotlib.pyplot as plt

# Design conditions for CAV (to determine m)
θOd = -1        # outdoor temperarture
θId = 18        # indoor temperature
θSd = 30        # summply air temperature
mid = 0.0057    # infiltration
Qsad = 216      # 2163.5
Qlad = 145.2
mid = 0.0057
UAd = 85.9

# constants
c = 1e3         # air specific heat J/kg K
l = 2496e3      # latent heat J/kg

# *****************************************
# ALL OUT AIR
# *****************************************


def ModelAllOutAir(m, θS, θIsp, φIsp, θO, φO, Qsa, Qla, mi, UA):
    """
    Model:
        All outdoor air
        CAV Constant Air Volume:
            mass flow rate given
            control of indoor condition (t2, w2)
    INPUTS:
        m     mass flow of supply dry air kg/s
        θS    supply air °C
        θIsp  indoor air setpoint °C
        φIsp -
        θO    outdoor temperature for design °C
        φO    outdoor relative humidity for design -
        Qsa   aux. sensible heat W
        Qla   aux. latente heat W
        mi    infiltration massflow rate kg/s
        UA    global conductivity bldg W/K

    OUTPUTS:
        x     vector 10 elements:
            θ0, w0, θ1, w1, t2, w2, QsHC, QlVH, QsTZ, QlTZ

    System:
        HC:     Heating Coil
        VH:     Vapor Humidifier
        TZ:     Thermal Zone
        BL:     Building
        Kw:     Controller - humidity
        Kt:     Controller - temperature
        o:      outdoor conditions

    10 Unknowns
        0, 1, 2 points (temperature, humidity ratio)
        QsHC, QlVH, QsTZ, QlTZ
        x = [θ0, w0, θ1, w1, θ2, w2, QsHC, QlVH, QsTZ, QlTZ]

    --o->HC--0->VH--1->TZ--2-->
         s       l     sl  |
         |       |     ||  |
         |       |     BL  |
         |       |         |
         |       |<----Kw--|-w2
         |<------------Kt--|-t2
    """
    Kt, Kw = 1e10, 1e10             # controller gain
    wO = psy.w(θO, φO)              # outdoor mumidity ratio
    wIsp = psy.w(θIsp, φIsp)        # indoor mumidity ratio

    # Model
    A = np.zeros((10, 10))          # coefficents of unknowns
    b = np.zeros(10)                # vector of inputs
    # HC heating coil
    A[0, 0], A[0, 6], b[0] = m * c, -1, m * c * θO
    A[1, 1], b[1] = m * l, m * l * wO
    # VA vapor humidifier
    A[2, 0], A[2, 2], b[2] = -m * c, m * c, 0
    A[3, 1], A[3, 3], A[3, 7], b[3] = -m * l, m * l, -1, 0
    # TZ thermal zone
    A[4, 2], A[4, 4], A[4, 8], b[4] = -m * c, m * c, -1, 0
    A[5, 3], A[5, 5], A[5, 9], b[5] = -m * l, m * l, -1, 0
    # BL building
    A[6, 4], A[6, 8], b[6] = UA + mi * c, 1, (UA + mi * c) * θO + Qsa
    A[7, 5], A[7, 9], b[7] = mi * l, 1, mi * l * wO + Qla
    # Kt indoor temperature controller
    A[8, 4], A[8, 6], b[8] = Kt, 1, Kt * θIsp
    # Kw indoor hum.ratio controller
    A[9, 5], A[9, 7], b[9] = Kw, 1, Kw * wIsp

    # Solution
    x = np.linalg.solve(A, b)
    return x


def AllOutAirCAV(θS=30, θIsp=18, φIsp=0.5, θO=-1, φO=1,
                 Qsa=0, Qla=0, mi=2.12, UA=935.83):
    """
    All out air
    CAV Constant Air Volume:
        mass flow rate calculated for design conditions
        maintained constant in all situations

    INPUTS:
        m     mass flow of supply dry air kg/s
        θS    supply air °C
        θIsp  indoor air setpoint °C
        φIsp -
        θO    outdoor temperature for design °C
        φO    outdoor relative humidity for design -
        Qsa   aux. sensible heat W
        Qla   aux. latente heat W
        mi    infiltration massflow rate kg/s
        UA    global conductivity bldg W/K

    System:
        HC:     Heating Coil
        VH:     Vapor Humidifier
        TZ:     Thermal Zone
        BL:     Building
        Kw:     Controller - humidity
        Kt:     Controller - temperature
        o:      outdoor conditions

    10 Unknowns
        0, 1, 2 points (temperature, humidity ratio)
        QsHC, QlVH, QsTZ, QlTZ

    --o->HC--0->VH--1->TZ--2-->
         /       /     ||  |
         |       |     BL  |
         |       |         |
         |       |<----Kw--|-w2
         |<------------Kt--|-t2
    """
    plt.close('all')
    wO = psy.w(θO, φO)  # hum. out

    # Mass flow rate for design conditions
    # θOd = -1                        # outdoor design conditions
    # mid = 2.18                      # infiltration design
    QsZ = UAd * (θOd - θId) + mid * c * (θOd - θId) + Qsad
    m = - QsZ / (c * (θSd - θId))
    print('Winter All_out_air CAV')
    print(f'm = {m: 5.3f} kg/s constant (from design conditions)')
    print(f'Design conditions θS = {θS: 3.1f} °C,'
          f'mi = {mid:3.1f} kg/s, θO = {θOd:3.1f} °C, '
          f'θI = {θIsp:3.1f} °C')

    # Model
    x = ModelAllOutAir(m, θS, θIsp, φIsp, θO, φO, Qsa, Qla, mi, UA)

    # Processes on psychrometric chart
    t = np.append(θO, x[0:5:2])
    w = np.append(wO, x[1:6:2])

    # Adjancy matrix: rows=lines; columns=points
    # Points       O    0   1   2       Elements
    A = np.array([[-1, 1, 0, 0],     # HC
                 [0, -1, 1, 0],      # VH
                 [0, 0, 1, -1]])     # TZ

    psy.chartA(t, w, A)

    t = pd.Series(t)
    w = 1000 * pd.Series(w)
    P = pd.concat([t, w], axis=1)       # points
    P.columns = ['θ [°C]', 'w [g/kg]']
    P.index = ['O', 'HC', 'VH', 'I']

    output = P.to_string(formatters={
        'θ [°C]': '{:,.2f}'.format,
        'w [g/kg]': '{:,.2f}'.format
    })
    print()
    print(output)

    Q = pd.Series(x[6:], index=['QsHC', 'QlVH', 'QsTZ', 'QlTZ'])
    # Q.columns = ['kW']
    pd.options.display.float_format = '{:,.2f}'.format
    print()
    print(Q.to_frame().T / 1000, 'kW')

    return x


def AllOutAirVAV(θSsp=30, θIsp=18, φIsp=0.5, θO=-1, φO=1,
                 Qsa=0, Qla=0, mi=2.12, UA=935.83):
    """
    All out air
    Heating & Vapor humidification
    VAV Variable Air Volume:
        mass flow rate calculated to have const. supply temp.

    INPUTS:
        θS    supply air °C
        θIsp  indoor air setpoint °C
        φIsp -
        θO    outdoor temperature for design °C
        φO    outdoor relative humidity for design -
        Qsa   aux. sensible heat W
        Qla   aux. latente heat W
        mi    infiltration massflow rate kg/s
        UA    global conductivity bldg W/K

    System:
        HC:     Heating Coil
        VH:     Vapor Humidifier
        TZ:     Thermal Zone
        BL:     Building
        Kw:     Controller - humidity
        Kt:     Controller - temperature
        o:      outdoor conditions

    10 Unknowns
        0, 1, 2 points (temperature, humidity ratio)
        QsHC, QlVH, QsTZ, QlTZ

    --o->HC--0->VH--F-----1-->TZ--2-->
         /       /  m     θ   //  θ,w
         s       l  |     |   sl  |
         |       |  |_Kθ1_|   BL  |
         |       |                |
         |       |<----Kw---------|-w2
         |<------------Kt---------|-θ2

        Mass-flow rate (VAV) I-controller:
        start with m = 0
        measure the supply temperature
        while -(θSsp - θS)>0.01, increase m (I-controller)
    """
    plt.close('all')
    wO = psy.w(θO, φO)            # outdoor mumidity ratio

    # Mass flow rate
    DθS, m = 2, 0                   # initial temp; diff; flow rate
    while DθS > 0.01:
        m = m + 0.01                # mass-flow rate I-controller
        # Model
        x = ModelAllOutAir(m, θSsp, θIsp, φIsp, θO, φO, Qsa, Qla, mi, UA)
        θS = x[2]
        DθS = -(θSsp - θS)
    print('Winter All_out_air VAV')
    print(f'm = {m: 5.3f} kg/s')
    # Processes on psychrometric chart
    t = np.append(θO, x[0:5:2])
    w = np.append(wO, x[1:6:2])
    # Points       o    0   1   2   Elements
    A = np.array([[-1, 1, 0, 0],    # HC
                 [0, -1, 1, 0],     # VH
                 [0, 0, 1, -1]])    # TZ

    psy.chartA(t, w, A)

    t = pd.Series(t)
    w = 1000 * pd.Series(w)
    P = pd.concat([t, w], axis=1)       # points
    P.columns = ['θ [°C]', 'w [g/kg]']

    output = P.to_string(formatters={
        'θ [°C]': '{:,.2f}'.format,
        'w [g/kg]': '{:,.2f}'.format
    })
    print()
    print(output)

    Q = pd.Series(x[6:], index=['QsHC', 'QlVH', 'QsTZ', 'QlTZ'])
    # Q.columns = ['kW']
    pd.options.display.float_format = '{:,.2f}'.format
    print()
    print(Q.to_frame().T / 1000, 'kW')

    return None


# *****************************************
# RECYCLED AIR
# *****************************************


def ModelRecAir(m, α, θS, θIsp, φIsp, θO, φO, Qsa, Qla, mi, UA):
    """
    Model:
        Heating and vapor humidification
        Recycled air
        CAV Constant Air Volume:
            mass flow rate calculated for design conditions
            maintained constant in all situations
    INPUTS:
        m     mass flow of supply dry air kg/s
        α mixing ratio of outdoor air
        θS    supply air °C
        θIsp  indoor air setpoint °C
        φIsp -
        θO    outdoor temperature for design °C
        φO  outdoor relative humidity for design -
        Qsa   aux. sensible heat W
        Qla   aux. latente heat W
        mi    infiltration massflow rate kg/s
        UA    global conductivity bldg W/K

    OUTPUTS:
        x     vector 12 elements:
            θ0, w0, θ1, w1, t2, w2, t3, w3, QsHC, QlVH, QsTZ, QlTZ

    System:
        MX:     Mixing Box
        HC:     Heating Coil
        VH:     Vapor Humidifier
        TZ:     Thermal Zone
        BL:     Buildings
        Kw:     Controller - humidity
        Kt:     Controller - temperature
        o:      outdoor conditions

    12 Unknowns
        0, 1, 2, 3 points (temperature, humidity ratio)
        QsHC, QlVH, QsTZ, QlTZ

    <-3--|<-------------------------|
         |                          |
    -o->MX--0->HC--1->VH--2->TZ--3-->
               /       /     ||  |
               |       |     BL  |
               |       |         |
               |       |<----Kw--|-w3
               |<------------Kt--|-t3
    """
    Kt, Kw = 1e10, 1e10             # controller gain
    wO = psy.w(θO, φO)            # hum. out
    wIsp = psy.w(θIsp, φIsp)      # hum. in set point

    # Model
    A = np.zeros((12, 12))          # coefficents of unknowns
    b = np.zeros(12)                # vector of inputs
    # MX mixing box
    A[0, 0], A[0, 6], b[0] = m * c, -(1 - α) * m * c, α * m * c * θO
    A[1, 1], A[1, 7], b[1] = m * l, -(1 - α) * m * l, α * m * l * wO
    # HC hearing coil
    A[2, 0], A[2, 2], A[2, 8], b[2] = m * c, -m * c, 1, 0
    A[3, 1], A[3, 3], b[3] = m * l, -m * l, 0
    # VH vapor humidifier
    A[4, 2], A[4, 4], b[4] = m * c, -m * c, 0
    A[5, 3], A[5, 5], A[5, 9], b[5] = m * l, -m * l, 1, 0
    # TZ thermal zone
    A[6, 4], A[6, 6], A[6, 10], b[6] = m * c, -m * c, 1, 0
    A[7, 5], A[7, 7], A[7, 11], b[7] = m * l, -m * l, 1, 0
    # BL building
    A[8, 6], A[8, 10], b[8] = (UA + mi * c), 1, (UA + mi * c) * θO + Qsa
    A[9, 7], A[9, 11], b[9] = mi * l, 1, mi * l * wO + Qla
    # Kt indoor temperature controller
    A[10, 6], A[10, 8], b[10] = Kt, 1, Kt * θIsp
    # Kw indoor humidity controller
    A[11, 7], A[11, 9], b[11] = Kw, 1, Kw * wIsp

    # Solution
    x = np.linalg.solve(A, b)
    return x


def RecAirCAV(α=0.5, θS=30, θIsp=18, φIsp=0.5, θO=-1, φO=1,
              Qsa=0, Qla=0, mi=2.12, UA=935.83):
    """
    CAV Constant Air Volume:
    mass flow rate calculated for design conditions
    maintained constant in all situations
    INPUTS:
        α mixing ratio of outdoor air
        θS    supply air °C
        θIsp  indoor air setpoint °C
        φIsp -
        θO    outdoor temperature for design °C
        φO  outdoor relative humidity for design -
        Qsa   aux. sensible heat W
        Qla   aux. latente heat W
        mi    infiltration massflow rate kg/s
        UA    global conductivity bldg W/K

    System:
        HC:     Heating Coil
        VH:     Vapor Humidifier
        TZ:     Thermal Zone
        Kw:     Controller - humidity
        Kt:     Controller - temperature
        o:      outdoor conditions

    12 Unknowns
        0, 1, 2, 3 points (temperature, humidity ratio)
        QsHC, QlVH, QsTZ, QlTZ

    <-3--|<-------------------------|
         |                          |
    -o->MX--0->HC--1->VH--2->TZ--3-->
               /       /     ||  |
               |       |     BL  |
               |       |         |
               |       |_____Kw__|_w3
               |_____________Kt__|_t3
    """
    plt.close('all')
    wO = psy.w(θO, φO)            # hum. out

    # Mass flow rate for design conditions
    # Supplay air mass flow rate
    # QsZ = UA*(θO - θIsp) + mi*c*(θO - θIsp)
    # m = - QsZ/(c*(θS - θIsp)
    # where
    # θOd, wOd = -1, 3.5e-3           # outdoor
    # θS = 30                       # supply air
    # mid = 2.18                     # infiltration
    QsZ = UA * (θOd - θIsp) + mid * c * (-1 - θIsp) + Qsa
    m = - QsZ / (c * (θS - θIsp))
    print('Winter Recirculated_air CAV')
    print(f'm = {m: 5.3f} kg/s constant (from design conditions)')
    print(f'Design conditions θS = {θS: 3.1f} °C,'
          f'mi = {mid:3.1f} kg/s, θO = {θOd:3.1f} °C, '
          f'θI = {θIsp:3.1f} °C')

    # Model
    x = ModelRecAir(m, α, θS, θIsp, φIsp, θO, φO, Qsa, Qla, mi, UA)
    # (m, θS, mi, θO, φO, α)

    # Processes on psychrometric chart
    # Points      o    0    1   2   3       Elements
    #             0    1    2   3   4
    A = np.array([[-1, 1, 0, 0, -1],        # MX
                 [0, -1, 1, 0, 0],          # HC
                 [0, 0, -1, 1, 0],          # VH
                 [0, 0, 0, -1, 1]])         # TZ
    t = np.append(θO, x[0:8:2])

    print(f'wO = {wO:6.5f}')
    w = np.append(wO, x[1:8:2])
    psy.chartA(t, w, A)

    t = pd.Series(t)
    w = 1000 * pd.Series(w)
    P = pd.concat([t, w], axis=1)       # points
    P.columns = ['θ [°C]', 'w [g/kg]']

    output = P.to_string(formatters={
        'θ [°C]': '{:,.2f}'.format,
        'w [g/kg]': '{:,.2f}'.format
    })
    print()
    print(output)

    Q = pd.Series(x[8:], index=['QsHC', 'QlVH', 'QsTZ', 'QlTZ'])
    pd.options.display.float_format = '{:,.2f}'.format
    print()
    print(Q.to_frame().T / 1000, 'kW')

    return None


def RecAirVAV(α=0.5, θSsp=30, θIsp=18, φIsp=0.5, θO=-1, φO=1,
              Qsa=0, Qla=0, mi=2.12, UA=935.83):
    """
    CAV Variable Air Volume:
    mass flow rate calculated s.t.
    he supply temp. is maintained constant in all situations
    INPUTS:
    INPUTS:
        m     mass flow of supply dry air kg/s
        α mixing ratio of outdoor air
        θS    supply air °C
        θIsp  indoor air setpoint °C
        φIsp -
        θO    outdoor temperature for design °C
        φO  outdoor relative humidity for design -
        Qsa   aux. sensible heat W
        Qla   aux. latente heat W
        mi    infiltration massflow rate kg/s
        UA    global conductivity bldg W/K

    System (CAV & m introduced by the Fan is cotrolled by θS )
        HC:     Heating Coil
        VH:     Vapor Humidifier
        TZ:     Thermal Zone
        F:      Supply air fan
        Kw:     Controller - humidity
        Kt:     Controller - temperature
        o:      outdoor conditions

    12 Unknowns
        0, 1, 2, 3 points (temperature, humidity ratio)
        QsHC, QlVH, QsTZ, QlTZ

    <----|<--------------------------------|
         |                                 |
    -o->MX--0->HC--1->VH--F-----2-->TZ--3-->
               /       /  |     |   ||  |
               |       |  |     |   BL  |
               |       |  |     |       |
               |       |  |_Kt2_|_t2    |
               |       |                |
               |       |_____Kw_________|_w3
               |_____________Kt_________|_t3

    Mass-flow rate (VAV) I-controller:
        start with m = 0
        measure the supply temperature
        while -(θSsp - θS)>0.01, increase m (I-controller)
    """
    plt.close('all')
    wO = psy.w(θO, φO)            # hum. out

    # Mass flow rate
    DθS, m = 2, 0                   # initial temp; diff; flow rate
    while DθS > 0.01:
        m = m + 0.01

        # Model
        x = ModelRecAir(m, α, θSsp, θIsp, φIsp, θO, φO,
                        Qsa, Qla, mi, UA)
        θS = x[4]
        DθS = -(θSsp - θS)

    print('Winter Rec_air VAV')
    print(f'm = {m: 5.3f} kg/s')

    # Processes on psychrometric chart
    # Points      o    0    1   2   3       Elements
    #             0    1    2   3   4
    A = np.array([[-1, 1, 0, 0, -1],        # MX
                 [0, -1, 1, 0, 0],          # HC
                 [0, 0, -1, 1, 0],          # VH
                 [0, 0, 0, -1, 1]])         # TZ
    t = np.append(θO, x[0:8:2])
    print(f'wO = {wO:6.5f}')
    w = np.append(wO, x[1:8:2])
    psy.chartA(t, w, A)

    t = pd.Series(t)
    w = 1000 * pd.Series(w)
    P = pd.concat([t, w], axis=1)           # points
    P.columns = ['θ [°C]', 'w [g/kg]']

    output = P.to_string(formatters={
        'θ [°C]': '{:,.2f}'.format,
        'w [g/kg]': '{:,.2f}'.format
    })
    print()
    print(output)

    Q = pd.Series(x[8:], index=['QsHC', 'QlVH', 'QsTZ', 'QlTZ'])
    pd.options.display.float_format = '{:,.2f}'.format
    print()
    print(Q.to_frame().T / 1000, 'kW')

    return None


# Uncomment to test a function
# AllOutAirCAV(θS=30, θIsp=18, φIsp=0.5, θO=-1, φO=1,
#              Qsa=2163.5, Qla=145.2, mi=0.0057, UA=85.9)
# AllOutAirCAV()
# AllOutAirVAV()
# RecAirCAV()
# RecAirVAV()
