#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 07:14:57 2021

@author: cghiaus
from cool10b.py
mo  outdoor mass flow rate as an input.

Checks as a direct problem Ex. 7.6
1. tout air neuf (pg. 2 & 9-10)
2. by-pass (pg. 19)
3. recycl & by-pass (pg. 19)

Cooling as a control & parameter optilizaton problem OOP
Recycling, Cooling & desumidification (with by-pass), reheating

GENERALITIES
========================================================================
Units
Temperatures: °C
Humidity ration: kg_vapor / kg_dry_air
Relative humidity: -
Heat flow rate: W
Mass flow rate: kg/s

Points on psychrometric chart (θ, w):
o) out      outdoor
0) M        mixed: fresh + recycled
1) s        efective coil surface temperature (ADP: Apparatus Dew Point)
2) C        coil leaving air temp. LAT (saturated and by-passed)
3) S        supply
4) I        indoor


System as a direct problem
--------------------------
<=4================================4============================
  mo     ||                        m                          ||
         4 (m-mo) =======0=======                             ||
         ||       ||  (1-β)m   ||                             ||
θo,φo==>[MX1]==0==||          [MX2]==2==[HC]==F==3==>[TZ]==4==||
 mo               ||           ||        /   /       //       |
                  ===0=[CC]==1===       s   m       sl        |
                       /\\   βm         |           ||        |
                      t  sl             |          [BL]<-mi   |
                      |                 |          //         |
                      |                 |         sl          |
                      |                 |                     |
                      |                 |<------[K]-----------|<-wI
                      |<------------------------[K]-----------|<-θI

Inputs:
θo, φo      outdoor temperature & humidity ratio
θIsp, φIsp  indoor temperature & humidity ratio set points
mi          infiltration mass flow rate
Qsa, Qla    auxiliary sensible and latent loads [kW]
Parameters:
m           mass flow rate of dry air
α           fraction of fresh air
β           by-pass factir od cooling coil
UA          overall heat transfer coefficient

Elements (16 equations):
MX1         mixing box (2 equations)
CC          cooling coil (4 equations)
MX2         mixing process (2 equations)
HC          heating coil (2 equations)
TZ          thermal zone (2 equations)
BL          building (2 equations)
Kθ          indoor temperature controller (1 equation)
Kw          indoor humidity controller (1 equation)
F           fan (m is a given parameter)

Outputs (16 unknowns):
0, ..., 4   temperature and humidity ratio (10 unknowns)
Qt, Qs, Ql  total, sensible and latent heat of CC (3 unknowns)
Qs          sensible heat of HC (1 unknown)
Qs, Ql      sensible and latent heat of TZ (2 unknowns)






CAV System with linear controllers for θI & φI
----------------------------------------------
 out        s              S          I
==0==>[CC]==1==>[HC]===F===2===>[TZ]==3==>
       /\\      /     /         //    ||
      t  sl    s     m         sl     ||
      |        |                      ||
      |        |                      ||
      |        |<------[K]------------||<-wI<-φI
      |<---------------[K]------------|<-θI

Inputs:
θo, φo      outdoor temperature & relative humidity
θI, φI      indoor air temperature and humidity
QsTZ        sensible heat load of TZ
QlTZ        latent heat load of TZ
Parameter:
m           mass flow rate of dry air

Elements (10 equations):
CC          cooling coil (4 equations)
HC          heating coil (2 equations)
TZ          thermal zone (2 equations)
F           fan (no equation, m is given)
KθI         indoor temperature controller (1 equation)
KwI         indoor humidity controller (1 equation)

Outputs (10 unknowns):
0, 1, 2     temperature and humidity ratio (6 unknowns)
QsCC, QlCC  sensible and latent heat of CC (2 unknowns)
QtCC        total heat load of CC (1 unknown)
QsHC        sensible heat load of HC (1 unknown)


VAV System with linear & least-squares controllers for  & θS
------------------------------------------------------------------

linear controller (Kθ & Kw) for θI, φI
non-linear controller (ls) for θS

<=4================================m==========================
       ||                                                   ||
       4 (m-mo) =======0=======                             ||
out    ||    M  ||  (1-β)m   ||    C            S        I  ||
mo===>[MX1]==0==||          [MX2]==2==[HC]==F===3=>[TZ]==4==||
                ||         s ||        /   /    |    //     |
                ===0=[CC]==1===       s   m     |   sl      |
                     /\\   βm         |   |     |   ||      |
                    t  sl             |   |<-ls-|  [BL]<-mi |
                    |                 |            //       |
                    |                 |           sl        |
                    |                 |                     |
                    |                 |<------[K]-----------|<-wI
                    |<------------------------[K]-----------|<-θI

Inputs:
θo, φo      outdoor temperature & relative humidity
θI, φI      indoor air temperature and humidity
θS          supply air temperature
QsTZ        sensible heat load of TZ
QlTZ        latent heat load of TZ

Elements (11 equations):
CC          cooling coil (4 equations)
HC          heating coil (2 equations)
TZ          thermal zone (2 equations)
F           fan (m is given)
KθI         indoor temperature controller (1 equation)
KwI         indoor humidity controller (1 equation)
lsθS        mass flow rate of dry air controller (1 non-linear equation)

Outputs (11 unknowns):
0, 1, 2     temperature and humidity ratio (6 unknowns)
QsCC, QlCC  sensible and latent heat of CC (2 unknowns)
QtCC        total heat load of CC (1 unknown)
QsHC        sensible heat load of HC (1 unknown)
Parameter:
m           mass flow rate of dry air (1 unknown)


CONTENTS (methods)
========================================================================
__init__    Initialization of CcTZ object.
lin_model   Solves the set of linear equations
            with saturation curve linearized around ts0
solve_lin   Solves iteratively the lin_model s.t. the error of
            humid. ratio between two iterrations is approx. zero
            (i.e. solves ws = f(θs) for saturation curve).
m_ls        Finds m s.t. θS = θSsp (solves θS - θSsp = 0 for m).
            Uses least-squares to find m that minimizes θS - θSsp
psy_chart   Draws psychrometric chart (imported from psychro)
CAV_wd      CAV to be used in Jupyter widgets.
            solve_lin and draws psy_chart.
VAV_wd      VAV to be used in Jupyter widgets.
            m_ls and draws psy_chart.
"""
import numpy as np
import pandas as pd
import psychro as psy

# constants
c = 1e3         # J/kg K, air specific heat
l = 2496e3      # J/kg, latent heat

# to be used in self.m_ls / least_squares
m_max = 100     # ks/s, max dry air mass flow rate
θs_0 = 5        # °C, initial guess for saturation temperature


class MxCcRhTzBl:
    """
    **HVAC composition**:
        mixing, cooling, reaheating, thermal zone of building, recycling
    """

    def __init__(self, parameters, inputs):
        m, mo, β, Kθ, Kw = parameters
        θo, φo, θIsp, φIsp, mi, UA, Qsa, Qla = inputs

        self.design = np.array([m, mo, β, Kθ, Kw,       # parameters
                                θo, φo, θIsp, φIsp,     # inputs air out, in
                                mi, UA, Qsa, Qla])      # --"--  building
        self.actual = np.array([m, mo, β, Kθ, Kw,
                                θo, φo, θIsp, φIsp,
                                mi, UA, Qsa, Qla])

    def lin_model(self, θs0):
        """
        Linearized model.
            Solves a set of 16 linear equations.
            Saturation curve is linearized in θs0.

        s-point (θs, ws):

        - is on a tangent to φ = 100 % in θs0;

        - is **not** on the saturation curve (Apparatus Dew Point ADP).


        Parameter from function call
        ----------------------------
        θs0     °C, temperature for which the saturation curve is liniarized

        Parameters from object
        ---------------------
        m, mo, θo, φo, θIsp, φIsp, β, mi, UA, Qsa, Qla = self.actual

        Equations (16)
        -------------
        +-------------+-----+----+-----+----+----+----+----+----+
        | Element     | MX1 | CC | MX2 | HC | TZ | BL | Kθ | Kw |
        +=============+=====+====+=====+====+====+====+====+====+
        | N° equations|  2  | 4  |  2  |  2 | 2  | 2  |  1 |  1 |
        +-------------+-----+----+-----+----+----+----+----+----+

        Returns (16 unknowns)
        ---------------------
        x : θM, wM, θs, ws, θC, wC, θS, wS, θI, wI,
            QtCC, QsCC, QlCC, QsHC, QsTZ, QlTZ
        """
        """
        <=4================================m==========================
               ||                                                   ||
               4 (m-mo) =======0=======                             ||
               ||       ||  (1-β)m   ||                             ||
       θo,φo=>[MX1]==0==||          [MX2]==2==[HC]==F==3==>[TZ]==4==||
         mo             ||           ||        /   /       //       |
                        ===0=[CC]==1===       s   m       sl        |
                             /\\   βm         |           ||        |
                            t  sl             |          [BL]<-mi   |
                            |                 |          //         |
                            |                 |         sl          |
                            |                 |                     |
                            |                 |<------[K]-----------+<-wI
                            |<------------------------[K]-----------+<-θI
        """
        m, mo, β, Kθ, Kw, θo, φo, θIsp, φIsp, mi, UA, Qsa, Qla = self.actual
        wo = psy.w(θo, φo)      # hum. out

        A = np.zeros((16, 16))  # coefficents of unknowns
        b = np.zeros(16)        # vector of inputs
        # MX1
        A[0, 0], A[0, 8], b[0] = m * c, -(m - mo) * c, mo * c * θo
        A[1, 1], A[1, 9], b[1] = m * l, -(m - mo) * l, mo * l * wo
        # CC
        A[2, 0], A[2, 2], A[2, 11], b[2] = (1 - β) * m * c, -(1 - β) * m * c,\
            1, 0
        A[3, 1], A[3, 3], A[3, 12], b[3] = (1 - β) * m * l, -(1 - β) * m * l,\
            1, 0
        A[4, 2], A[4, 3], b[4] = psy.wsp(θs0), -1,\
            psy.wsp(θs0) * θs0 - psy.w(θs0, 1)
        A[5, 10], A[5, 11], A[5, 12], b[5] = -1, 1, 1, 0
        # MX2
        A[6, 0], A[6, 2], A[6, 4], b[6] = β * m * c, (1 - β) * m * c,\
            - m * c, 0
        A[7, 1], A[7, 3], A[7, 5], b[7] = β * m * l, (1 - β) * m * l,\
            - m * l, 0
        # HC
        A[8, 4], A[8, 6], A[8, 13], b[8] = m * c, -m * c, 1, 0
        A[9, 5], A[9, 7], b[9] = m * l, -m * l, 0
        # TZ
        A[10, 6], A[10, 8], A[10, 14], b[10] = m * c, -m * c, 1, 0
        A[11, 7], A[11, 9], A[11, 15], b[11] = m * l, -m * l, 1, 0
        # BL
        A[12, 8], A[12, 14], b[12] = (UA + mi * c), 1, (UA + mi * c) * θo + Qsa
        A[13, 9], A[13, 15], b[13] = mi * l, 1, mi * l * wo + Qla
        # Kθ indoor temperature controller
        A[14, 8], A[14, 10], b[14] = Kθ, 1, Kθ * θIsp
        # Kw indoor humidity ratio controller
        A[15, 9], A[15, 13], b[15] = Kw, 1, Kw * psy.w(θIsp, φIsp)
        x = np.linalg.solve(A, b)
        return x

    def solve_lin(self, θs0):
        """
        Finds saturation point on saturation curve ws = f(θs).
            Solves iterativelly *lin_model(θs0)*:
            θs -> θs0 until ws = psy(θs, 1).

        Parameters
        ----------
        θs0     initial guess saturation temperature

        Method from object
        ---------------------
        *self.lin_model(θs0)*

        Returns (16 unknowns)
        ---------------------
        x of *self.lin_model(self, θs0)*
        """
        Δ_ws = 10e-3  # kg/kg, initial difference to start the iterations
        while Δ_ws > 0.01e-3:
            x = self.lin_model(θs0)
            Δ_ws = abs(psy.w(x[2], 1) - x[3])   # psy.w(θs, 1) = ws
            θs0 = x[2]                          # actualize θs0
        return x

    def m_ls(self, value, sp):
        """
        Mass flow rate m controls supply temperature θS or indoor humidity wI.
            Finds m which solves value = sp, i.e. minimizes ε = value - sp.
            Uses *scipy.optimize.least_squares* to solve the non-linear system.

        Parameters
        ----------
        value   string: 'θS' od 'wI' type of controlled variable
        sp      float: value of setpoint

        Calls
        -----
        *ε(m)*  gives (value - sp) to be minimized for m

        Returns (16 unknowns)
        ---------------------
        x           given by *self.lin_model(self, θs0)*
        """
        from scipy.optimize import least_squares

        def ε(m):
            """
            Gives difference ε = (values - sp) function of m
                ε  calculated by self.solve_lin(ts0)
                m   bounds=(0, m_max); m_max hard coded (global variable)

            Parameters
            ----------
            m : mass flow rate of dry air

            From object
                Method: self.solve.lin(θs0)
                Variables: self.actual <- m (used in self.solve.lin)
            Returns
            -------
            ε = value - sp: difference between value and its set point
            """
            self.actual[0] = m
            x = self.solve_lin(θs_0)
            if value == 'θS':
                θS = x[6]       # supply air
                return abs(sp - θS)
            elif value == 'φI':
                wI = x[9]       # indoor air
                return abs(sp - wI)
            else:
                print('ERROR in ε(m): value not in {"θS", "wI"}')

        m0 = self.actual[0]     # initial guess
        if value == 'φI':
            self.actual[4] = 0
            sp = psy.w(self.actual[7], sp)
        # gives m for min(θSsp - θS); θs_0 is the initial guess of θs
        res = least_squares(ε, m0, bounds=(0, m_max))

        if res.cost < 0.1e-3:
            m = float(res.x)
            # print(f'm = {m: 5.3f} kg/s')
        else:
            print('RecAirVAV: No solution for m')

        self.actual[0] = m

        x = self.solve_lin(θs_0)
        return x

    def β_ls(self, value, sp):
        """
        Bypass β controls supply temperature θS or indoor humidity wI.
            Finds β which solves value = sp, i.e. minimizes ε = value - sp.
            Uses *scipy.optimize.least_squares* to solve the non-linear system.

        Parameters
        ----------
        value   string: 'θS' od 'wI' type of controlled variable
        sp      float: value of setpoint

        Calls
        -----
        *ε(m)*  gives (value - sp) to be minimized for m

        Returns (16 unknowns)
        ---------------------
        x           given by *self.lin_model(self, θs0)*
        """
        from scipy.optimize import least_squares

        def ε(β):
            """
            Gives difference ε = (values - sp) function of β
                ε  calculated by self.solve_lin(ts0)
                β   bounds=(0, 1)

            Parameters
            ----------
            β : by-pass factor of the cooling coil

            From object
                Method: self.solve.lin(θs0)
                Variables: self.actual <- m (used in self.solve.lin)
            Returns
            -------
            ε = value - sp: difference between value and its set point
            """
            self.actual[2] = β
            x = self.solve_lin(θs_0)
            if value == 'θS':
                θS = x[6]       # supply air
                return abs(sp - θS)
            elif value == 'φI':
                wI = x[9]       # indoor air
                return abs(sp - wI)
            else:
                print('ERROR in ε(β): value not in {"θS", "wI"}')

        β0 = self.actual[2]     # initial guess
        β0 = 0.1
        if value == 'φI':
            self.actual[4] = 0
            sp = psy.w(self.actual[7], sp)
        # gives m for min(θSsp - θS); θs_0 is the initial guess of θs
        res = least_squares(ε, β0, bounds=(0, 1))

        if res.cost < 1e-5:
            β = float(res.x)
            # print(f'm = {m: 5.3f} kg/s')
        else:
            print('RecAirVBP: No solution for β')

        self.actual[2] = β
        x = self.solve_lin(θs_0)
        return x

    def psy_chart(self, x, θo, φo):
        """
        Plot results on psychrometric chart.

        Parameters
        ----------
        x : θM, wM, θs, ws, θC, wC, θS, wS, θI, wI,
            QtCC, QsCC, QlCC, QsHC, QsTZ, QlTZ
                    results of self.solve_lin or self.m_ls
        θo, φo      outdoor point

        Returns
        -------
        None.

        """
        # Processes on psychrometric chart
        wo = psy.w(θo, φo)
        # Points: O, s, S, I
        θ = np.append(θo, x[0:10:2])
        w = np.append(wo, x[1:10:2])
        # Points       0   1  2  3  4  5       Elements
        A = np.array([[-1, 1, 0, 0, 0, 1],      # MR
                      [0, -1, 1, 0, 0, 0],      # CC
                      [0, 0, -1, 1, -1, 0],     # MX
                      [0, 0, 0, -1, 1, 0],      # HC
                      [0, 0, 0, 0, -1, 1]])     # TZ
        psy.chartA(θ, w, A)

        θ = pd.Series(θ)
        w = 1000 * pd.Series(w)         # kg/kg -> g/kg
        P = pd.concat([θ, w], axis=1)   # points
        P.columns = ['θ [°C]', 'w [g/kg]']

        output = P.to_string(formatters={
            'θ [°C]': '{:,.2f}'.format,
            'w [g/kg]': '{:,.2f}'.format})
        print()
        print(output)

        Q = pd.Series(x[10:], index=['QtCC', 'QsCC', 'QlCC', 'QsHC',
                                     'QsTZ', 'QlTZ'])
        # Q.columns = ['kW']
        pd.options.display.float_format = '{:,.2f}'.format
        print()
        print(Q.to_frame().T / 1000, 'kW')
        return None

    def CAV_wd(self, θo=32, φo=0.80, θIsp=26, φIsp=0.5,
               mi=1.35, UA=675, QsBL=34_000, QlBL=4_000):
        """
        Constant air volume (CAV) to be used in Jupyter with widgets

        Parameters: given in Jupyetr widget
        ----------

        Returns
        -------
        None.
        """
        # To use fewer variables in Jupyter widget:
        # select what to be updated in self.actual, e.g.:
        # self.actual[[0, 1, 2, 5, 6]] = m, θo, φo, 1000 * QsTZ, 1000 * QlTZ
        self.actual[5:] = np.array([θo, φo, θIsp, φIsp,
                                    mi, UA, QsBL, QlBL])
        # self.actual[5:] = θo, φo, θIsp, φIsp, mi, UA, Qsa, Qla

        θ0 = 40
        x = self.solve_lin(θ0)
        # print(f'm = {self.actual[0]: .3f} kg/s,\
        #       mo = {self.actual[1]: .3f} kg/s')
        print('m = {m: .3f} kg/s, mo = {mo: .3f} kg/s'.format(
            m=self.actual[0], mo=self.actual[1]))
        self.psy_chart(x, self.actual[5], self.actual[6])

    def VAV_wd(self, value='θS', sp=18, θo=32, φo=0.5, θIsp=24, φIsp=0.5,
               mi=1.35, UA=675, QsBL=34_000, QlBL=4_000):
        """
        Variable air volume (VAV) to be used in Jupyter with widgets

        Parameters
        ----------
        value       {"θS", "wI"}' type of value controlled
        sp          set point for the controlled value
        θo, φo, θIsp, φIsp, mi, UA, Qsa, Qla
                    given by widgets in Jupyter Lab

        Returns
        -------
        None.
        """
        """
        value='θS' (KwI = 0)

        <=4================================4===========================
                ||                         m                          ||
                4 (m-mo) =======0=======                              ||
                ||    M  ||  (1-β)m   ||    C            S         I  ||
        θo,φo=>[MX1]==0==||          [MX2]==2==[HC]==F===3==>[TZ]==4==||
         mo              ||         s ||        /   /    |    //      |
                         ===0=[CC]==1===       s   m     |   sl       |
                              /\\   βm         |   |     |   ||       |
                             t  sl             |   |     |  [BL]<-mi  |
                             |                 |   |     |   //       |
                             |                 |   |     |  sl        |
                             |                 |   |     |            |
                             |                 |   |--ls-|<-θS        |
                             |                 |<-----[K]-------------|<-wI
                             |<-----------------------[K]-------------|<-θI

        value='wI' (KwI = 0)

        <=4================================4===========================
                ||                         m                          ||
                4 (m-mo) =======0=======                              ||
                ||    M  ||  (1-β)m   ||    C            S         I  ||
        θo,φo=>[MX1]==0==||          [MX2]==2==[HC]==F===3==>[TZ]==4==||
         mo              ||         s ||        /   /         //      |
                         ===0=[CC]==1===       s   m         sl       |
                              /\\   βm         |   |         ||       |
                             t  sl             |   |        [BL]<-mi  |
                             |                 |   |         //       |
                             |                 |   |        sl        |
                             |                 |   |                  |
                             |                 |   |--ls--------------|<-wI
                             |                 |<-----[K]-------------|<-wI
                             |<-----------------------[K]-------------|<-θI

        """
        # Design values
        self.actual[5:] = θo, φo, θIsp, φIsp, mi, UA, QsBL, QlBL

        x = self.m_ls(value, sp)
        print('m = {m: .3f} kg/s, mo = {mo: .3f} kg/s'.format(
            m=self.actual[0], mo=self.actual[1]))
        self.psy_chart(x, θo, φo)

    def VBP_wd(self, value='θS', sp=18, θo=32, φo=0.5, θIsp=24, φIsp=0.5,
               mi=1.35, UA=675, Qsa=34_000, Qla=4_000):
        """
        Variabl by-pass (VBP) to be used in Jupyter with widgets

        Parameters
        ----------
        value       {"θS", "wI"}' type of value controlled
        sp          set point for the controlled value
        θo, φo, θIsp, φIsp, mi, UA, Qsa, Qla
                    given by widgets in Jupyter Lab

        Returns
        -------
        None.
        """
        """
        value='θS'

        <=4================================4============================
                ||                         m                          ||
                4 (m-mo) =======0=======                              ||
                ||    M  ||  (1-β)m   ||    C            S         I  ||
        θo,φo=>[MX1]==0==||          [MX2]==2==[HC]==F===3==>[TZ]==4==||
         mo              ||         s ||        /   /    |    //      |
                         ===0=[CC]==1===       s   m     |   sl       |
                              /\\   βm         |         |   ||       |
                             t  sl  |          |         |  [BL]<-mi  |
                             |      |          |         |  //        |
                             |      |          |         | sl         |
                             |      |          |         |            |
                             |      |<-------- | -----ls-|<-θSsp      |
                             |                 |<-----[K]-------------|<-wI
                             |<-----------------------[K]-------------|<-θI


        value='φI' (KwI = 0)

        <=4================================4============================
                ||                         m                          ||
                4 (m-mo) =======0=======                              ||
                ||    M  ||  (1-β)m   ||    C            S         I  ||
        θo,φo=>[MX1]==0==||          [MX2]==2==[HC]==F===3==>[TZ]==4==||
         mo              ||         s ||        /   /         //      |
                         ===0=[CC]==1===       s   m         sl       |
                              /\\   βm         |             ||       |
                             t  sl  |          |            [BL]<-mi  |
                             |      |          |            //        |
                             |      |          |           sl         |
                             |      |          |                      |
                             |      |<-------- | ------ls-------------|<-φI
                             |                 |<-----[K]-------------|<-wI
                             |<-----------------------[K]-------------|<-θI
        """
        # Design values
        self.actual[5:] = θo, φo, θIsp, φIsp, mi, UA, Qsa, Qla

        x = self.β_ls(value, sp)
        print('m = {m: .3f} kg/s, mo = {mo: .3f} kg/s, β = {β: .3f}'.format(
            m=self.actual[0], mo=self.actual[1], β=self.actual[2]))
        self.psy_chart(x, θo, φo)

        return x


# TESTS: uncomment
# Kθ, Kw = 1e10, 0     # Kw can be 0
# β = 0

# m, mo = 3.093, 0.94179
# θo, φo = 32, 0.5
# θIsp, φIsp = 26, 0.5

# mi = 15e3 / (l * (psy.w(θo, φo) - psy.w(θIsp, φIsp)))   # kg/s
# UA = 45e3 / (θo - θIsp) - mi * c                        # W
# Qsa, Qla = 34_000, 4_000     # W

# print(f'QsTZ = {(UA + mi * c) * (θo - θIsp): ,.1f} W')
# print(f'QlTZ = {mi * l * (psy.w(θo, φo) - psy.w(θIsp, φIsp)): ,.1f} W')


# TEST cool.ipynb
# =========================================
# 1.2 Create air handing unit (AHU) object
# Kθ, Kw = 1e10, 0     # Kw can be 0
# β = 0.2              # by-pass factor

# m, mo = 3.1, 1.      # kg/s, mass flow rate: supply & outdoor (fresh) air
# θo, φo = 32., 0.8    # outdoor conditions
# θIsp, φIsp = 26., 0.5        # set point for indoor condition

# mi = 1.35            # kg/s, mass flow rate of infiltration air
# UA = 675.            # W/K, overall heat coefficient of the building
# QsBL, QlBL = 34000., 4000.     # W, auxiliary loads: sensible & latent

# parameters = m, mo, β, Kθ, Kw
# inputs = θo, φo, θIsp, φIsp, mi, UA, QsBL, QlBL
# cool = MxCcRhTzBl(parameters, inputs)

# 2. system wthout reheating
# 2.1. CAV constant air volume
# print('\nCAV: m given')
# cool.CAV_wd()
# cool.CAV_wd(θo, φo, θIsp, φIsp, mi, UA, QsBL, QlBL)

# # 2.
# print('\nVAV: control θS')
# θSsp = 11.45
# # cool.VAV_wd('θS', θSsp)
# cool.VAV_wd('θS', θSsp, θo, φo, θIsp, φIsp, mi, UA, Qsa, Qla)
# # 3.
# print('\nVAV: control φI')
# cool.VAV_wd('φI', 0.5, θo, φo, θIsp, φIsp, mi, UA, Qsa, Qla)
# # 4.
# print('\nVBP: control φI')
# wIsp = psy.w(26, 0.5)
# m = 3.5
# cool.actual[0] = m
# # cool.VBP_wd('wI', wIsp)
# φI = 0.4
# x = cool.VBP_wd('φI', φI, θo, φo, θIsp, φIsp, mi, UA, Qsa, Qla)
# # cool.VBP_wd('wI', wIsp, θo, φo, θIsp, φIsp, mi, UA, Qsa, Qla)

# """
# Two solutions as a function of β0:
#     β0 = 0.1 ... 0.4 -> β = 0.101
#     β0 = 0.5 ... 0.9 -> β = 0.743
# Explanation
# The thermal zone characteristics cuts the saturation curve in two points.
# """

# """
# Controlling θS with β: NO SOLUTION
# Explanation
# Cannot have imposed m, θS, θI and QsTZ.
# Given θI and QsTZ:
#     either m --> θS
#     or θS --> m
# """
# # print('\nVBP: control θS')
# # θSsp = 11.77
# # m = 3.162
# # cool.actual[0] = m
# # cool.VBP_wd('θS', θSsp, θo, φo, θIsp, φIsp, mi, UA, Qsa, Qla)


# Kw = 1e10
# cool.actual[4] = Kw
# cool.CAV_wd(θo, φo, θIsp, φIsp, mi, UA, Qsa, Qla)
