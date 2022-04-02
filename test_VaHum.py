"""
Created on Tue Apr 21 10:54:52 2020

@author: cghiaus
test winter_VaHum.py
https://www.spyder-ide.org/blog/introducing-unittest-plugin/
"""
# from winter_VaHum import ModelAllOutAir, AllOutAirCAV, ModelRecAir
import winter_VaHum as vh

import numpy as np


def test_vhModelAllOutAir():
    """
    --o->HC--0->VH--1->TZ--2-->
         /       /     ||  |
         |       |     BL  |
         |       |         |
         |       |<----Kw--|-w2
         |<------------Kt--|-t2
    """
    # Input data
    m, tS = 4.84, 30            # supply air mass flow kg/s; temp. °C
    tIsp, phiIsp = 18, 0.50     # indoor set point temp. °C; RH -
    tO, phiO = -1, 1            # outdoor temp. °C; RH -
    Qsa, Qla = 0, 0             # aux. sens. ; latent loads W
    mi, UA = 2.12, 935.83       # bldg. inflitr. k/s ; conductivity W/K

    # Expected results
    y = np.array([30, 3.5076e-3,            # point 0 (t, w)
                  30, 7.670e-3,             # point 1 (t, w)
                  18, 6.4025e-3,            # point 2 (t, w)
                  150020.8, 50290.27,       # QsHC, QlVH
                  -58060.78, -15318.30])    # QsTZ, QlTZ
    np.testing.assert_almost_equal(
        vh.ModelAllOutAir(m, tS, tIsp, phiIsp, tO, phiO, Qsa, Qla, mi, UA),
        y, 1)


def test_vhAllOutAirCAV():
    """
    --o->HC--0->VH--1->TZ--2-->
         /       /     ||  |
         |       |     BL  |
         |       |         |
         |       |<----Kw--|-w2
         |<------------Kt--|-t2
    """
    # Expected results
    y = np.array([29.77, 3.5076e-3,         # point 0 (t0, w0)
                  29.77, 7.6465e-3,         # point 1 (t1, w1)
                  18, 6.4025e-3,            # point 2 (t2, w2)
                  151795.37, 50965.12,      # QsHC, QlVH
                  -58060.78, -15318.30])    # QsTZ, QlTZ
    np.testing.assert_almost_equal(
        vh.AllOutAirCAV(tS=30, tIsp=18, phiIsp=0.5, tO=-1, phiO=1,
                     Qsa=0, Qla=0, mi=2.12, UA=935.83),
        y, 2)


def test_vhModelRecAir():
    """
    <-3--|<-------------------------|
         |                          |
    -o->MX--0->HC--1->VH--2->TZ--3-->
               /       /     ||  |
               |       |     BL  |
               |       |         |
               |       |<----Kw--|-w3
               |<------------Kt--|-t3
    """
    # Input data
    m = 4.84                    # supply air mass flow kg/s
    alpha = 1                   # 100 % outdoor air
    tS = 30                     # supply air °C
    tIsp, phiIsp = 18, 0.50     # indoor set point temp. °C; RH -
    tO, phiO = -1, 1            # outdoor temp. °C; RH -
    Qsa, Qla = 0, 0             # aux. sens. ; latent loads W
    mi, UA = 2.12, 935.83       # bldg. inflitr. k/s ; conductivity W/K

    # Expected results
    y = np.array([-1, 3.5076e-3,            # point 0 (t0, w0)
                  29.996, 3.5076e-3,        # point 1 (t1, w1)
                  29.996, 7.6611e-3,        # point 2 (t2, w2)
                  18, 6.3960e-3,            # point 3 (t3, w3)
                  150020.65, 50176.49,            # QsHC, QlVH
                  -58060.72, -15283.64])          # QsTZ, QlTZ
    np.testing.assert_almost_equal(
        vh.ModelRecAir(m, alpha, tS, tIsp, phiIsp, tO, phiO, Qsa, Qla, mi, UA),
        y, 2)
