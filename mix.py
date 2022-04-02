#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 09:20:28 2020
Updated on Sat Apr  2 15:37:19 2022

@author: cghiaus
To be used with Jupyter Notebook: T03_Mix.ipynb

Plots:
- point 3 if mixing with condensation;
- point 2 if mixing w/o condensation
by selecting from two diffrent models.

Two models are calculated

    |                   |
    1                   1
    |                   |
-0->MX-2->AD-3->    -0->MX-2->

A function to
solve the linear system of equations A * x = b:

m * c * θ2 = α * c * θ0 + (1 - α) * m * c * θ1  # [MX] sensible
m * l * w2 = α * l * w0 + (1 - α) * m * l * w1  # [MX] latent
c * θ2 + l * w2  - c * θ3 - l * w3 = 0          # [AH] h const.
ws'(θ3) * θ3 - w3 = ws'(θ30) * θ30 - ws(θ30)    # [AH] saturation curve

for: x = [θ2, w2, θ3, w3]. (Note: point 3 is on saturation curve)

The systems is solved iterativelly till abs(psy.w(θ30, 1) - x[3]) is small.

Plot results on psychrometric chart
"""
import numpy as np
import psychro as psy


# A * x = b ==> x = inv(A) * b
# constants
c = 1e3                     # air specific heat J/kg K
l = 2496e3                  # latent heat J/kg


def mixing(m=1, θ0=3, φ0=0.8, θ1=32, φ1=0.95, α=0.5):
    """
    Adiabatic mixing.
    If the point is in oversaturation, then adiabatic condensation.
    """
    w0 = psy.w(θ0, φ0)
    w1 = psy.w(θ1, φ1)

    def MX():
        """
        Mixing with given ratio
        """
        A = np.zeros((2, 2))            # coefficents of unknowns
        b = np.zeros(2)                 # vector of inputs
        A[0, 0] = m * c
        b[0] = α * m * c * θ0 + (1 - α) * m * c * θ1

        A[1, 1] = m * l
        b[1] = α * m * l * w0 + (1 - α) * m * l * w1

        x = np.linalg.solve(A, b)
        return x

    def MX_AD():
        """
        Mixing with given ration
        Adiabatic humidification / condensation
        """
        A = np.zeros((4, 4))            # coefficents of unknowns
        b = np.zeros(4)                 # vector of inputs
        A[0, 0] = m * c
        b[0] = α * m * c * θ0 + (1 - α) * m * c * θ1

        A[1, 1] = m * l
        b[1] = α * m * l * w0 + (1 - α) * m * l * w1

        A[2, 0], A[2, 1], A[2, 2], A[2, 3] = c, l, -c, -l
        b[2] = 0

        A[3, 2], A[3, 3] = psy.wsp(θs0), -1
        b[3] = psy.wsp(θs0) * θs0 - psy.w(θs0, 1)

        x = np.linalg.solve(A, b)
        return x

    θs0, Δ_θs = (θ0 + θ1) / 2, 2             # initial guess saturation temp.

    x = MX_AD()
    if x[1] > psy.w(x[0], 1):
        # Model MX & AD
        while Δ_θs > 0.01:
            x = MX_AD()
            Δ_θs = abs(θs0 - x[2])
            θs0 = x[2]

        print(f'θ2 = {x[0]:5.2f} °C, w2 = {1000*x[1]:5.2f} g/kg')
        print(f'θ3 = {x[2]:5.2f} °C, w3 = {1000*x[3]:5.2f} g/kg')

        # Processes on psychrometric chart
        # Points        o   i  0  1     Elements
        #               0   1  1  3
        A = np.array([[-1, -1, 1, 0],   # MX
                      [0, 0, -1, 1]])   # AD
        θ = np.append([θ0, θ1], x[0:4:2])
        w = np.append([w0, w1], x[1:4:2])
        psy.chartA(θ, w, A)
    else:
        x = MX()
        print(f'θ2 = {x[0]:5.2f} °C, w2 = {1000*x[1]:5.2f} g/kg')
        print('---')
        # Processes on psychrometric chart
        # Points        o   i  0        Elements
        #               0   1  2
        A = np.array([[-1, -1, 1]])     # MX
        θ = np.array([θ0, θ1, x[0]])
        w = np.array([w0, w1, x[1]])
        psy.chartA(θ, w, A)
    return


# Uncommet next line in order to test
# x = mixing(m=1, θ0=0, φ0=0.8, θ1=32, φ1=0.95, α=0.5)
# >>> θ2 = 16.00 °C, w2 = 16.03 g/kg
# >>> θ3 = 19.80 °C, w3 = 14.51 g/kg
# These are the outputs given by function
# test_mean_temperature() du test_mix07.py

# Uncomment next line in order to test
# x = mixing(m=1, θ0=10, φ0=0.8, θ1=32, φ1=0.95, α=0.5)
