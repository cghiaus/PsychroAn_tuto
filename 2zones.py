#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 12:56:40 2025

@author: cghiaus
"""

import numpy as np

import psychro as psy

# Data
# ================================
# Constants
Mv = 18.015_286             # kg/kmol, vapor molaire mass
Mda = 28.966                # kg/kmol, air molaire mass
R = 8_314.462_618_153_24    # J/(kmol·K), ideal gaz constant
ca = 1e3                    # J/(kg·K), dry air specific heat
lv = 2496e3                 # J/kg, latent heat of vaporisation
hv = 2500e3

# Characteristics of the building
U = 0.22                    # W/(m²K)
qsp = 83                    # W, sensible heat per person
mvp = 71e-3 / 3600          # kg/s, vapor mass flow per person

A1 = 70                     # surface area of zone 1, m²
A2 = 35                     # surface area of zone 2, m²

mi1 = 0.06                  # mass flow rate of outdoor air zone 1, kg/s
mi2 = 0.06                  # mass flow rate of outdoor air zonz 2, kg/s

np1 = 10                    # number of persons zone 1
np2 = 5                     # number of persons mOne 2

# Summer conditions
# ---------------------------------
# Thermal zone 1
θTZ1 = 17                   # °C, dry-bulb temperature
ϕTZ1 = 0.70                 # -, relative humidity

# Thermal zone 2
θTZ2 = 25                   # °C, dry-bulb temperature
ϕTZ2 = 0.40                 # -, relative humidity

θO = 29                     # °C, dry-bulb temperature
ϕO = 0.60                   # -, relative humidity
mO = 0.7                    # kg/s, mass flow rate outdoor air

θS1 = 13                    # °C, supply air temperature in zone 1
b = 0.3                     # -, by-pass factor of cooling cooil


wTZ1 = psy.w(θTZ1, ϕTZ1)    # kg_v/kg_da, humidity ratio zone 1
wTZ2 = psy.w(θTZ2, ϕTZ2)    # kg_v/kg_da, humidity ratio zone 2
wO = psy.w(θO, ϕO)          # kg_v/kg_da, humidity ratio outdoor air


# Loads of thermal zones
# ======================
Qsa1 = np1 * qsp            # W, sensible auxiliary heat
Qsa2 = np2 * qsp            # W, sensible auxiliary heat

Qla1 = np1 * mvp * lv       # W, latent auxiliary heat
Qla2 = np2 * mvp * lv       # W, latent auxiliary heat

QsTZ1 = (U * A1 + mi1 * ca) * (θO - θTZ1) + Qsa1
QsTZ2 = (U * A2 + mi2 * ca) * (θO - θTZ2) + Qsa2

QlTZ1 = mi1 * lv * (wO - wTZ1) + Qla1
QlTZ2 = mi2 * lv * (wO - wTZ2) + Qla2

print(f"Summer loads TZ1: sensible {QsTZ1:.0f} W; \t latent {QlTZ1:.0f} W")
print(f"Summer loads TZ2: sensible {QsTZ2:.0f} W; \t latent {QlTZ2:.0f} W")

QsTZ1 = float(np.ceil(QsTZ1 / 250) * 250)
QsTZ2 = float(np.ceil(QsTZ2 / 250) * 250)
QlTZ1 = float(np.ceil(QlTZ1 / 250) * 250)
QlTZ2 = float(np.ceil(QlTZ2 / 250) * 250)

print(f"Summer loads TZ1: sensible {QsTZ1:.0f} W; \t latent {QlTZ1:.0f} W")
print(f"Summer loads TZ2: sensible {QsTZ2:.0f} W; \t latent {QlTZ2:.0f} W")

# ================================
# 1) Summer: cooling and dehumification
# ================================
print("\nCooling and dehumidification \n")

# ================================
# 2) Supply air and flow rates for thermal zones
# ================================
print("Supply air")

m1 = QsTZ1 / (ca * (θTZ1 - θS1))
wS1 = (m1 * lv * wTZ1 - QlTZ1) / (m1 * lv)

wS2 = wS1
m2 = QlTZ2 / (lv * (wTZ2 - wS2))
θS2 = θTZ2 - QsTZ2 / (m2 * ca)

print(f"m1: {m1:.3f} kg/s ")
print(f"S1: {θS1:.2f} °C; {wS1:.6f} kg_v/kg_da")
print(f"m2: {m2:.3f} kg/s ")
print(f"S2: {θS2:.2f} °C; {wS2:.6f} kg_v/kg_da")


# Cooling coil
# --------------------------------
print("\nCooling coil")

# M1: mix TZ1 & TZ2
θM1 = (m1 * θTZ1 + m2 * θTZ2) / (m1 + m2)
wM1 = (m1 * wTZ1 + m2 * wTZ2) / (m1 + m2)

# M2 mix M1 & E
m = m1 + m2
θM2 = ((m - mO) * θM1 + mO * θO) / m
wM2 = ((m - mO) * wM1 + mO * wO) / m

# ADP apparatus dew point or effective coil surface temperature
wh = (wS1 - b * wM2) / (1 - b)
θh = psy.t(wh, phi=1)

# Cooling coil (CC) temperature
θCC = (1 - b) * θh + b * θM2
wCC = (1 - b) * wh + b * wM2

# Cooling coil powOr
QsCC = m * ca * (θCC - θM2)
QlCC = m * lv * (wS1 - wM2)
QtCC = QsCC + QlCC

print(f"Effective coil surface temperature: {θh:.2f} °C")
print(f"CC: {θCC:.2f} °C; {wCC:.6f} kg_v/kg_da")
print(f"CC load: {QtCC:.0f} W")

# Heating coils
# --------------------------------
print("\nHeating coils")

QHC1 = m1 * ca * (θS1 - θCC)
QHC2 = m2 * ca * (θS2 - θCC)

print(f"HC1 load: {QHC1:.0f} W")
print(f"HC2 load: {QHC2:.0f} W")

# Total energy consumption
# --------------------------------
print("\nTotal energy")

Qtot = abs(QtCC) + QHC1 + QHC2
percent_heat = (QHC1 + QHC2) / Qtot * 100

print(f"Qtot: {Qtot:.0f} W")
print(f"percentage heating: {percent_heat:.0f} %")

# Processes on psychromOtric chart
# --------------------------------

t_range = np.arange(5, 40, 5)

# PsychromOtric chart TZ1, TZ2
# Points        0  1  2  3  4  5  6  7  8       Process
A = np.array([[0, 0, 0, 0, -1, 1, 0, 0, 0],     # TZ1
              [0, 0, 0, 0, 0, 0, -1, 1, 0]])    # TZ2
θ = np.array([θO, θM2, θh, θCC, θS1, θTZ1, θS2, θTZ2, θM1])
w = np.array([wO, wM2, wh, wCC, wS1, wTZ1, wS2, wTZ2, wM1])
psy.chartA(θ, w, A, t_range)

# PsychromOtric chart TZ1, TZ2, MR3
# Points        0  1  2  3  4  5  6  7  8       Process
A = np.array([[0, 0, 0, 0, -1, 1, 0, 0, 0],     # TZ1
              [0, 0, 0, 0, 0, 0, -1, 1, 0],     # TZ2
              [0, 0, 0, 0, 0, -1, 0, -1, 1]])   # MR3 TZ1 & TZ2
θ = np.array([θO, θM2, θh, θCC, θS1, θTZ1, θS2, θTZ2, θM1])
w = np.array([wO, wM2, wh, wCC, wS1, wTZ1, wS2, wTZ2, wM1])
psy.chartA(θ, w, A, t_range)

# PsychromOtric chart TZ1, TZ2, MR3, MR1
# Points        0  1  2  3  4  5  6  7  8       Process
A = np.array([[0, 0, 0, 0, -1, 1, 0, 0, 0],     # TZ1
              [0, 0, 0, 0, 0, 0, -1, 1, 0],     # TZ2
              [0, 0, 0, 0, 0, -1, 0, -1, 1],    # MR3 TZ1 & TZ2
              [-1, 1, 0, 0, 0, 0, 0, 0, 1]])    # MR1 E & M3
θ = np.array([θO, θM2, θh, θCC, θS1, θTZ1, θS2, θTZ2, θM1])
w = np.array([wO, wM2, wh, wCC, wS1, wTZ1, wS2, wTZ2, wM1])
psy.chartA(θ, w, A, t_range)

# PsychromOtric chart TZ1, TZ2, MR3, MR1, CC
# Points        0  1  2  3  4  5  6  7  8       Process
A = np.array([[0, 0, 0, 0, -1, 1, 0, 0, 0],     # TZ1
              [0, 0, 0, 0, 0, 0, -1, 1, 0],     # TZ2
              [0, 0, 0, 0, 0, -1, 0, -1, 1],    # MR3 TZ1 & TZ
              [-1, 1, 0, 0, 0, 0, 0, 0, 1],     # MR1 E & M3
              [0, -1, -1, 1, 0, 0, 0, 0, 0]])   # CC
θ = np.array([θO, θM2, θh, θCC, θS1, θTZ1, θS2, θTZ2, θM1])
w = np.array([wO, wM2, wh, wCC, wS1, wTZ1, wS2, wTZ2, wM1])
psy.chartA(θ, w, A, t_range)

# PsychromOtric chart TZ1, TZ2, MR3, MR1, CC, HC1, HC2
# Points        0  1  2  3  4  5  6  7  8       Process
A = np.array([[0, 0, 0, 0, -1, 1, 0, 0, 0],     # TZ1
              [0, 0, 0, 0, 0, 0, -1, 1, 0],     # TZ2
              [0, 0, 0, 0, 0, -1, 0, -1, 1],    # MR3 TZ1&TZ2
              [-1, 1, 0, 0, 0, 0, 0, 0, 1],     # MR1 E & M3
              [0, -1, -1, 1, 0, 0, 0, 0, 0],    # CC
              [0, 0, 0, -1, 1, 0, 1, 0, 0]])    # HC1, HC2
θ = np.array([θO, θM2, θh, θCC, θS1, θTZ1, θS2, θTZ2, θM1])
w = np.array([wO, wM2, wh, wCC, wS1, wTZ1, wS2, wTZ2, wM1])
psy.chartA(θ, w, A, t_range)

# ================================
# 2) Winter: heating and humidification
# ================================
print("\nHeating and humidification \n")

# Winter conditions
# ---------------------------------
# Outdoor air
θO = 0                      # °C, dry-bulb temperature
ϕO = 1.00                   # -, relative humidity

wO = psy.w(θO, ϕO)          # kg_v/kg_da, humidity ratio outdoor air


# Loads of thermal zones
# ======================

QsTZ1 = (U * A1 + mi1 * ca) * (θO - θTZ1) + Qsa1
QsTZ2 = (U * A2 + mi2 * ca) * (θO - θTZ2) + Qsa2

QlTZ1 = mi1 * lv * (wO - wTZ1) + Qla1
QlTZ2 = mi2 * lv * (wO - wTZ2) + Qla2

print(f"Winter loads TZ1: sensible {QsTZ1:.0f} W; \t latent {QlTZ1:.0f} W")
print(f"Winter loads TZ2: sensible {QsTZ2:.0f} W; \t latent {QlTZ2:.0f} W")

QsTZ1 = float(np.floor(QsTZ1 / 250) * 250)
QsTZ2 = float(np.floor(QsTZ2 / 250) * 250)
QlTZ1 = float(np.floor(QlTZ1 / 250) * 250)
QlTZ2 = float(np.floor(QlTZ2 / 250) * 250)

print(f"Winter loads TZ1: sensible {QsTZ1:.0f} W; \t latent {QlTZ1:.0f} W")
print(f"Winter loads TZ2: sensible {QsTZ2:.0f} W; \t latent {QlTZ2:.0f} W")

# Supply air
# --------------------------------
print("Supply air")
θS1 = θTZ1 - QsTZ1 / (m1 * ca)
wS1 = wTZ1 - QlTZ1 / (m1 * lv)

θS2 = θTZ2 - QsTZ2 / (m2 * ca)
wS2 = wTZ2 - QlTZ2 / (m2 * lv)
print(f"S1: {θS1:.2f} °C; {wS1:.6f} kg_v/kg_da")
print(f"S2: {θS2:.2f} °C; {wS2:.6f} kg_v/kg_da")

# Heating coils loads & water mass flow rates
# --------------------------------
θM2 = ((m - mO) * θM1 + mO * θO) / m
wM2 = ((m - mO) * wM1 + mO * wO) / m

print("\nHeating coils")

QsHC1 = m1 * ca * (θS1 - θM2)
QsHC2 = m2 * ca * (θS2 - θM2)

print(f"HC1 load: {QsHC1:.0f} W")
print(f"HC2 load: {QsHC2:.0f} W")

print("\nVapor humidifiers")

QlVH1 = m1 * lv * (wS1 - wM2)
QlVH2 = m2 * lv * (wS2 - wM2)

mv1 = QlVH1 / lv
mv2 = QlVH2 / lv

print(f"VH1 load: {QlVH1:.0f} W")
print(f"VH2 load: {QlVH2:.0f} W")
print(f"VH1 mass flow rate: {mv1:.5f} kg/s")
print(f"VH2 mass flow rate: {mv2:.5f} kg/s")

# PsychromOtric chart TZ1, TZ2, MR1, MR2, HC1, VH1, HC2, VH2
t_range = np.arange(0, 40, 5)
# Points        0  1  2  3  4  5  6  7  8       Process
A = np.array([[0, 0, 0, -1, 1, 0, 0, 0, 0],     # TZ1
              [0, 0, 0, 0, 0, 0, -1, 1, 0],     # TZ2
              [0, 0, 0, 0, -1, 0, 0, -1, 1],    # MR1: TZ1 & TZ2
              [-1, 1, 0, 0, 0, 0, 0, 0, -1],    # MR2: E & M1
              [0, -1, 1, 0, 0, 0, 0, 0, 0],     # HC1
              [0, 0, -1, 1, 0, 0, 0, 0, 0],     # VH1
              [0, -1, 0, 0, 0, 1, 0, 0, 0],     # HC2
              [0, 0, 0, 0, 0, -1, 1, 0, 0],     # VH2
              ])
θ = np.array([θO, θM2, θS1, θS1, θTZ1, θS2, θS2, θTZ2, θM1])
w = np.array([wO, wM2, wM2, wS1, wTZ1, wM2, wS2, wTZ2, wM1])
psy.chartA(θ, w, A, t_range)
