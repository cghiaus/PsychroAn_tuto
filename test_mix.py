#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Apr 30 07:55:42 2020

@author: cghiaus
test mix07.py
"""
import mix07 as mx
import numpy as np


def test_mean_temperature():
    """
    The name of the function needs to start with test_
    Returns
    -------
    None.

    """
    # Some inputs
    θ0, θ1 = 0, 32
    # Expected output
    ye = (θ0 + θ1) / 2

    # Obtained output: the temperature after MX
    yo = mx.mixing(1, θ0, 0.8, θ1, φ1=0.95, α=0.5)[0]

    np.testing.assert_almost_equal(yo, ye, 1)


def test_mix():
    # Some inputs
    t0, t1 = 0, 32
    # Expected output
    ye = (t0 + t1) / 2
    # Obtained output
    yo = mx.mixing(1, t0, 0.8, t1, φ1=0.95, α=0.5)[0]

    np.testing.assert_almost_equal(yo, ye, 1)


# test_mean_temperature()
