"""
Created on Wed May 25 21:01:26 2022
@author: cghiaus

Psychrometry
pvs(t)      pressure of saturated vapor
v(t, r)     specific volume
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Constants
Mv = 18.015_286             # [kg/kmol] vapor molaire mass
Mda = 28.966                # [kg/kmol] air molaire mass
R = 8_314.462_618_153_24    # [J/(kmol*K)] ideal gaz constant


def p(Z):
    # Atmospheric pressure function of altitude (-500 .. 10000 m):
    return 101325 * (1 - 2.25577e-5 * Z)**5.2559     # [Pa]


def pvs(t):
    """
    Saturation vapor pressure as a function of tempetature
    t [°C]
    """
    import numpy as np
    T = t + 273.15      # [K] Temperature
    # pws(T) [Pa] saturation pressure over liquid water
    # for temp range [0 200] °C eq. (6)
    C8 = -5.800_220_6e3
    C9 = 1.391_499_3e0
    C10 = -4.864_023_9e-2
    C11 = 4.176_476_8e-5
    C12 = -1.445_209_3e-8
    C13 = 6.545_967_3e0
    pws = np.exp(
        C8 / T + C9 + C10 * T + C11 * T**2 + C12 * T**3 + C13 * np.log(T))
    return pws


def w(t, phi, Z=0):
    """
    Humidity ratio as a function of temperature and relative humidity
    t : temperature [°C]
    phi : relative humidity [-]
    Z : altitude [m]; default value = 0
    """
    w = Mv / Mda * phi * pvs(t) / (p(Z) - phi * pvs(t))
    return w


def wsp(ts, Z=0):
    """
    Derivative of the saturation curve for temperature ts
    wsp = Humidity ratio (w) at saturation (s) - derivative (prime)
    Parameters
    ----------
    ts : temperature on saturation curve [°C]
    p  : pressure [Pa]
    Returns
    -------
    wsp : value of the derivative of the function w(ts) Tetens eq.
    Murray F.W. (1967) On the computation of saturation vapour pressure.
    J. Applied Meteorology 6: 203-204
    """
    # equilibrum pressure at saturation eq.(6)
    a = 17.2693882
    b = 273.16 - 35.86
    C = 610.78
    es = C * np.exp(a * ts / (ts + b))
    # ws = Mv/Mda*610.78*exp_t/(p - 610.78*exp_t)
    wp = Mv / Mda * a * b * p(Z) * es / (
        (ts + b)**2 * (p(Z) - es)**2)
    return wp


def v(t, w, Z=0):
    """
    Specific volume as a function of température and humidity ratio
    for a given altitude (default 0 m)
    t : temperature [°C]
    w : humidity ratio [kg/kg_da]
    Z : altitude [m]; default value = 0
    """
    v = R / Mv * (Mv / Mda + w) * (t + 273.15) / p(Z)
    return v


def phi(t, w, Z=0):
    """
    Relative humidity as a function of temperature and humidity ratio
    for a given altitude (default 0 m)
    t : temperature [°C]
    w : humidity ratio [kg/kg_da]
    Z : altitude [m]; default value = 0
    """
    phi = p(Z) / pvs(t) * w / (Mv / Mda + w)
    return phi


def t(w, phi, Z=0, tol=1e-6, max_iter=100):
    """
    Temperature as a function of humidity ratio and relative humidity
    (Inverse function of w(t, phi, Z))

    w : humidity ratio [kg/kg]
    phi : relative humidity [-]
    Z : altitude [m]; default value = 0
    tol : tolerance for convergence
    max_iter : maximum number of iterations

    Returns temperature t in °C
    """

    def equation(t_guess):
        p_ws = pvs(t_guess)
        return (Mv / Mda) * phi * p_ws / (p(Z) - phi * p_ws) - w

    # Initial guess (20°C in Kelvin)
    t_initial_guess = 20.0

    # Solve numerically
    t_solution = fsolve(equation, t_initial_guess,
                        xtol=tol, maxfev=max_iter)[0]

    return t_solution


def chart(t, w,
          t_range=np.arange(-10, 50, 0.01),
          w_range=np.arange(0, 0.030, 0.00001)):
    """
    Parameters
    ----------
    t_range : temperature vector t = np.arange(-10, 50, 0.1)
    w_range : humidity ration vector w = np.arange(0, 0.030, 0.0001)
    Returns
    -------
    None. Psycrometric chart
    """

    import matplotlib.pyplot as plt
    import psychro as psy

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.yaxis.tick_right()
    plt.xlabel(r"Temperature $\theta$ / °C")
    ax.yaxis.set_label_position("right")
    plt.ylabel(r"Humidity ratio w / kg·kg⁻¹")
    plt.grid(True)
    plt.plot(t_range, psy.w(t_range, 100), linewidth=2)    # saturation curve

    # Plot relative humidity curves
    for phi in np.arange(0, 100, 20):
        w4t = psy.w(t_range, phi)
        plt.plot(t_range, w4t, linewidth=0.5)
        s_phi = "%3.0f" % phi
        ax.annotate(s_phi + ' %', xy=(t_range[-1] - 3, w4t[-1]))

    plt.plot(t, w, linewidth=3)    # processes
    plt.draw()
    plt.show()
    return None


def chartA(t, wv, A,
           t_range=np.arange(-10, 50, 5),
           w_range=np.arange(0, 0.030, 0.01)):
    """
    Plots processes on psychrometric chart.

    Parameters
    ----------
    t : np.array, no. equal to no. points in the psy-chart
        temperatures, °C
    wv: np.array, wv.shape = t.shape
        weight vapor, kg/kg_da
    A : np.array [no. processes, no. points = no. temperatures]
        adjancy matrix: -1 flow our of node, 1 flow in node, 0 no connection
    t_range : np.arange
        range of temperature
        the default is np.arange(-10, 50, 0.1).
        temperature vector t = np.arange(-10, 50, 0.1)
    w_range : np.arange
        humidity ration vector
        The default is np.arange(0, 0.030, 0.01).
    Returns
    -------
    None.

    Example
    -------
    # Points        0   1  2  3         Process
    A = np.array([[-1, -1, 1, 0],       # MX mixing
                  [0, 0, -1, 1]])       # AD adiabatic
    θ = np.array([θ0, θ1, θ2, θ3])
    w = np.array([w0, w1, w2, w3])
    psy.chartA(θ, w, A)
    """
    # import matplotlib.pyplot as plt
    # import psychro as psy

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.yaxis.tick_right()
    plt.xlabel(r'Temperature $\theta$ / °C')
    ax.yaxis.set_label_position("right")
    plt.ylabel(r'Humidity ratio w / kg·kg⁻¹')
    plt.grid(True)
    plt.plot(t_range, w(t_range, 1), linewidth=2)    # saturation curve

    # Plot relative humidity curves
    for phi in np.arange(0, 1, 0.2):
        w4t = w(t_range, phi)
        plt.plot(t_range, w4t, linewidth=0.5)
        phi100 = phi * 100
        s_phi = "%3.0f" % phi100
        ax.annotate(s_phi + ' %', xy=(t_range[-1] - 3, w4t[-1]))

    for k in range(0, A.shape[0]):
        tk = np.nonzero(A[k, :])
        wk = np.nonzero(A[k, :])
        plt.plot(t[tk], wv[wk], linewidth=3)    # processes
        # plot no. point
        for j in range(0, np.shape(tk)[1]):
            plt.text(t[tk][j], wv[tk][j], str(tk[0][j]))
    plt.draw()
    plt.show()
    return None
