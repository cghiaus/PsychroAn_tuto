"""
Created on Sun Apr  5 08:14:37 2020

@author: cghiaus
https://problemsolvingwithpython.com
Psycrometry
pvs(t)      pressure of saturated vapor
v(t, r)     specific volume

"""
import numpy as np


def pvs(t):
    """
    Saturation vapor pressure as a function of tempetature
    t [°C]
    """
    import numpy as np
    T = t + 273.15      # [K] Temperature
    # pws(T) saturation pressure over liquid water
    # for temp range [0 200] °C eq. (6)
    C8 = -5.8002206e3
    C9 = 1.3914993e0
    C10 = -4.8640239e-2
    C11 = 4.1764768e-5
    C12 = -1.4452093e-8
    C13 = 6.5459673e0
    y = np.exp(C8/T + C9 + C10*T + C11*T**2 + C12*T**3 + C13*np.log(T))  # Pa
    return y


def v(t, w, Z=0):
    """
    Specific volum as a function of température and humidity ratio
    for a given altitude (default 0 m)
    t : temperature [°C]
    w : humidity ratio [kg/kg_da]
    Z : altitude [m]; default value = 0
    """
    Mv = 18.01528       # [kg/kmol] vapor molaire mass
    Mda = 28.9645       # [kg/kmol] air molaire mass
    R = 8320            # [J/(kmol*K)] ideal gaz constant

    # Static pressure function of altitude;
    p = 101325*(1 - 2.25577e-5 * Z)**5.2559     # [Pa]
    v = R/Mv*(Mv/Mda + w)*(t + 273.15)/p
    return v


def w(t, phi, Z=0):
    """
    Humidity ratio as function of temperature and relative humidity
    t : temperature [°C]
    phi : relative humidity [-]
    Z : altitude [m]; default value = 0
    """
    phi = phi       # psi [-]
    Mv = 18.01528       # [kg/kmol] vapor molaire mass
    Mda = 28.9645       # [kg/kmol] air molaire mass

    # Static pressure function of altitude;
    p = 101325*(1 - 2.25577e-5 * Z)**5.2559     # [Pa]
    w = Mv/Mda*phi*pvs(t)/(p - phi*pvs(t))
    return w


def wsp(ts, p=101325):
    """
    Derivative of the saturation curve for temperature ts
    Parameters
    ----------
    ts : temperature on saturation curve [°C]
    p  : pressure [Pa]

    Returns
    -------
    wsp : value of the derivative of the function w(ts) Tetens eq.
    https://en.wikipedia.org/wiki/Tetens_equation
    """
    Mv = 18.01528       # [kg/kmol] vapor molaire mass
    Mda = 28.9645       # [kg/kmol] air molaire mass
    exp_t = np.exp(17.2694*ts/(ts + 238.3))
    # ws = Mv/Mda*610.78*exp_t/(p - 610.78*exp_t)
    wp = Mv/Mda*p*2.51354e6*exp_t/((ts + 238.3)**2*(p - exp_t)**2)
    return wp


def chart(t, w,
          t_range=np.arange(-10, 50, 0.1),
          w_range=np.arange(0, 0.030, 0.0001)):
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
    plt.xlabel(r'Temperature $\theta$ [°C]')
    ax.yaxis.set_label_position("right")
    plt.ylabel(r'Humidity ratio w [kg/kg]')
    plt.grid(True)
    plt.plot(t_range, psy.w(t_range, 100), linewidth=2)    # saturation curve

    # Plot relative humidity curves
    for phi in np.arange(0, 100, 20):
        w4t = psy.w(t_range, phi)
        plt.plot(t_range, w4t, linewidth=0.5)
        s_phi = "%3.0f" % phi
        ax.annotate(s_phi+' %', xy=(t_range[-1]-3, w4t[-1]))

    plt.plot(t, w, linewidth=3)    # processes
    return None


def chartA(t, wv, A,
           t_range=np.arange(-10, 50, 5),
           w_range=np.arange(0, 0.030, 0.01)):
    """
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

    """
    import matplotlib.pyplot as plt
    import psychro as psy

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.yaxis.tick_right()
    plt.xlabel(r'Temperature $\theta$ [°C]')
    ax.yaxis.set_label_position("right")
    plt.ylabel(r'Humidity ratio w [kg/kg]')
    plt.grid(True)
    plt.plot(t_range, psy.w(t_range, 1), linewidth=2)    # saturation curve

    # Plot relative humidity curves
    for phi in np.arange(0, 1, 0.2):
        w4t = psy.w(t_range, phi)
        plt.plot(t_range, w4t, linewidth=0.5)
        phi100 = phi*100
        s_phi = "%3.0f" % phi100
        ax.annotate(s_phi+' %', xy=(t_range[-1]-3, w4t[-1]))

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
