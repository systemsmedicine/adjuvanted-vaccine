# -*- coding: utf-8 -*-
"""
@author: cparrarojas

This module contains the PDE model of affinity maturation and RMSE cost
functions for parameter estimation.

The integration of the PDEs is adapted from scipy.integrate's reference guide:
    https://docs.scipy.org/doc/scipy/reference/tutorial/integrate.html

Implementation of the cost function with data for the parameter fitting
procedure is adapted from:
    https://stackoverflow.com/questions/32302654/scipy-differential-evolution
"""

# We import the necessary packages
import numpy as np
from scipy.integrate import odeint

# Global parameters of the model
dt = 0.05  # Time-step in days
timeStop = 133  # Time, in days, at which the simulation stops
t = np.arange(0, timeStop + dt, dt)  # Time domain
t_boost = 21.0  # Day of immunisation boost

Nx = 1500  # Number of points in x-axis discretisation
dx = 1.0/Nx  # Step-size in x
grid = np.arange(0.5*dx, 1.0 + 0.5*dx, dx)  # Numpy array of x values

sigma_N = 1.0  # Naive B cell stimulation in 1/days
sigma_B = 3.0  # GC B cell stimualtion in 1/days
r = 3.0  # Stimulated B cell proliferation rate in 1/days
rtilde = sigma_B*r/(sigma_B + r) # Effective B cell proliferation rate in 1/days
sigma_M = 1.0  # Memory B cell stimulation rate in 1/days
kappa = 5000.0  # GC carrying capacity
c = rtilde/kappa  # Competition factor
nu = 0.0001  # Variance of normal distribution representing mutation
D = 0.5*rtilde*nu  # Diffusion coefficient
g_B = 1.0/4.5  # B cell base decay rate in 1/days
g_B_1 = 161.222 # B cell decay attenuation
delta = 0.1  # B cell differentiation probability
k_Ab = 1.0  # Ab production rate in 1/days
ktilde = sigma_B*k_Ab/(sigma_B + k_Ab) # Effective Ab production rate in 1/days
g_Ab = 0.1  # Ab decay rate in 1/days
epsilon = 10**10 # B cell 'enhancement factor'
Htilde = 100 # Height of naive B cell pool

# Base binding affinity
def Q0(x, dmax):

    """
    Computes binding affinity between a given B cell and a virus protein at
    antigenic distance x.
    dmax -> affinity cutoff, x > dmax results in zero binding affinity.
    """

    if x <= dmax:
        return epsilon**(-x)
    else:
        return 0.0


# We vectorise the affinity function in order to apply it element-wise
vQ0 = np.vectorize(Q0)


# PDE system
def affinityMaturation(y, t, boost, H, baseQ, Q, ktilde, mu, dx):

    """
    PDEs for B cells and Ab. Takes a numpy array y = (B0, Ab0, B1, Ab1,...).
    Returns a numpy array dydt = (dBdt0, dAbdt0, dBdt1, dAbdt1,...).
    Q is the numpy array of reaction affinities Q_B(x)
    """

    # We obtain views of B and Ab by slicing y
    B = y[:: 2]
    Ab = y[1 :: 2]

    # Activity level G(t)
    if t<boost:
        g = np.exp(-mu*t)
    else:
        g = np.exp(-mu*(t - boost))

    # We create the return variable dydt
    dydt = np.empty_like(y)

    # We slice dydt to get dBdt and dAbdt
    dBdt = dydt[:: 2]
    dAbdt = dydt[1 :: 2]

    # We separate the extreme points of B, Q and H
    B1 = B[1 : -1]
    Q1 = Q[1 : -1]
    H1 = H[1 : -1]
    baseQ1 = baseQ[1 : -1]

    # Compute dB/dt and dAb/dt.  Non-flux boundary conditions are used for B.
    dBdt[0] = (g*(sigma_N*H[0] + (rtilde + sigma_M*Q[0])*B[0]
               + D*(2.0*B[1] - 2.0*B[0])/dx**2) - B[0]*(g_B/(1.0
               + g_B_1*baseQ[0]) + c*B[0]))

    dBdt[1:-1] = (g*(sigma_N*H1 + (rtilde + sigma_M*Q1)*B1
                  + D*np.diff(B,2)/dx**2) - B1*(g_B/(1.0 + g_B_1*baseQ1)
                   + c*B1))

    dBdt[-1] = (g*(sigma_N*H[-1] + (rtilde + sigma_M*Q[-1])*B[-1]
                + D*(2.0*B[-2] - 2.0*B[-1])/dx**2) - B[-1]*(g_B/(1.0
                + g_B_1*baseQ[-1]) + c*B[-1]))

    dAbdt[:] = Q*delta*ktilde*B - g_Ab*Ab

    return dydt


# Cost functions
def rmse(parameters, *data):

    """
    Integrates the PDEs and computes the root mean square of erros between the
    model and the data for a given a set of parameter values.
    Non-adjuvanted case.
    Parameters:
    gammaNA and gammaHA are the immunigenicities of NA and HA, respectively.
    mu is the activity level decay rate in 1/days.
    dmax is the affinity cutoff.
    """

    # We take the parameter values and data from the arguments.
    gammaNA, gammaHA, mu, dmax = parameters
    x, y = data

    # We create the naive B cell 'pool'.
    H = Htilde*0.5*(np.sign(grid - 0.99*dmax) + np.sign(1.0 - 0.99*dmax - grid))


    # We compute the numpy array of base and reaction affinities for the given
    # immunogenicities using the vectorised version of Q0.
    baseQ = vQ0(np.abs(grid), dmax) + vQ0(np.abs(1 - grid), dmax)
    Q = gammaNA*vQ0(np.abs(grid), dmax) + gammaHA*vQ0(np.abs(1 - grid), dmax)

    # Numpy array of values for B, Ab in the form (B0, Ab0, B1, Ab1...).
    # The initial condition corresponds to # B(x, 0) = 0.0, Ab(x, 0) = 1.0
    y0 = np.zeros(2*Nx)
    y0[1 :: 2] = np.ones(Nx)

    # We integrate the PDEs
    sol = odeint(affinityMaturation, y0, t, args=(t_boost, H, baseQ, Q, ktilde,
                 mu, dx), ml=2, mu=2)

    # Numpy array with the errors squared
    eSq = np.array([(np.sum(sol[np.argwhere(t == x[i])[0][0]][1 :: 2])*dx
                   - y[i])**2 for i in range(len(x))])

    return np.mean(eSq)**0.5


def rmse_adj(parameters, *data):

    """
    Integrates the PDEs and computes the root mean square of erros between the
    model and the data for a given a set of parameter values.
    Adjuvanted case: enhanced immunogenicities, cross-reactivity and antibody
    production.
    Parameters:
    betaNA and betaHA are the enhancement factors for the immunogenicities,
    with gamma -> gamma*beta for both NA and HA, and betaNA > betaHA
    reflecting enhanced cross-reactivity.
    betaAb is the enhancement factor for the production of antibodies, with
    k_Ab -> k_Ab*betaAb.
    """

    # We take the parameter values and data from the arguments.

    betaNA, betaHA, betaAb = parameters
    x, y, base_params = data
    gammaNA, gammaHA, mu, dmax, baseQ, H = base_params

    # We compute the numpy array of reaction affinities for the given
    # immunogenicities and adjuvanticity parameters using the vectorised version
    # of Q0.
    Q = (gammaNA*betaNA*vQ0(np.abs(grid), dmax)
         + gammaHA*betaHA*vQ0(np.abs(1 - grid), dmax))

    # Enhanced antibody production
    ktilde_adj = ktilde*betaAb

    # Numpy array of values for B, Ab in the form (B0, Ab0, B1, Ab1...).
    # The initial condition corresponds to # B(x, 0) = 0.0, Ab(x, 0) = 1.0
    y0 = np.zeros(2*Nx)
    y0[1 :: 2] = np.ones(Nx)

    # We integrate the PDEs with the new value of k_Ab
    sol = odeint(affinityMaturation, y0, t, args=(t_boost, H, baseQ, Q,
                 ktilde_adj, mu, dx), ml=2, mu=2)

    # Numpy array with the errors squared
    eSq = np.array([(np.sum(sol[np.argwhere(t == x[i])[0][0]][1 :: 2])*dx
                   - y[i])**2 for i in range(len(x))])

    return np.mean(eSq)**0.5


def rmse_log(parameters, *data):

    """
    Integrates the PDEs and computes the root mean square of erros between the
    model and the data for a given a set of parameter values, in log scale.
    Non-adjuvanted case.
    Parameters:
    gammaNA and gammaHA are the immunigenicities of NA and HA, respectively.
    mu is the activity level decay rate in 1/days.
    dmax is the affinity cutoff.
    """

    # We take the parameter values and data from the arguments.
    gammaNA, gammaHA, mu, dmax = parameters
    x, y = data

    # We create the naive B cell 'pool'.
    H = Htilde*0.5*(np.sign(grid - 0.99*dmax) + np.sign(1.0 - 0.99*dmax - grid))


    # We compute the numpy array of base and reaction affinities for the given
    # immunogenicities using the vectorised version of Q0.
    baseQ = vQ0(np.abs(grid), dmax) + vQ0(np.abs(1 - grid), dmax)
    Q = gammaNA*vQ0(np.abs(grid), dmax) + gammaHA*vQ0(np.abs(1 - grid), dmax)

    # Numpy array of values for B, Ab in the form (B0, Ab0, B1, Ab1...).
    # The initial condition corresponds to # B(x, 0) = 0.0, Ab(x, 0) = 1.0
    y0 = np.zeros(2*Nx)
    y0[1 :: 2] = np.ones(Nx)

    # We integrate the PDEs
    sol = odeint(affinityMaturation, y0, t, args=(t_boost, H, baseQ, Q, ktilde,
                 mu, dx), ml=2, mu=2)

    # Numpy array with the errors squared
    eSq = np.array([(np.log2(np.sum(sol[np.argwhere(t == x[i])[0][0]][1 :: 2])
                   *dx) - np.log2(y[i]))**2 for i in range(len(x))])

    return np.mean(eSq)**0.5


def rmse_adj_log(parameters, *data):

    """
    Integrates the PDEs and computes the root mean square of erros between the
    model and the data for a given a set of parameter values, in log2 scale.
    Adjuvanted case: enhanced immunogenicities, cross-reactivity and antibody
    production.
    Parameters:
    betaNA and betaHA are the enhancement factors for the immunogenicities,
    with gamma -> gamma*beta for both NA and HA, and betaNA > betaHA
    reflecting enhanced cross-reactivity.
    betaAb is the enhancement factor for the production of antibodies, with
    k_Ab -> k_Ab*betaAb.
    """

    # We take the parameter values and data from the arguments.
    betaNA, betaHA, betaAb = parameters
    x, y, base_params = data
    gammaNA, gammaHA, mu, dmax, baseQ, H = base_params

    # We compute the numpy array of reaction affinities for the given
    # immunogenicities and adjuvanticity parameters using the vectorised version
    # of Q0.
    Q = (gammaNA*betaNA*vQ0(np.abs(grid), dmax)
         + gammaHA*betaHA*vQ0(np.abs(1 - grid), dmax))

    # Enhanced antibody production
    ktilde_adj = ktilde*betaAb

    # Numpy array of values for B, Ab in the form (B0, Ab0, B1, Ab1...).
    # The initial condition corresponds to # B(x, 0) = 0.0, Ab(x, 0) = 1.0
    y0 = np.zeros(2*Nx)
    y0[1 :: 2] = np.ones(Nx)

    # We integrate the PDEs with the new value of k_Ab
    sol = odeint(affinityMaturation, y0, t, args=(t_boost, H, baseQ, Q,
                 ktilde_adj, mu, dx), ml=2, mu=2)

    # Numpy array with the errors squared
    eSq = np.array([(np.log2(np.sum(sol[np.argwhere(t == x[i])[0][0]][1 :: 2])
                   *dx) - np.log2(y[i]))**2 for i in range(len(x))])

    return np.mean(eSq)**0.5
