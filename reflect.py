# -*- coding: utf-8 -*-
"""
Scott Dufferwiel

These functions calculate the reflectivity or transmission
of a DBR based microcavity using the transfer matrix method

"""
import numpy as np
import scipy as sp

def snell(n1,n2,theta1):
    theta1 = theta1*np.pi/180
    return sp.arcsin(np.real_if_close(n1*np.sin(theta1) / n2))

def interface_r(polarisation, n_i, n_f, th_i, th_f):
    """
    reflection amplitude at interface

    polarisation is s/p

    n_i, n_f are refractive indicies for incident and final media

    th_i, th_f are angles (in radians) for incident and final fields
    """
    if polarisation == 's':
        return ((n_i * np.cos(th_i) - n_f * np.cos(th_f)) / (n_i * np.cos(th_i) + n_f * np.cos(th_f)))
    elif polarisation == 'p':
        return ((n_f * np.cos(th_i) - n_i * np.cos(th_f)) / (n_f * np.cos(th_i) + n_i * np.cos(th_f)))
    else:
        raise ValueError("Polarisation value must be 's' or 'p'")

def interface_t(polarisation, n_i, n_f, th_i, th_f):
    """
    transmission amplitude (from Fresnel equations)

    polarisation is s/p

    n_i, n_f are refractive indicies for incident and final media

    th_i, th_f are propagation angle for incident and final fields

    """
    if polarisation == 's':
        return 2 * n_i * np.cos(th_i) / (n_i * np.cos(th_i) + n_f * np.cos(th_f))
    elif polarisation == 'p':
        return 2 * n_i * np.cos(th_i) / (n_f * np.cos(th_i) + n_i * np.cos(th_f))
    else:
        raise ValueError("Polarisation must be 's' or 'p'")

"""
To calculate reflectivity we loop through structure and multiply the matching matrices
i.e. match air to A -> propagate across A -> match A to B --> propagate across B ->
match B to substrate
"""

def tmm_reflect(polarisation, n, d, theta, wl, a):
    """
    polarisation is either 's' or 'p'
    n_list is the list of refractive indices across structure
    d_list is the list of layer thicknesses
    theta is the incident angle
    wl is the wavelength in which the reflectivity is being calculated
    """
    nair = 1                # assuming air on above top DBR
    n_substrate = 1.5

    T = np.identity(2)     # create identity matrices for s and p polarisations
    theta_i = np.pi/180 * theta # convert initial incident angle to radians

    E_field = np.matrix([[1],[0]])
    z = np.arange(0,np.sum(d),0.01)
    E_out = np.zeros(0)

    for x in range(0, np.size(n)):
        zr = 0
        zl = 0

        nl = n[x]           # set index to left of interface
        if x == np.size(n)-1:
            nr = n_substrate    # if final layer set index to right to substrate index
        else:
            nr = n[x+1]         # otherwise set index to layer to right of interface

        theta_f = snell(nl,nr,theta_i)  # calculate theta in material to right of interface

        kzr = 2*np.pi*nr*np.cos(theta_f)/wl  # calculates wavevector to right of interface
        kzl = 2*np.pi*nl*np.cos(theta_i)/wl # calculates wavevvector to the left of interface

        # Calculate reflectivity and transmission amplitudes at interface
        r = interface_r(polarisation, nl, nr, theta_i, theta_f)
        t = interface_t(polarisation, nl, nr, theta_i, theta_f)

        # match matrix for interface
        M = (1/t)*np.matrix([[1,r],[r,1]])

        # propagation matrix
        if x == np.size(n)-1:
            P = np.matrix([[np.exp(-1*np.sqrt(-1+0j)*kzr*0), 0],[0, np.exp(1*np.sqrt(-1+0j)*kzr*0)]])
        else:
            P = np.matrix([[np.exp(-1*np.sqrt(-1+0j)*kzr*d[x+1]), 0],[0, np.exp(1*np.sqrt(-1+0j)*kzr*d[x+1])]])

        T = T*M*P

    theta_i = theta_f

    r = T[1,0]/T[0,0]
    t = 1/T[0,0]

    R = r*np.conj(r)
    if polarisation == 's':
        T = t * np.conj(t) * (n_substrate * np.cos(theta) / nair * np.cos(theta))
    elif polarisation == 'p':
        T = t * np.conj(t) *  (n_substrate * np.conj(np.cos(theta)) / nair * np.conj(np.cos(theta)))

    return {'R': R.real, 'T': T.real, 'E': E_out, 'z' : z}
