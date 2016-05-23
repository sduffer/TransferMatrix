# -*- coding: utf-8 -*-
"""
@Scott Dufferwiel, University of Sheffield 2016

This script contains the functions which produce the arrays of refractive indices
and thicknesses of a microcavity required for the transfer matrix program
"""
import numpy as np
import scipy as sp

def top_DBR(n1,n2,Npairs,wl):
    d1 = wl/(4.*n1)
    d2 = wl/(4.*n2)

    pair_n = (n1,n2)
    pair_d = (d1,d2)
    top_n = pair_n
    top_d = pair_d

    for x in range(0, Npairs-1):
        a = np.append(top_n, pair_n)
        b = np.append(top_d, pair_d)
        top_n = a
        top_d = b

    return (top_n,top_d)

def cavity(n,q,wl):
    cavity_n = n
    cavity_d = (q*wl)/(2.*n)
    return (cavity_n, cavity_d)

def bottom_DBR(n1,n2,Npairs,wl):
    d1 = wl/(4.*n1)
    d2 = wl/(4.*n2)

    pair_n = (n1,n2)
    pair_d = (d1,d2)
    bot_n = pair_n
    bot_d = pair_d

    for x in range(0, Npairs-1):
        a = np.append(bot_n,pair_n)
        b = np.append(bot_d,pair_d)
        bot_n = a
        bot_d = b

    return (bot_n,bot_d)
