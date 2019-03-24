#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
LUE model, vaguely following MODIS algorithm

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (24.03.2017)"
__email__ = "mdekauwe@gmail.com"

import constants as c
from utils import quadratic

def calc_lue(fPAR, PAR, T_scal, D_scal, emax=1.23):
    """
    Simple light use efficiency model to calculate GPP

    Monteith (1972) suggested that the NPP of well-watered and fertilized annual
    crop plants was linearly related to the amount of solar energy the plants
    absorbed over a growing season. This logic combined the meteorological
    constraint of available sunlight reaching a site with the ecological
    constraint of the amount of leaf area absorbing the solar energy,
    while avoiding many complexities of canopy micrometeorology and
    carbon balance theory.

    Parameters:
    ----------
    p : list
        contains all the model parameters
    fPAR : float
        fraction of photosynthetically active radiation [0-1].
    PAR : float
        photosynthetically active radiation [MJ m-2 d-1].
    T_scal : float
        Temperature scalar [0-1]
    D_scal : float
        Vapour pressure deficit scalar [0-1]
    emax : float
        light use efficiency constant, g C m−2 MJ−1 APAR

    Returns:
    --------
    GPP : float
        Gross primary productivity [g C m-2 d-1]

    References:
    ----------
    * Monteith, J. (1972) Solar radiation and productivity in tropical ecosystems.
      Journal of Applied Ecology, 9(3), 747-766.
    """

    e = emax * T_scal * D_scal
    GPP = fPAR * PAR * e

    return ( GPP )
