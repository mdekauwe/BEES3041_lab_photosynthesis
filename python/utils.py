#!/usr/bin/env python
"""
Model parameters

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (24.03.2019)"
__email__ = "mdekauwe@gmail.com"

import numpy as np
import constants as c

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    """
    Check if two values are close ...

    Parameters:
    ----------
    a : float
        value to compare
    b : float
        value to compare
    rel_tol : float
        relative tolerance, max allowed diff btw a and b
    abs_tol : float
        minimum absolute tolerance, must be at least zero

    Returns:
    -------
    boolean : logical
        Return True if the values a and b are close to each other and False
        otherwise.
    """
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def vpd_to_rh(vpd, tair, pressure):
    """
    Convert from VPD to RH.

    Parameters:
    ----------
    tair : float
        air temperature (deg C)
    vpd : float
        Vapour pressure deficit (kPa, needs to be in Pa, see conversion
        below)

    Returns:
    --------
    rh : float
        relative humidity (fraction)

    """
    rh = 1.0 - (vpd * c.KPA_2_PA) / calc_esat(tair)

    return max(0.0, min(1.0, rh))

def calc_esat(tair):
    """
    Saturation vapour pressure or saturation partial pressure of water vapour

    Calculated following Tetens forumula based on Buck.

    Parameters:
    ----------
    tair : float
        air temperature (deg C)

    Returns:
    --------
    esat : float
        Saturation vapor pressure (Pa K-1)

    References:
    * Buck, A. (1981) New equations for computing vapor pressure and
      enhancement factor. Journal of Applied Meteorology, 20, 1527-1532

    but also see...
    * Stull 2000 Meteorology for Scientist and Engineers
    * Jones, H. G. (1992) Plants and microclimate. Second edition, pg 110
      (note error in a - wrong units)

    """
    a = 613.75 # correct units
    b = 17.502
    c = 240.97
    esat = a * np.exp( (b * tair) / (c + tair) )

    # Buck...
    #kpa_2_pa = 1000.
    #a = 0.61121 # kPa
    #b = 17.502
    #c = 240.97 # deg C
    #esat = a * (math.exp(b * tair / (c + tair))) * kpa_2_pa

    return esat

def get_dewpoint(tair, rh):
    """
    The air is saturated when it reaches maximum water holding capacity at a
    given temperature, the dew point

    Formula is apparently relatively accurate for relative humidity values
    above 50%.

    Parameters:
    ----------
    tair : float
        air temperature (deg C)
    RH : float
        relative humidity (percent)

    Returns:
    --------
    Td : float
        Dew point temp (deg C)

    Reference:
    ----------
    * Lawrence, Mark G., 2005: The relationship between relative humidity and
      the dewpoint temperature in moist air: A simple conversion and
      applications. Bull. Amer. Meteor. Soc., 86, 225-233.
      doi: http;//dx.doi.org/10.1175/BAMS-86-2-225
    """
    Td = tair - ((100.0 - rh) / 5.)

    return Td

def quadratic(a=None, b=None, c=None, large=False):
    """ Solve the quadratic equation

    Parameters:
    ----------
    a : float
        co-efficient
    b : float
        co-efficient
    c : float
        co-efficient
    large : logical
        True -> finds largest root
        False -> finds smallest root

    Returns:
    -------
    val : float
        positive root
    """
    d = b**2.0 - 4.0 * a * c # discriminant

    if np.any(d <= 0.0):
        print("imaginary root found")

    if large:
        root = (-b + np.sqrt(d)) / (2.0 * a)
    else:
        root = (-b - np.sqrt(d)) / (2.0 * a)

    root = np.where(np.logical_and(np.isclose(a, 0.0), b > 0.0), -c / b, root)
    root = np.where(np.logical_and(np.isclose(a, 0.0), np.isclose(b, 0.0)),
                    0.0, root)

    return root
