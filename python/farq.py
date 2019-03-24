#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Model to calculate net leaf-level C3 photosynthesis.

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (24.03.2017)"
__email__ = "mdekauwe@gmail.com"

import sys
import numpy as np
import os
import math
import constants as c
from utils import quadratic

class FarquharC3(object):
    """
    Calculate photosyntheis following the Farquhar, von Caemmerer & Berry
    (1980) model of C3 photosynthesis.

    The model mechanistically represents the effects of PAR, temperature and
    [CO2] on photosynthesis. The rate of photosynthesis in a leaf is
    calculated as the minimum of two states, ignoring a third state called
    the triose phosphate (TPU) or "export" limitation. This condition occurs
    under high levels of irradiance but is thought to rarely occur in field
    conditions.

    In one state, the rate of photosynthesis is predicted by the
    properties of Ribulose-1,5-bisphosphate carboxylase/oxygenase (Rubisco),
    assuming a saturating supply of substrate, RuBP. This state is called
    Rubisco-limited photosynthesis and typically occurs when [CO2] is low.
    The model assumes that Rubisco-limited A follows a Michaelis–Menten
    response function modified to account for a competitive inhibitor, oxygen.
    As is typical of any Michaelis–Menten reaction, increasing the limiting
    substrate (CO2), the amount of enzyme present (Rubisco), or decreasing
    the competitive inhibitor (O2) will yield higher reaction rates.

    In the other state, photosynthetic rates are predicted assuming that the
    rate of regeneration of RuBP is limiting and so RuBP is used at a constant
    rate; this is called RuBP-regeneration-limited photosynthesis. This
    typically occurs at higher [CO2]. RuBP-regeneration-limited
    photosynthesis includes the conditions where light intensity limits the
    rate of photosynthesis. RuBP-regeneration-limited photosynthesis
    increases as [CO2] increases because increasing [CO2] causes more RuBP
    to be carboxylated at the expense of oxygenation.

    The model has two major parameters, the potential rate of electron
    transport (Jmax) and the maximum rate of Rubisco activity (Vcmax).
    Following Farquhar et al., photosynthesis is solved as the minimum of two
    limiting rates (Ac and Aj), assuming the conductance between the
    intercellular space and the site of carboxylation is negligible, i.e. an
    infinite mesophyll conductance.

    NB. All calculations in Kelvins...

    References:
    -----------
    * De Pury and Farquhar, G. D. (1997) Simple scaling of photosynthesis from
      leaves to canopies without the errors of big-leaf models. Plant Cell and
      Environment, 20, 537-557.
    * Farquhar, G.D., Caemmerer, S. V. and Berry, J. A. (1980) A biochemical
      model of photosynthetic CO2 assimilation in leaves of C3 species. Planta,
      149, 78-90.
    * Medlyn, B. E., Dreyer, E., Ellsworth, D., Forstreuter, M., Harley, P.C.,
      Kirschbaum, M.U.F., Leroux, X., Montpied, P., Strassemeyer, J.,
      Walcroft, A., Wang, K. and Loustau, D. (2002) Temperature response of
      parameters of a biochemically based model of photosynthesis. II.
      A review of experimental data. Plant, Cell and Enviroment 25, 1167-1179.
    """

    def __init__(self, peaked_Jmax=False, peaked_Vcmax=False,
                 model_Q10=False, gs_model=None):
        """
        Parameters
        ----------
        peaked_Jmax : logical
            use peaked or non-peaked arrhenius func
        peaked_Vcmax : logical
            use peaked or non-peaked arrhenius func
        model_Q10 : logical
            use Q10 to calculate Rd
        gs_model : string
            medlyn or leuning model

        """
        self.peaked_Jmax = peaked_Jmax
        self.peaked_Vcmax = peaked_Vcmax
        self.model_Q10 = model_Q10
        self.gs_model = gs_model

    def photosynthesis(self, p, Cs=None, Tleaf=None, PAR=None, vpd=None,
                       mult=None, scalex=None):
        """
        Parameters
        ----------
        p : struct
            contains all the model params
        Cs : float
            leaf surface CO2 concentration [umol mol-1]
        Tleaf : float
            leaf temp [deg K]
        PAR : float
            photosynthetically active radiation [umol m-2 s-1].
        vpd : float
            vapour pressure deficit [kPa].
        mult : float
            if passing a user defined gs model, i.e. not Medlyn or Leuning.
        scalex : float
            big/two-leaf scalar

        Returns:
        --------
        An : float
            Net leaf assimilation rate [umol m-2 s-1]
        gsc : float
            stomatal conductance to CO2 [mol m-2 s-1]
        """

        # calculate temp dependancies of Michaelis-Menten constants for CO2, O2
        Km = self.calc_michaelis_menten_constants(p, Tleaf)

        # CO2 compensation point in the absence of mitochondrial respiration
        # [umol mol-1]
        gamma_star = self.arrh(p.gamstar25, p.Eag, Tleaf)

        # Calculate the maximum rate of Rubisco activity (Vcmax), accounting
        # for temperature dependancies
        if self.peaked_Vcmax:
            Vcmax = self.peaked_arrh(p.Vcmax25, p.Eav, Tleaf, p.deltaSv, p.Hdv)
        else:
            Vcmax = self.arrh(p.Vcmax25, p.Eav, Tleaf)

        # Calculate the potential rate of electron transport (Jmax), accounting
        # for temperature dependancies
        if self.peaked_Jmax:
            Jmax = self.peaked_arrh(p.Jmax25, p.Eaj, Tleaf, p.deltaSj, p.Hdj)
        else:
            Jmax = self.arrh(p.Jmax25, p.Eaj, Tleaf)

        # Leaf mitochondrial respiration in the light or day respiration
        # (umol m-2 s-1).
        if p.Rd25 is not None:
            Rd = self.calc_resp(Tleaf, p.Q10, p.Rd25, p.Ear)
        else:
            # Following Collatz et al. (1991), assume Rd 1.5% of Vcmax
            Rd = 0.015 * Vcmax

        # Scaling from single leaf to canopy, see Wang & Leuning 1998 appendix C
        if scalex is not None:
            Rd *= scalex
            Vcmax *= scalex
            Jmax *= scalex

        # Rate of electron transport, which is a function of absorbed PAR
        J = self.calc_electron_transport_rate(p, PAR, Jmax)
        Vj = J / 4.0

        (g0, gs_over_a) = self.calc_stomatal_coeff(p, Cs, vpd)

        # Solution when Rubisco activity is limiting
        Cic = self.solve_ci(g0, gs_over_a, Rd, Cs, gamma_star, Vcmax, Km)

        # Solution when electron transport rate is limiting
        Cij = self.solve_ci(g0, gs_over_a, Rd, Cs, gamma_star, Vj,
                              2.0*gamma_star)

        # Catch for low PAR, issue with Vj
        Cic = np.where(np.logical_or(np.isclose(PAR, 0.0), np.isclose(Vj, 0.0)),
                       Cs, Cic)
        Cij = np.where(np.logical_or(np.isclose(PAR, 0.0), np.isclose(Vj, 0.0)),
                       Cs, Cij)

        # Rate of photosynthesis when Rubisco activity is limiting
        Ac = self.assim(Cic, gamma_star, a1=Vcmax, a2=Km)

        # Rate of photosynthesis when RuBP regeneration is limiting
        Aj = self.assim(Cij, gamma_star, a1=Vj, a2=2.0*gamma_star)

        # Catch for negative Ci and instances where Ci > Cs
        Ac = np.where(np.logical_or(Cic <= 0.0, Cic > Cs), 0.0, Ac)
        Aj = np.where(np.logical_or(Cic <= 0.0, Cic > Cs), 0.0, Aj)

        # When below light-compensation points, assume Ci=Ca.
        Aj = np.where(Aj <= Rd + 1E-09,
                      self.assim(Cs, gamma_star, a1=Vj, a2=2.0*gamma_star), Aj)

        # Hyperbolic minimum of Ac and Aj to smooth over discontinuity when
        # moving from electron # transport limited to rubisco limited
        # photosynthesis
        A = -quadratic(a=1.0-1E-04, b=Ac+Aj, c=Ac*Aj, large=True)

        # Net photosynthesis rate (umol m-2 s-1)
        An = A - Rd

        # Calculate conductance to CO2 (mol m-2 s-1)
        gsc = np.maximum(g0, g0 + gs_over_a * An)

        # Calculate conductance to water (mol H2O m-2 s-1)
        gsw = gsc * c.GSC_2_GSW

        # calculate the real Ci
        Ci = np.where(np.logical_and(gsc > 0.0, An > 0.0), Cs - An / gsc, Cs)

        An = np.where(np.isclose(Cs, 0.0), 0.0 - Rd, An)
        gsc = np.where(np.isclose(Cs, 0.0), 0.0, gsc)
        Ci = np.where(np.isclose(Cs, 0.0), Cs, Ci)

        return (An, Ac, Aj, gsw, Rd)

    def calc_michaelis_menten_constants(self, p, Tleaf):
        """ Michaelis-Menten constant for O2/CO2, Arrhenius temp dependancy

        Parameters:
        ----------
        p : struct
            contains all the model params
        Tleaf : float
            leaf temperature [deg K]

        Returns:
        --------
        Km : float
            Michaelis–Menten constant [umol mol-1]
        """
        # Michaelis– Menten coefficients for carboxylation, Kc (umol mol−1)
        Kc = self.arrh(p.Kc25, p.Ec, Tleaf)

        # Michaelis– Menten coefficients for oxygenation, Ko (mmol mol−1)
        Ko = self.arrh(p.Ko25, p.Eo, Tleaf)

        # The Michaelis–Menten constant (umol mol−1)
        Km = Kc * (1.0 + p.Oi / Ko)

        return Km

    def calc_electron_transport_rate(self, p, PAR, Jmax):
        """
        Electron transport rate for a given absorbed irradiance

        Parameters
        ----------
        p : struct
            contains all the model params
        PAR : float
            photosynthetically active radiation [umol m-2 s-1].
        Jmax : float
            potential rate of electron transport

        Reference:
        ----------
        * Farquhar G.D. & Wong S.C. (1984) An empirical model of stomatal
          conductance. Australian Journal of Plant Physiology 11, 191-210,
          eqn A but probably clearer in:
        * Leuning, R. et a., Leaf nitrogen, photosynthesis, conductance and
          transpiration: scaling from leaves to canopies, Plant Cell Environ.,
          18, 1183– 1200, 1995. Leuning 1995, eqn C3.
        """
        A = p.theta_J
        B = -(p.alpha * PAR + Jmax);
        C = p.alpha * PAR * Jmax;

        J = quadratic(a=A, b=B, c=C, large=False)

        return J

    def solve_ci(self, g0, gs_over_a, rd, Cs, gamma_star, gamma, beta):
        """
        Solve intercellular CO2 concentration using quadric equation, following
        Leuning 1990, see eqn 15a-c, solving simultaneous solution for Eqs 2, 12
        and 13

        Parameters
        ----------
        g0 : float
            residual stomatal conductance as net assimilation rate reaches zero
            (mol m-2 s-1)
        gs_over_a : float
            gs / A
        rd : float
            Rspiration rate [umol m-2 s-1]
        Cs : float
            leaf surface CO2 concentration [umol mol-1]
        gamma_star : float
            CO2 compensation point - base rate at 25 deg C / 298 K [umol mol-1]
        gamma : float
            if calculating Cic, this will be Vcmax
            if calculating Cij, this will be Vj
        beta : float
            if calculating Cic, this will be Km
            if calculating Cij, this will be 2.0*gamma_star

        Reference:
        ----------
        * Leuning (1990) Modelling Stomatal Behaviour and Photosynthesis of
          Eucalyptus grandis. Aust. J. Plant Physiol., 17, 159-75.
        """

        A = g0 + gs_over_a * (gamma - rd)

        arg1 = (1. - Cs * gs_over_a) * (gamma - rd)
        arg2 = g0 * (beta - Cs)
        arg3 = gs_over_a * (gamma * gamma_star + beta * rd)
        B = arg1 + arg2 - arg3

        arg1 = -(1.0 - Cs * gs_over_a)
        arg2 = (gamma * gamma_star + beta * rd)
        arg3 =  g0 * beta * Cs
        C = arg1 * arg2 - arg3

        Ci = quadratic(a=A, b=B, c=C, large=True)

        return Ci

    def calc_stomatal_coeff(self, p, Cs, vpd):
        """
        Stomatal coefficent, gs_over_a

        Parameters:
        ----------
        p : list
            contains all the model parameters
        Cs : float
            CO2 concentration at the leaf surface [umol mol-1]
        vpd : float
            vapour pressure deficit at the leaf surface [kPa]
        g1 : float
            stomatal slope parameter [kpa^0.5; Medlyn model, unitless othewise]

        Returns:
        --------
        g0 : float
            residual stomatal conductance as net assimilation rate reaches zero
            (mol m-2 s-1)
        gs_over_a : float
            stomatal coefficient
        """
        if self.gs_model == "leuning":
            if np.isclose(p.g0, 0.0):
                # I want a zero g0, but zero messes up the convergence,
                # numerical fix
                g0 = 1E-09

            # convert to conductance to CO2
            g0 = g0 * c.GSW_2_GSC
            gs_over_a = p.g1 / (Cs - gamma_star) / (1.0 + vpd / p.D0)

            # conductance to CO2
            gs_over_a *= c.GSW_2_GSC
            ci_over_ca = 1.0 - 1.6 * (1.0 + vpd / p.D0) / p.g1

        elif self.gs_model == "medlyn":
            if np.isclose(p.g0, 0.0):
                # I want a zero g0, but zero messes up the convergence,
                # numerical fix
                g0 = 1E-09

            # Medlyn moment can't have v.low VPD vals
            vpd = np.where(vpd < 0.05, 0.05, vpd)

            # 1.6 (from corrigendum to Medlyn et al 2011) is missing here,
            # because we are calculating conductance to CO2!
            gs_over_a = np.where(np.isclose(Cs, 0.0), 0.0,
                                (1.0 + p.g1 / np.sqrt(vpd)) / Cs)

            #ci_over_ca = p.g1 / (p.g1 + np.sqrt(vpd))

        elif self.gs_model == "user_defined":
            # convert to conductance to CO2
            g0 = p.g0 / c.GSC_2_GSW
            gs_over_a = mult / c.GSC_2_GSW

        return (g0, gs_over_a)

    def arrh(self, k25, Ea, Tk):
        """
        Temperature dependence of kinetic parameters is described by an
        Arrhenius function.

        Parameters:
        ----------
        k25 : float
            rate parameter value at 25 degC or 298 K
        Ea : float
            activation energy for the parameter [J mol-1]
        Tk : float
            leaf temperature [deg K]

        Returns:
        -------
        kt : float
            temperature dependence on parameter

        References:
        -----------
        * Medlyn et al. 2002, PCE, 25, 1167-1179.
        """
        return k25 * np.exp((Ea * (Tk - 298.15)) / (298.15 * c.RGAS * Tk))

    def peaked_arrh(self, k25, Ea, Tk, deltaS, Hd):
        """
        Temperature dependancy approximated by peaked Arrhenius eqn,
        accounting for the rate of inhibition at higher temperatures.

        Parameters:
        ----------
        k25 : float
            rate parameter value at 25 degC or 298 K
        Ea : float
            activation energy for the parameter [J mol-1]
        Tk : float
            leaf temperature [deg K]
        deltaS : float
            entropy factor [J mol-1 K-1)
        Hd : float
            describes rate of decrease about the optimum temp [J mol-1]

        Returns:
        -------
        kt : float
            temperature dependence on parameter

        References:
        -----------
        * Medlyn et al. 2002, PCE, 25, 1167-1179.

        """
        arg1 = self.arrh(k25, Ea, Tk)
        arg2 = 1.0 + np.exp((298.15 * deltaS - Hd) / (298.15 * c.RGAS))
        arg3 = 1.0 + np.exp((Tk * deltaS - Hd) / (Tk * c.RGAS))

        return arg1 * arg2 / arg3

    def assim(self, Ci, gamma_star, a1, a2):
        """
        Calculation of photosynthesis with the limitation defined by the
        variables passed as a1 and a2, i.e. if we are calculating vcmax or
        jmax limited assimilation rates.

        Parameters:
        ----------
        Ci : float
            intercellular CO2 concentration.
        gamma_star : float
            CO2 compensation point in the abscence of mitochondrial respiration
            [umol m-2 s-1]
        a1 : float
            variable depends on whether the calculation is light or rubisco
            limited.
        a2 : float
            variable depends on whether the calculation is light or rubisco
            limited.

        Returns:
        -------
        assimilation_rate : float
            assimilation rate assuming either light or rubisco limitation.
            [umol m-2 s-1]
        """
        return a1 * (Ci - gamma_star) / (a2 + Ci)

    def calc_resp(self, Tleaf=None, Q10=None, Rd25=None, Ear=None, Tref=25.0):
        """ Calculate leaf respiration accounting for temperature dependence.

        Parameters:
        ----------
        Tleaf : float
            leaf temp [deg K]
        Q10 : float
            ratio of respiration at a given temperature divided by respiration
            at a temperature 10 degrees lower
        Rd25 : float
            Estimate of respiration rate at the reference temperature 25 deg C
            or or 298 K
        Ear : float
            activation energy for the parameter [J mol-1]
        Tref : float
            reference temperature

        Returns:
        -------
        Rd : float
            leaf respiration

        References:
        -----------
        Tjoelker et al (2001) GCB, 7, 223-230.
        """
        if self.model_Q10:
            Rd = Rd25 * Q10**(((Tleaf - c.DEG_2_KELVIN) - Tref) / 10.0)
        else:
            Rd = self.arrh(Rd25, Ear, Tleaf)

        return Rd

    def adj_for_low_temp(self, param, Tk, lower_bound=0.0, upper_bound=10.0):
        """
        Function allowing Jmax/Vcmax to be forced linearly to zero at low T

        Parameters:
        ----------
        param : float
            value to adjust
        Tk : float
            air temperature (Kelvin)
        """
        Tc = Tk - c.DEG_2_KELVIN

        if Tc < lower_bound:
            param = 0.0
        elif Tc < upper_bound:
            param *= (Tc - lower_bound) / (upper_bound - lower_bound)

        return param
