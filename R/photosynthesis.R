####
##  Model to calculate net leaf-level C3 photosynthesis.
##
##  author: Martin De Kauwe
##  date: 13th March, 2019
##  email: mdekauwe@gmail.com
####

source("R/utils.R")
source("R/constants.R")

calc_photosynthesis <-function(p, Tleaf, PAR, Cs, vpd, peaked_Vcmax=TRUE,
                               peaked_Jmax=TRUE) {
  #
  #   Calculate photosyntheis following the Farquhar, von Caemmerer & Berry
  #   (1980) model of C3 photosynthesis.
  #
  #   The model mechanistically represents the effects of PAR, temperature and
  #   [CO2] on photosynthesis. The rate of photosynthesis in a leaf is
  #   calculated as the minimum of two states, ignoring a third state called
  #   the triose phosphate (TPU) or "export" limitation. This condition occurs
  #   under high levels of irradiance but is thought to rarely occur in field
  #   conditions.
  #
  #   In one state, the rate of photosynthesis is predicted by the
  #   properties of Ribulose-1,5-bisphosphate carboxylase/oxygenase (Rubisco),
  #   assuming a saturating supply of substrate, RuBP. This state is called
  #   Rubisco-limited photosynthesis and typically occurs when [CO2] is low.
  #   The model assumes that Rubisco-limited A follows a Michaelis–Menten
  #   response function modified to account for a competitive inhibitor, oxygen.
  #   As is typical of any Michaelis–Menten reaction, increasing the limiting
  #   substrate (CO2), the amount of enzyme present (Rubisco), or decreasing
  #   the competitive inhibitor (O2) will yield higher reaction rates.
  #
  #   In the other state, photosynthetic rates are predicted assuming that the
  #   rate of regeneration of RuBP is limiting and so RuBP is used at a constant
  #   rate; this is called RuBP-regeneration-limited photosynthesis. This
  #   typically occurs at higher [CO2]. RuBP-regeneration-limited
  #   photosynthesis includes the conditions where light intensity limits the
  #   rate of photosynthesis. RuBP-regeneration-limited photosynthesis
  #   increases as [CO2] increases because increasing [CO2] causes more RuBP
  #   to be carboxylated at the expense of oxygenation.
  #
  #   The model has two major parameters, the potential rate of electron
  #   transport (Jmax) and the maximum rate of Rubisco activity (Vcmax).
  #   Following Farquhar et al., photosynthesis is solved as the minimum of two
  #   limiting rates (Ac and Aj), assuming the conductance between the
  #   intercellular space and the site of carboxylation is negligible, i.e. an
  #   infinite mesophyll conductance.
  #
  #   Args:
  #   -----
  #   p : list
  #      contains all the model parameters
  #   Tleaf : float
  #      leaf temperature [deg K]
  #   PAR : float
  #      photosynthetically active radiation [umol m-2 s-1].
  #   Cs : float
  #     CO2 concentration at the leaf surface [umol mol-1]
  #   vpd : float
  #     vapour pressure deficit at the leaf surface [kPa]
  #   peaked_Vcmax : logical
  #     Use the peaked Arrhenius function (if true)
  #   peaked_Jmax : logical
  #     Use the peaked Arrhenius function (if true)
  #
  #   Returns:
  #   --------
  #   An : float
  #     Net leaf assimilation rate [umol m-2 s-1]
  #   gsw:  float
  #     stomatal conductance to water [mol m-2 s-1]
  #
  #   References:
  #   ----------
  #   * Farquhar, G.D., Caemmerer, S. V. and Berry, J. A. (1980) A biochemical
  #     model of photosynthetic CO2 assimilation in leaves of C3 species.
  #     Planta, 149, 78-90.
  #   * Medlyn, B. E., Dreyer, E., Ellsworth, D., Forstreuter, M., Harley, P.C.,
  #     Kirschbaum, M.U.F., Leroux, X., Montpied, P., Strassemeyer, J.,
  #     Walcroft, A., Wang, K. and Loustau, D. (2002) Temperature response of
  #     Parameters of a biochemically based model of photosynthesis. II.
  #     A review of experimental data. Plant, Cell and Enviroment 25, 1167-1179.
  #
  
  # g1 and g0 are in units of H20, g0 must be converted to CO2 (but not g1, see below)
  if (any((is_close(p$g0, 0.0)))) {
    # Numerical fix for a zero g0
    p.g0 <- 1E-09
  } else {
    p.g0 = p.g0 * GSW_2_GSC
  }
  
  # calculate temp dependancies of Michaelis-Menten constants for CO2, O2
  Km <- calc_michaelis_menten_constants(p, Tleaf)

  # CO2 compensation point in the absence of mitochondrial respiration
  # [umol mol-1]
  gamma_star <- arrh(p$gamstar25, p$Eag, Tleaf)
 
  # Calculate the maximum rate of Rubisco activity (Vcmax), accounting for
  # temperature dependancies
  if (peaked_Vcmax) {
    Vcmax <- peaked_arrh(p$Vcmax25, p$Eav, Tleaf, p$deltaSv, p$Hdv)
  } else {
    Vcmax <- arrh(p$Vcmax25, p$Eav, Tleaf)
  }

  # Calculate the potential rate of electron transport (Jmax), accounting for
  # temperature dependancies
  if (peaked_Jmax) {
    Jmax <- peaked_arrh(p$Jmax25, p$Eaj, Tleaf, p$deltaSj, p$Hdj)
  } else {
    Jmax <- arrh(p$Jmax25, p$Eaj, Tleaf)
  }

  # Leaf mitochondrial respiration in the light or day respiration
  # (umol m-2 s-1). Following Collatz et al. (1991), assume Rd 1.5% of Vcmax
  Rd <- 0.015 * Vcmax
  
  # Rate of electron transport, which is a function of absorbed PAR
  J <- calc_electron_transport_rate(p, PAR, Jmax)
  Vj <- J / 4.0

  gs_over_a <- calc_stomatal_coeff(p, Cs, vpd)
 
  # Catch for low PAR, issue with Vj
  if ( any(is_close(PAR, 0.0) | is_close(Vj, 0.0)) ) {
    Cic <- Cs
    Cij <- Cs
  } else {
    # Solution when Rubisco activity is limiting
    Cic <- solve_ci(p, gs_over_a, Rd, Cs, gamma_star, Vcmax, Km)
   
    # Solution when electron transport rate is limiting
    Cij <- solve_ci(p, gs_over_a, Rd, Cs, gamma_star, Vj, 2.0*gamma_star)
  }
  
  # Catch for negative Ci and instances where Ci > Cs
  if ( any((Cic <= 0.0) | (Cic > Cs)) ) {
    # Rate of photosynthesis when Rubisco activity is limiting
    Ac <- 0.0
  } else {
    # Rate of photosynthesis when Rubisco activity is limiting
    Ac <- assim(Cic, gamma_star, Vcmax, Km)

    # Rate of photosynthesis when RuBP regeneration is limiting
    Aj <- assim(Cij, gamma_star, Vj, 2.0*gamma_star)
  }

  # When below light-compensation points, assume Ci=Ca.
  if ( any(Aj <= Rd + 1E-09) ) {
    Cij <- Cs

    # Rate of photosynthesis when RuBP regeneration is limiting
    Aj <- assim(Cij, gamma_star, Vj, 2.0*gamma_star)
  }
  
  # Hyperbolic minimum of Ac and Aj to smooth over discontinuity when moving
  # from electron # transport limited to rubisco limited photosynthesis
  A <- -mapply(quadratic, 1.0-1E-04, Ac+Aj, Ac*Aj, large=TRUE)

  # Net photosynthesis rate (umol m-2 s-1)
  An <- A - Rd

  # Calculate conductance to CO2 (mol m-2 s-1)
  gsc <- max(p$g0, p$g0 + gs_over_a * An)

  # Calculate conductance to water (mol m-2 s-1)
  gsw <- gsc * GSC_2_GSW

  return ( list(An=An, Ac=Ac, Aj=Aj, gsc=gsc, Vcmax=Vcmax, Cic=Cic) )

}

calc_michaelis_menten_constants <- function(p, Tleaf) {
  #
  #   Michaelis-Menten constant for O2/CO2, Arrhenius temp dependancy
  #
  #   Args:
  #   -----
  #   Tleaf : float
  #     leaf temperature [deg K]
  #   Kc25 : float
  #     Michaelis-Menten coefficents for carboxylation by Rubisco at
  #     25degC [umol mol-1] or 298 K
  #   Kc25 : float
  #     Michaelis-Menten coefficents for oxygenation by Rubisco at
  #     25degC [mmol mol-1]. Note value in Bernacchie 2001 is in mmol!!
  #     or 298 K
  #   Ec : float
  #     Activation energy for carboxylation [J mol-1]
  #   Eo : float
  #     Activation energy for oxygenation [J mol-1]
  #   Oi : float
  #     The intercellular concentration of O2 [mmol mol-1]
  #
  #   Returns:
  #   --------
  #     Km : float
  #       Michaelis-Menten constant [umol mol−1]
  #

  # Michaelis– Menten coefficients of Rubisco, Kc (umol mol−1)
  Kc <- arrh(p$Kc25, p$Ec, Tleaf)

  # Michaelis– Menten coefficients of oxygenation, Ko (mmol mol−1)
  Ko <- arrh(p$Ko25, p$Eo, Tleaf)

  # The Michaelis–Menten constant (umol mol−1)
  Km <- Kc * (1.0 + p$Oi / Ko)

  return ( Km )
}

arrh <- function(k25, Ea, Tk) {
  #
  #   Temperature dependence of kinetic PARameters is described by an
  #   Arrhenius function.
  #
  #   Args:
  #   -----
  #   k25 : float
  #     rate PARameter value at 25 degC or 298 K
  #   Ea : float
  #     activation energy for the PARameter [J mol-1]
  #   Tk : float
  #     leaf temperature [deg K]
  #
  #   Returns:
  #   -------
  #   kt : float
  #     temperature dependence on PARameter
  #
  #   References:
  #   -----------
  #   * Medlyn et al. 2002, PCE, 25, 1167-1179.
  #
  #
  return ( k25 * exp((Ea * (Tk - 298.15)) / (298.15 * RGAS * Tk)) )
}

peaked_arrh <- function(k25, Ea, Tk, deltaS, Hd) {
  #
  #   Temperature dependancy approximated by peaked Arrhenius eqn,
  #   accounting for the rate of inhibition at higher temperatures.
  #
  #   Args:
  #   -----
  #   k25 : float
  #     rate PARameter value at 25 degC or 298 K
  #   Ea : float
  #      activation energy for the PARameter [J mol-1]
  #   Tk : float
  #      leaf temperature [deg K]
  #   deltaS : float
  #     entropy factor [J mol-1 K-1)
  #   Hd : float
  #     describes rate of decrease about the optimum temp [J mol-1]
  #
  #   Returns:
  #   -------
  #   kt : float
  #     temperature dependence on PARameter
  #
  #   References:
  #   -----------
  #   * Medlyn et al. 2002, PCE, 25, 1167-1179.
  #

  arg1 <- arrh(k25, Ea, Tk)
  arg2 <- 1.0 + exp((298.15 * deltaS - Hd) / (298.15 * RGAS))
  arg3 <- 1.0 + exp((Tk * deltaS - Hd) / (Tk * RGAS))

  return ( arg1 * arg2 / arg3 )
}


assim <- function(Ci, gamma_star, a1, a2) {
  #
  #   Calculation of assimilation rate with the limitation defined by the
  #   variables passed as a1 and a2, i.e. if we are calculating vcmax or
  #   jmax limited assimilation rates.
  #
  #   Args:
  #   -----
  #   Ci : float
  #     intercellular CO2 concentration [umol mol-1]
  #   gamma_star : float
  #     CO2 compensation point in the abscence of mitochondrial respiration
  #     [umol m-2 s-1]
  #   a1 : float
  #     variable depends on whether the calculation is light or rubisco
  #     limited.
  #   a2 : float
  #     variable depends on whether the calculation is light or rubisco
  #     limited.
  #
  #   Returns:
  #   -------
  #   assimilation_rate : float
  #     assimilation rate assuming either light or rubisco limitation.
  #     [umol m-2 s-1]
  #

  return ( a1 * (Ci - gamma_star) / (a2 + Ci) )
}

calc_electron_transport_rate <-function(p, PAR, Jmax) {
  #
  #   Electron transport rate for a given absorbed irradiance
  #
  #   Args:
  #   -----
  #   p : list
  #     contains all the model parameters
  #   PAR : float
  #     photosynthetically active radiation [umol m-2 s-1].
  #   Jmax : float
  #     potential rate of electron transport  [umol m-2 s-1].
  #   theta_J : float
  #     Curvature of the light response (-)
  #   alpha : float
  #     Leaf quantum yield (initial slope of the A-light response curve)
  #     [mol mol-1]
  #
  #   Reference:
  #   ----------
  #   * Farquhar G.D. & Wong S.C. (1984) An empirical model of stomatal
  #     conductance. Australian Journal of Plant Physiology 11, 191-210,
  #     eqn A but probably clearer in:
  #   * Leuning, R. et a., Leaf nitrogen, photosynthesis, conductance and
  #     transpiration: scaling from leaves to canopies, Plant Cell Environ.,
  #     18, 1183– 1200, 1995. Leuning 1995, eqn C3.
  #

  A <- p$theta_J
  B <- -(p$alpha * PAR + Jmax)
  C <- p$alpha * PAR * Jmax

  J <- mapply(quadratic, A, B, C, large=FALSE)

  return ( J )
}

calc_stomatal_coeff <- function(p, Cs, vpd) {
  #
  #   Stomatal coefficent, hardwired for Medlyn gs model
  #
  #   Args:
  #   -----
  #   p : list
  #     contains all the model parameters
  #   Cs : float
  #     CO2 concentration at the leaf surface [umol mol-1]
  #   vpd : float
  #     vapour pressure deficit at the leaf surface [kPa]
  #   g1 : float
  #     stomatal slope parameter [kpa^0.5]

  # 1.6 (from corrigendum to Medlyn et al 2011) is missing here,
  # because we are calculating conductance to CO2!
  #
  #   References:
  #   -----------
  #   * Medlyn et al. 2002, PCE, 25, 1167-1179.
  #
  
  # Medlyn moment can't have v.low VPD vals
  vpd <- ifelse(vpd<0.05, 0.05, vpd)
  
  # 1.6 (from corrigendum to Medlyn et al 2011) is missing here,
  # because we are calculating conductance to CO2!
  if (any((is_close(Cs, 0.0)))) {
    gs_over_a = 0.0
  } else {
    gs_over_a <- (1.0 + p$g1 / sqrt(vpd)) / Cs
  }

  return ( gs_over_a )
}

solve_ci <- function(p, gs_over_a, Rd, Cs, gamma_star, gamma, beta) {
  #
  #   Solve intercellular CO2 concentration using quadric equation, following
  #   Leuning 1990, see eqn 15a-c, solving simultaneous solution for Eqs 2, 12
  #   and 13
  #
  #   Args:
  #   -----
  #   p : list
  #     contains all the model parameters
  #   gs_over_a : float
  #     stomatal coefficient
  #   Rd : float
  #     day respiration rate [umol m-2 s-1]
  #   Cs : float
  #     leaf surface CO2 concentration [umol mol-1]
  #   gamma_star : float
  #     CO2 compensation point - base rate at 25 deg C / 298 K [umol mol-1]
  #   gamma : float
  #     if calculating Cic, this will be Vcmax
  #     if calculating Cij, this will be Vj
  #   beta : float
  #     if calculating Cic, this will be Km
  #     if calculating Cij, this will be 2.0*gamma_star
  #
  #   Reference:
  #   ----------
  #   * Leuning (1990) Modelling Stomatal Behaviour and Photosynthesis of
  #     Eucalyptus grandis. Aust. J. Plant Physiol., 17, 159-75.
  #

  A = p$g0 + gs_over_a * (gamma - Rd)

  arg1 <- (1. - Cs * gs_over_a) * (gamma - Rd)
  arg2 <- p$g0 * (beta - Cs)
  arg3 <- gs_over_a * (gamma * gamma_star + beta * Rd)
  B <- arg1 + arg2 - arg3

  arg1 <- -(1.0 - Cs * gs_over_a)
  arg2 <- (gamma * gamma_star + beta * Rd)
  arg3 <- p$g0 * beta * Cs
  C <- arg1 * arg2 - arg3

  Ci <- mapply(quadratic, A, B, C, large=TRUE)
  
  return ( Ci )

}
