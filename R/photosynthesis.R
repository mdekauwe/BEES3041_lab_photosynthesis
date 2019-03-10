####
##  Model to calculate net leaf-level C3 photosynthesis.
##
##  Rate of photosynthesis in a leaf is calculated as the minimum of two states,
##  ignoring a third state called the triose phosphate (TPU) or "export"
##  limitation. This condition occurs under high levels of irradiance but is
##  thought to rarely occur in field conditions.
##
##  In one state, the rate of photosynthesis is predicted by the
##  properties of Ribulose-1,5-bisphosphate carboxylase/oxygenase (Rubisco),
##  assuming a saturating supply of substrate, RuBP. This state is called
##  Rubisco-limited photosynthesis and typically occurs when [CO2] is low. The
##  The model assumes that Rubisco-limited A follows a Michaelis–Menten
##  response function modified to account for a competitive inhibitor, oxygen.
##  As is typical of any Michaelis–Menten reaction, increasing the limiting
##  substrate (CO2), the amount of enzyme present (Rubisco), or decreasing
##  the competitive inhibitor (O2) will yield higher reaction rates.
##
##  In the other state, photosynthetic rates are predicted assuming that the
##  rate of regeneration of RuBP is limiting and so RuBP is used at a constant
##  rate; this is called RuBP-regeneration-limited photosynthesis. This
##  typically occurs at higher [CO2]. RuBP-regeneration-limited photosynthesis
##  includes the conditions where light intensity limits the rate of
##  photosynthesis. RuBP-regeneration-limited photosynthesis increases as [CO2]
##  increases because increasing [CO2] causes more RuBP to be carboxylated at
##  the expense of oxygenation.
##
##  References:
##  ----------
##  * Farquhar, G.D., Caemmerer, S. V. and Berry, J. A. (1980) A biochemical
##    model of photosynthetic CO2 assimilation in leaves of C3 species. Planta,
##    149, 78-90.
##  * Medlyn, B. E., Dreyer, E., Ellsworth, D., Forstreuter, M., Harley, P.C.,
##    Kirschbaum, M.U.F., Leroux, X., Montpied, P., Strassemeyer, J.,
##    Walcroft, A., Wang, K. and Loustau, D. (2002) Temperature response of
##    parameters of a biochemically based model of photosynthesis. II.
##    A review of experimental data. Plant, Cell and Enviroment 25, 1167-1179.
##
##
##  author: Martin De Kauwe
##  date: 10th March, 2019
##  email: mdekauwe@gmail.com
####

source("utils.R")

calc_photosynthesis <-function(p, Tleaf) {
  #
  #
  #
  #


  # calculate temp dependancies of MichaelisMenten constants for CO2, O2
  Km <- calc_michaelis_menten_constants(p, Tleaf)

  print(Km)

  An <- 0.0
  
  return (An)
}

calc_michaelis_menten_constants <- function(p, Tleaf) {
  #
  # Michaelis-Menten constant for O2/CO2, Arrhenius temp dependancy
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
  #     Oi : float
  #     åintercellular concentration of O2 [mmol mol-1]
  #
  #   Returns:
  #   --------
  #     Km : float
  #       Michaelis-Menten constant
  #

  Kc <- arrh(p$Kc25, p$Ec, Tleaf)
  Ko <- arrh(p$Ko25, p$Eo, Tleaf)

  Km <- Kc * (1.0 + p$Oi / Ko)

  return ( Km )
}

calc_electron_transport_rate <- function(p, Par, Jmax) {
  #
  # Electron transport rate for a given absorbed irradiance
  #
  #   Args:
  #   -----
  #   p : struct
  #     contains all the model params
  #   Par : float
  #     photosynthetically active radiation [umol m-2 s-1].
  #   Jmax : float
  #     potential rate of electron transportzw
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

  A <- p["theta_J"]
  B <- -(p["alpha"] * Par + Jmax);
  C <- p["alpha"] * Par * Jmax;

  J <- quadratic(A, B, C, large=False)

  return ( J )
}

arrh <- function(k25, Ea, Tk) {
  #
  # Temperature dependence of kinetic parameters is described by an
  # Arrhenius function.
  #
  #   Args:
  #   -----
  #   k25 : float
  #     rate parameter value at 25 degC or 298 K
  #   Ea : float
  #     activation energy for the parameter [J mol-1]
  #   Tk : float
  #     leaf temperature [deg K]
  #
  #   Returns:
  #   -------
  #   kt : float
  #     temperature dependence on parameter
  #
  #   References:
  #   -----------
  #   * Medlyn et al. 2002, PCE, 25, 1167-1179.
  #
  #
  return ( k25 * exp((Ea * (Tk - 298.15)) / (298.15 * c.RGAS * Tk)) )
}

peaked_arrh <- function(k25, Ea, Tk, deltaS, Hd) {
  #
  # Temperature dependancy approximated by peaked Arrhenius eqn,
  # accounting for the rate of inhibition at higher temperatures.
  #
  #   Args:
  #   -----
  #   k25 : float
  #     rate parameter value at 25 degC or 298 K
  #   Ea : float
  #      activation energy for the parameter [J mol-1]
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
  #     temperature dependence on parameter
  #
  #   References:
  #   -----------
  #   * Medlyn et al. 2002, PCE, 25, 1167-1179.
  #

  arg1 <- arrh(k25, Ea, Tk)
  arg2 <- 1.0 + exp((298.15 * deltaS - Hd) / (298.15 * c.RGAS))
  arg3 <- 1.0 + exp((Tk * deltaS - Hd) / (Tk * c.RGAS))

  return ( arg1 * arg2 / arg3 )
}

assim <- function(Ci, gamma_star, a1, a2) {
  #
  # Calculation of assimilation rate with the limitation defined by the
  # variables passed as a1 and a2, i.e. if we are calculating vcmax or
  # jmax limited assimilation rates.
  #
  #   Args:
  #   -----
  #   Ci : float
  #     intercellular CO2 concentration.
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
