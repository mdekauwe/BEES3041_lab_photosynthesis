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
##    PARameters of a biochemically based model of photosynthesis. II.
##    A review of experimental data. Plant, Cell and Enviroment 25, 1167-1179.
##
##
##  author: Martin De Kauwe
##  date: 10th March, 2019
##  email: mdekauwe@gmail.com
####

source("R/utils.R")
source("R/constants.R")

calc_photosynthesis <-function(p, Tleaf, PAR, Cs, vpd, peaked_Vcmax=TRUE,
                               peaked_Jmax=TRUE) {
  #
  #
  #
  #


  # calculate temp dependancies of MichaelisMenten constants for CO2, O2
  Km <- calc_michaelis_menten_constants(p, Tleaf)

  # Effect of temp on CO2 compensation point
  gamma_star = arrh(p$gamstar25, p$Eag, Tleaf)

  # Calculate temperature dependancies on Vcmax
  if (peaked_Vcmax) {
    Vcmax = peaked_arrh(p$Vcmax25, p$Eav, Tleaf, p$deltaSv, p$Hdv)
  } else {
    Vcmax = arrh(p$Vcmax25, p$Eav, Tleaf)
  }

  # Calculate temperature dependancies on Jmax
  if (peaked_Jmax) {
    Jmax = peaked_arrh(p$Jmax25, p$Eaj, Tleaf, p$deltaSj, p$Hdj)
  } else {
    Jmax = arrh(p$Jmax25, p$Eaj, Tleaf)
  }

  # Leaf mitochondrial respiration in the light or day respiration
  # (umol m-2 s-1)
  Rd = 0.015 * Vcmax

  # Rate of electron transport, which is a function of absorbed PAR
  J <- calc_electron_transport_rate(p, PAR, Jmax)
  Vj <- J / 4.0

  gs_over_a <- calc_stomatal_coeff(p, Cs, vpd)
  print(gs_over_a)

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

arrh <- function(k25, Ea, Tk) {
  #
  # Temperature dependence of kinetic PARameters is described by an
  # Arrhenius function.
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
  # Temperature dependancy approximated by peaked Arrhenius eqn,
  # accounting for the rate of inhibition at higher temperatures.
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

calc_electron_transport_rate <-function(p, PAR, Jmax) {
  #
  # Electron transport rate for a given absorbed irradiance
  #
  #   Args:
  #   -----
  #   p : struct
  #     contains all the model Params
  #   PAR : float
  #     photosynthetically active radiation [umol m-2 s-1].
  #  Jmax : float
  #     potential rate of electron transport
  #  theta_J : float
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

  J <- quadratic(A, B, C, large=FALSE)

  return ( J )
}

calc_stomatal_coeff <- function(p, Cs, vpd) {
  #
  # Stomatal coefficent, hardwired for Medlyn gs model
  #
  #   Args:
  #   -----
  #   p : struct
  #     contains all the model Params
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

  if (is_close(Cs, 0.0)) {
    gs_over_a <- 0.0
  } else {
    gs_over_a <- (1.0 + p$g1 / sqrt(vpd)) / Cs
  }

  return ( gs_over_a )
}

quadratic <- function(a, b, c, large=FALSE) {
  #
  # minimilist quadratic solution as root for J solution should always
  # be positive, so I have excluded other quadratic solution steps. I am
  # only returning the smallest of the two roots
  #
  #   Args:
  #   -----
  #   a : float
  #     co-efficient
  #   b : float
  #     co-efficient
  #   c : float
  #     co-efficient
  #
  #   Returns:
  #   -------
  #   val : float
  #     positive root
  #

  # discriminant
  d <- b**2.0 - 4.0 * a * c

  if (d < 0.0) {
    stop("imaginary root found")
  }

  if (large) {

    if ( (is_close(a, 0.0)) & (b > 0.0) ) {
      root <- -c / b
    } else if ( (is_close(a, 0.0)) & (is_close(b, 0.0)) ) {
        root <- 0.0
        if (c != 0.0) {
          stop("Cant solve quadratic")
        }
    } else {
        root <- (-b + sqrt(d)) / (2.0 * a)
    }

  } else {

    if ( (is_close(a, 0.0)) & (b > 0.0) ) {
      root <- -c / b
    } else if ( (is_close(a, 0.0)) & (is_close(b, 0.0)) ) {
      root <- 0.0
      if (c != 0.0) {
        stop('Cant solve quadratic')
      }
    } else {
      root <- (-b - sqrt(d)) / (2.0 * a)
    }

  }

  return ( root )
}
