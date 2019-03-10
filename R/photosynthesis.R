
source("utils.R")

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
