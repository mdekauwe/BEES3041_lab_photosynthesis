
source("R/photosynthesis.R")
source("R/parameters.R")
source("R/constants.R")

# Build a params object
p <- list(Kc25=Kc25, Ko25=Ko25, Eo=Eo, Ec=Ec, Oi=Oi, gamstar25=gamstar25,
          Eag=Eag, Vcmax25=Vcmax25, Eav=Eav, deltaSv=deltaSv, Hdv=Hdv,
          Jmax25=Jmax25, Eaj=Eaj, deltaSj=deltaSj, Hdj=Hdj, theta_J=theta_J,
          alpha=alpha)

# Met variables ...
Tleaf <- 25.0 + DEG_2_KELVIN
PAR <- 1800.0
An <- calc_photosynthesis(p, Tleaf, PAR, peaked_Vcmax=TRUE, peaked_Jmax=TRUE)
