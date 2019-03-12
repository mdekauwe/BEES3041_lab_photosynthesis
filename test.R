
source("R/photosynthesis.R")
source("R/parameters.R")
source("R/constants.R")

# Build a params object
p <- list(Kc25=Kc25, Ko25=Ko25, Eo=Eo, Ec=Ec, Oi=Oi, gamstar25=gamstar25, Eag=Eag, Vcmax25=Vcmax25, 
          Eav=Eav, deltaSv=deltaSv, Hdv=Hdv)

# Met variables ...
Tleaf <- 25.0 + DEG_2_KELVIN

An <- calc_photosynthesis(p, Tleaf)

