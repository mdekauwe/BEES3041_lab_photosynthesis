
library(ggplot2)
source("R/photosynthesis.R")
source("R/parameters.R")
source("R/constants.R")


# Met variables ...
Tleaf <- 25.0 + DEG_2_KELVIN
PAR <- 1800.0
Cs <- 400.
vpd <- 1.5

Tleaf <- seq(1, 50.0, 0.5) + DEG_2_KELVIN
Vcmax <- peaked_arrh(p$Vcmax25, p$Eav, Tleaf, p$deltaSv, p$Hdv)
Jmax <- peaked_arrh(p$Jmax25, p$Eaj, Tleaf, p$deltaSj, p$Hdj)

df <- data.frame(Tleaf, Vcmax, Jmax)

ggplot(df, aes(Tleaf-DEG_2_KELVIN)) +
  geom_line(aes(y=Vcmax, colour="Vcmax")) +
  geom_line(aes(y=Jmax, colour="Jmax")) +
  ylab(expression("Paramater" ~ (mu * mol ~  m^{-2}  *  s^{-1}))) +
  xlab(expression('Temperature ('*~degree*C*')')) + 
  theme_classic()
  


out <- calc_photosynthesis(p, Tleaf, PAR, Cs, vpd, peaked_Vcmax=TRUE,
                           peaked_Jmax=TRUE)

print(out$gsc)
print(out$An)

