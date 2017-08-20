
# MAIN

# Clear workspace
rm(list = ls(all = TRUE))

# Imports
source("dat/phageParams.R") # needed for ODE
source("src/phageSolveODE.R") # needed for ODE
source("src/phagePlots.R") # needed for plots

# Run everything
param = getParameters()
sols = runOneDay(param)

# Plot all:
gCompartments = plotCompGraphs(sols)
gTotals = plotTotalGraphs(sols)

ggsave(filename="compartments.png",path="plt/",plot=gCompartments,width = 7, height = 6, units = c("in"), dpi = 300)
ggsave(filename="totals.png",path="plt/",plot=gTotals,width = 4, height = 6, units = c("in"), dpi = 300)
