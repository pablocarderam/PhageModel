
# MAIN

# Clear workspace
rm(list = ls(all = TRUE))

# Imports
source("dat/phageParams.R") # needed for ODE
source("src/phageSolveODE.R") # needed for ODE
source("src/phagePlots.R") # needed for plots

# Run everything
param = getParameters()
setup_list = setup(param) # returns list with [[1]] list of starting vectors for each of the four treatments and [[2]] the time vector
solsCulture = runCulture(param,setup_list[[1]],setup_list[[2]]) # solves model numerically

transfer_times_uni = c(24,48,72,96,120,144,168,192) # Transfer times in hours. Add one extra transfer time at end as an endpoint
solsTransfers_uni = runTransferExperiment(param,transfer_times_uni)

transfer_times_exp = c(24,48,72,120,144,168,192) # Transfer times in hours. Add one extra transfer time at end as an endpoint
solsTransfers_exp = runTransferExperiment(param,transfer_times_exp)

# Plot all:
gCompartments = plotCompGraphs(solsCulture) # graphs each compartment in model for every treatment
gTotals = plotTotalGraphs(solsCulture,c(0)) # compares treatments by graphing total populations and resistant fractions

gTotalsTransf_uni = plotTotalGraphs(solsTransfers_uni,c(0,transfer_times_uni)) # compares treatments by graphing total populations and resistant fractions
gTotalsTransf_exp = plotTotalGraphs(solsTransfers_exp,c(0,transfer_times_exp)) # compares treatments by graphing total populations and resistant fractions

ggsave(filename="compartments.png",path="plt/",plot=gCompartments,width = 7, height = 6, units = c("in"), dpi = 300)
ggsave(filename="totals.png",path="plt/",plot=gTotals,width = 4, height = 6, units = c("in"), dpi = 300)
ggsave(filename="totalsTransf_Uni.png",path="plt/",plot=gTotalsTransf_uni,width = 4, height = 6, units = c("in"), dpi = 300)
ggsave(filename="totalsTransf_Exp.png",path="plt/",plot=gTotalsTransf_exp,width = 4, height = 6, units = c("in"), dpi = 300)
