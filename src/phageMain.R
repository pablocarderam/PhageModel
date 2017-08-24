
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

transfer_times = c(24,48,72,96,120,144,168,192) # Transfer times in hours. Add one extra transfer time at end as an endpoint
solsTransfers = runTransferExperiment(param,transfer_times)

solsCulture = lapply(solsCulture, function (sol) {
  sol = sol[seq(1, nrow(sol), 10),] # get every 10th row to render dashed lines correctly
  return(sol)
})

solsTransfers = lapply(solsTransfers, function (sol) {
  sol = sol[seq(1, nrow(sol), 10),] # get every 10th row to render dashed lines correctly
  return(sol)
})

# Plot all:
gCompartments = plotCompGraphs(solsCulture) # graphs each compartment in model for every treatment
gTotalsBoth = plotAllTotals(solsCulture,solsTransfers,transfer_times) # compares treatments by graphing total populations and resistant fractions

ggsave(filename="Compartments.png",path="plt/",plot=gCompartments,width = 7, height = 6, units = c("in"), dpi = 300)
ggsave(filename="Totals.png",path="plt/",plot=gTotalsBoth,width = 8, height = 6, units = c("in"), dpi = 300)
