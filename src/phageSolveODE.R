
# ###
# Solves ODE system for bacterial resistance to phages
# Coded by Pablo CR, pablocarderam@gmail.com
# Aug 2017
# ###

# Imports
library(deSolve) # contains numerical ODE solvers

setup = function(pParam) {
  with(as.list(pParam),{ # replace parameters into starting conditions
    rho = k                 # Bacteria in inoculation broth ( CFU/mL )
    N0 = rho*d              # Starting total bacteria ( CFU/mL )
    
    S0 = N0*(1-2*m-m^2)     # Starting susceptible bacteria ( CFU/mL )
    R10 = N0*m              # Starting bacteria resistant to phage 1 ( CFU/mL )
    R20 = N0*m              # Starting bacteria resistant to phage 1 ( CFU/mL )
    R120 = N0*m^2           # Starting bacteria resistant to phages 1 and 2 ( CFU/mL )
    I10 = 0                 # Starting bacteria infected with phage 1 ( CFU/mL )
    I20 = 0                 # Starting bacteria infected with phage 2 ( CFU/mL )
    V10 = N0*moi1           # Starting phage 1 virions ( PFU/mL )
    V20 = N0*moi2           # Starting phage 1 virions ( PFU/mL )
    
    # Initial conditions vector
    X0_2Phage = c(S=S0,R1=R10,R2=R20,R12=R120,I1=I10,I2=I20,V1=V10,V2=V20) # one vector for the dual phage case
    X0_1Phage = X0_2Phage # one vector for the single phage case
    X0_1Phage["V2"] = 0 # in which there is no Phage 2
    
    ## Time Vector
    t_vec = seq(0,t_max,t_step) # time vector
    
    initialConditions = list(X0_1Phage,X0_1Phage,X0_2Phage,X0_2Phage) # list with initial conditions in order required
    returnList = list(initialConditions,t_vec) # list to be returned by function
    
    return(returnList)
  })
}

### NUMERICAL SOLUTION
## Derivative Function Definition
dX = function(t,X,pParam){
  # Get state variables
  S = X[1]
  R1 = X[2]
  R2 = X[3]
  R12 = X[4]
  I1 = X[5]
  I2 = X[6]
  V1 = X[7]
  V2 = X[8]
  
  with(as.list(pParam),{ # replace parameters into differential equations
    
    # Differential equations
    dS   = r*(1-alp*a)*S  *(1-(S+R1+R2+R12+I1+I2)/k) + (m-m^2)*R1 + (m-m^2)*R2 + m^2*R12 - (2*m-m^2)*S  - bet1*S*V1 - bet2*S*V2
    dR1  = r*(1-alp*a)*R1 *(1-(S+R1+R2+R12+I1+I2)/k) + (m-m^2)*S  + m^2*R2 + (m-m^2)*R12 - (2*m-m^2)*R1 - bet2*R1*V2
    dR2  = r*(1-alp*a)*R2 *(1-(S+R1+R2+R12+I1+I2)/k) + (m-m^2)*S  + m^2*R1 + (m-m^2)*R12 - (2*m-m^2)*R2 - bet1*R2*V1
    dR12 = r*(1-alp*a)*R12*(1-(S+R1+R2+R12+I1+I2)/k) + (m-m^2)*R1 + (m-m^2)*R2 + m^2*S   - (2*m-m^2)*R12
    dI1  = bet1*(S+R2)*V1 - 1/lam1*I1
    dI2  = bet2*(S+R1)*V2 - 1/lam2*I2
    dV1  = b1*1/lam1*I1 - bet1*(S+R2+I1)*V1 - mu1*V1
    dV2  = b2*1/lam2*I2 - bet2*(S+R1+I2)*V2 - mu2*V2
    
    # Make list with state variable differentials
    dy = c(dS,dR1,dR2,dR12,dI1,dI2,dV1,dV2)
    list(dy)
  })
}

runCulture = function(pParam,initial_conditions,t_vec) {
  ### FINISH SETUP
  ## Initial Conditions
  with(as.list(pParam),{ # replace parameters into starting conditions
    paramAB = pParam # create separate vector for antibiotic case
    paramAB["a"] = 1 # add antibiotic
    
    ## Calculate System Solution Numerically (fourth order Runge Kutta implementation)
    sol_1Phage = as.data.frame(lsoda(initial_conditions[[1]],t_vec,dX,pParam)) # one solution with one phage
    sol_1Phage_AB = as.data.frame(lsoda(initial_conditions[[2]],t_vec,dX,paramAB)) # one solution with one phage + antibiotic
    sol_2Phage = as.data.frame(lsoda(initial_conditions[[3]],t_vec,dX,pParam)) # one solution with two phages
    sol_2Phage_AB = as.data.frame(lsoda(initial_conditions[[4]],t_vec,dX,paramAB)) # one solution with two phages + antibiotc
    
    sols = list(sol_1Phage,sol_1Phage_AB,sol_2Phage,sol_2Phage_AB) # store all solutions
    
    ## Calculate Additional Time Series
    sols = lapply(sols, function(sol) {
      sol$N = sol$S+sol$R1+sol$R2+sol$R12+sol$I1+sol$I2 # Total bacteria population
      
      sol$fracRes1 = (sol$R12 + sol$R1)/(sol$N) # Fraction of bacteria resistant to Phage 1
      sol$fracNonres1 = 1-sol$fracRes1 # Fraction of bacteria not resistant to Phage 1
      
      if(sum(sol$V2) == 0) { # if treatment is only one phage
        sol$fracResT = sol$fracRes1 # Fraction of bacteria resistant to treatment
      }
      else { # if treatment is two phages
        sol$fracResT = (sol$R12)/(sol$N) # Fraction of bacteria resistant to both phages
      }
      sol$fracNonresT = 1-sol$fracResT # Fraction of bacteria not resistant to treatment
      
      return(sol)
    })
    
    return(sols)
  })
}

runTransferExperiment = function(pParam, transfer_times) {
  param = pParam # meh
  setup_list = setup(param) # returns list with [[1]] list of starting vectors for each of the four treatments and [[2]] the time vector
  init_cond = setup_list[[1]] # list with initial conditions for each of the four treatments
  prev_max = 0 # stores previous maximum time value
  
  # initialize solution dfs
  sols = list()
  
  for(t in transfer_times) { # propagate the culture once every transfer
    print(paste("Started transfer at time",t))
    t_max = t # set max time as the next time there will be a transfer
    t_vec = seq(prev_max,t_max,param["t_step"]) # time vector
    newSols = runCulture(param,init_cond,t_vec) # solve system for this culture transfer
    
    if(length(sols) == 0) { # if solutions are empty
      sols = newSols # these are the first solutions
    }
    else { # if not,
      for(i in 1:length(sols)) { # for every solution df
        sols[[i]] = rbind.data.frame(sols[[i]],newSols[[i]]) # add this transfer's results to it
      }
    }
    
    init_cond = lapply(newSols, function(sol) { 
      # set these cultures' final conditions, multiplied by the dilution factor of the transfer,
        # to be the new starting conditions
      i_cond = c(S=sol$S[[nrow(sol)]]*param[["v"]],
                 R1=sol$R1[[nrow(sol)]]*param[["v"]],
                 R2=sol$R2[[nrow(sol)]]*param[["v"]],
                 R12=sol$R12[[nrow(sol)]]*param[["v"]],
                 I1=sol$I1[[nrow(sol)]]*param[["v"]],
                 I2=sol$I2[[nrow(sol)]]*param[["v"]],
                 V1=sol$V1[[nrow(sol)]]*param[["v"]],
                 V2=sol$V2[[nrow(sol)]]*param[["v"]]
                 ) 
      
      return(i_cond)
    })
    
    prev_max = t # set previous max time to this iteration's max time
  }
  
  return(sols)
}
