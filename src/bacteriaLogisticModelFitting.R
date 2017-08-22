
# ###
# Estimates growth rates from absorbance data and carrying capacity from dilution plate data
# Coded by Pablo CR, pablocarderam@gmail.com
# Aug 2017
# ###

# Clear workspace
rm(list = ls(all = TRUE))

library(readr) # reads files
library(car) # for initial parameter estimate
library(minpack.lm) # More robust to bad starting parameters than nls 
library(ggplot2) # C'mon, who uses R's native graphs?
source("src/dep/fancyScientific.R") # Needed for graphing

fitLogisticEquation = function(indepVar,depVar1,depVar2,depVar1_Err,depVar2_Err,nameDepVar1,nameDepVar2,est_K,title,title_x,title_y,col1,col2) {
  # Initial estimate of logistic parameters, 
  initial_est_asymp = est_K # assume initial estimate of upper asymptote
  initial_par_est=coef(lm(logit(depVar1/initial_est_asymp)~indepVar)) # initial estimate by regression as per https://bscheng.com/2014/05/07/modeling-logistic-growth-data-in-r/
  # nlsLM uses Levenberg-Marquardt algorithm for optimization (switches between Gauss-Newton and gradient descent, avoids singular gradient matrix problems)
  log_fit = nlsLM(depVar1~K/(1+exp(-(phi + r*indepVar))),start=list(K=initial_est_asymp,phi=initial_par_est[1],r=initial_par_est[2]),trace=TRUE)
  
  # set parameters
  K1=coef(log_fit)[1]
  phi1=coef(log_fit)[2]
  r1=coef(log_fit)[3]
  x=seq(min(indepVar),max(indepVar),length=10000) # construct a range of x values bounded by the data
  y1=K1/(1+exp(-(phi1+r1*x))) # predicted growth
  
  dat = data.frame(indepVar,depVar1,depVar2)
  K1_displayed = K1
  K1_label = "K = "
  if(K1>10^5 - 1) {
    K1_displayed = log10(K1)
    K1_label = "log(K) = "
  }
  
  # Initial estimate of logistic parameters, by regression as per https://bscheng.com/2014/05/07/modeling-logistic-growth-data-in-r/
  initial_par_est=coef(lm(logit(depVar2/initial_est_asymp)~indepVar)) 
  # nlsLM uses Levenberg-Marquardt algorithm for optimization (switches between Gauss-Newton and gradient descent, avoids singular gradient matrix problems)
  log_fit = nlsLM(depVar2~K/(1+exp(-(phi + r*indepVar))),start=list(K=initial_est_asymp,phi=initial_par_est[1],r=initial_par_est[2]),trace=TRUE)
  
  # Now plot the fit:
  # set parameters
  K2=coef(log_fit)[1]
  phi2=coef(log_fit)[2]
  r2=coef(log_fit)[3]
  y2=K2/(1+exp(-(phi2+r2*x))) # predicted growth 
  
  
  fit=data.frame(x,y1,y2) # create the fit data frame #And add a nice plot
  dat = data.frame(indepVar,depVar1,depVar2,depVar1_Err,depVar2_Err)
  names(dat) = c("indepVar","depVar1","depVar2","depVar1_E","depVar2_E")
  K2_displayed = K2
  K2_label = "K = "
  if(K2>10^5 - 1) {
    K2_displayed = log10(K2)
    K2_label = "log(K) = "
  }
  
  df = cbind.data.frame(indepVar,depVar1,depVar2)
  names(df) = c("time",nameDepVar1,nameDepVar2) # Add labels for legend
  leg = c(col1,col2) # select colors
  df = melt(df,id.vars="time") # reshape for ggplot
  
  dfErr = cbind.data.frame(indepVar,depVar1_Err,depVar2_Err)
  names(dfErr) = c("time",nameDepVar1,nameDepVar2) # Add labels for legend
  dfErr = melt(dfErr,id.vars="time") # reshape for ggplot
  
  df = cbind.data.frame(df,dfErr$value)
  names(df)[4] = "Err" # Add labels for legend
  
  graph = ggplot(data=df,aes(x=time,y=value,colour=variable))+
    geom_point(aes(x=time,y=value,colour=variable),size=2)+
    geom_errorbar(aes(ymin=value-Err,ymax=value+Err),width=0.6,size=0.75,alpha=0.75)+
    geom_line(data=fit,aes(x=x,y=y1), size=1,color=col1,alpha=0.5)+
    geom_line(data=fit,aes(x=x,y=y2), size=1,color=col2,alpha=0.5)+
    scale_colour_manual("",values=leg)+ # "" Removes legend title
    theme_bw(7)+labs(x=title_x,y=title_y)+
    ggtitle(title)
  
  return(list(graph,K1,r1,phi1,K2,r2,phi2))
}

GrowthCurveDataAbs = read_csv("dat/GrowthCurveDataAbsorbance.csv")[1:10,] # Last datapoint discarded for goodness of fit
GrowthCurveDataDil = read_csv("dat/GrowthCurveDataDilution.csv")

# This used to get estimate of r using absorbance curves:
out1 = fitLogisticEquation(GrowthCurveDataAbs$Time,GrowthCurveDataAbs$NoPhage,GrowthCurveDataAbs$Antibiotic,GrowthCurveDataAbs$NoPhage_StdErr,GrowthCurveDataAbs$Antibiotic_StdErr,"No Antibiotic   ","Antibiotic",0.6,"","Time (h)","Absorbance (AU)",'#407ffc','#f1ac00')
abs = out1[[1]] # get graph
penalty = 1 - round(out1[[6]]/out1[[3]],5) # get penalty of antibiotic on reproduction
abs = abs + annotate(parse=TRUE,"text", x=0.2*max(GrowthCurveDataAbs$Time), y=0.9*max(GrowthCurveDataAbs$NoPhage), size=3, color="#505050", label="'Specific growth rate ('*italic(r)*')'" )
abs = abs + annotate("text", x=0.2*max(GrowthCurveDataAbs$Time), y=0.8*max(GrowthCurveDataAbs$NoPhage), size=3, color="#505050", label="fractional decrease in presence" )
abs = abs + annotate("text", x=0.2*max(GrowthCurveDataAbs$Time), y=0.7*max(GrowthCurveDataAbs$NoPhage), size=3, color="#505050", label=paste("of antibiotic: ", penalty) )
abs = abs + theme(legend.position = "bottom",legend.text=element_text(size=7))

# This used to get estimate of CFU at K:
out2 = fitLogisticEquation(GrowthCurveDataDil$Time,GrowthCurveDataDil$NoPhage/10,GrowthCurveDataDil$NoPhage/10,GrowthCurveDataDil$NoPhage_StdErr/10,GrowthCurveDataDil$NoPhage_StdErr/10,"No Antibiotic","No Antibiotic",10^8,"","Time (h)","CFU/mL",'#407ffc','#407ffc')
dil = out2[[1]] # get graph
dil = dil + annotate("text",x=0.2*max(GrowthCurveDataDil$Time),y=0.9*max(GrowthCurveDataDil$NoPhage/10),parse=TRUE,size=3,color="#505050",label=paste("italic(K)*' = '*",fancy_scientific(out2[[2]]),"*' CFU/mL'" ) )
dil = dil + annotate("text",x=0.2*max(GrowthCurveDataDil$Time),y=0.8*max(GrowthCurveDataDil$NoPhage/10),parse=TRUE,size=3,color="#505050",label=paste("italic(r)*' = '*",round(out2[[3]],5),"*' CFU (mL h)'^{-1}" ) )
dil = dil + annotate("text",x=0.2*max(GrowthCurveDataDil$Time),y=0.7*max(GrowthCurveDataDil$NoPhage/10),parse=TRUE,size=3,color="#505050",label=paste("italic(φ)*' = '*",round(out2[[4]],5)) )
dil = dil + scale_y_continuous(labels=fancy_scientific)
dil = dil + theme(legend.position = "none")

gDil = ggplotGrob(dil)
gAbs = ggplotGrob(abs)

lay = rbind(c(1),c(2)) # layout

gFit = arrangeGrob(gDil,gAbs, layout_matrix = lay) # combine both plots
ggsave(filename="LogisticFit.png",path="plt/",plot=gFit,width = 5, height = 5, units = c("in"), dpi = 300)
