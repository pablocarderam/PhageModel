
# ###
# Estimates growth rates from absorbance data and carrying capacity from dilution plate data
# Coded by Pablo CR, pablocarderam@gmail.com
# Aug 2017
# ###

# Clear workspace
rm(list = ls(all = TRUE))

library(readr) # reads files
library(ggplot2) # C'mon, who uses R's native graphs?
source("src/dep/fancyScientific.R") # Needed for graphing
source("dat/phageParams.R") # needed for ODE
source("src/phageSolveODE.R") # needed for ODE

## Colorblind-friendly palette with black, 
# taken from Winston Chang's http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
# developed by Okabe M. and Ito K. http://jfly.iam.u-tokyo.ac.jp/color/
cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Shape palette
shapePalette = 1:20#c("a","b","c","d","e","f","g","h","i","j")

plotExp = function(dat,series_labs,series_cols,model) {
  pts = dat$time*24 # graph in hours
  err = dat$time*24
  for(i in 1:length(series_labs)) { # add each time series
    pts = cbind.data.frame(pts,1-dat[,2*i]) # we want resistance, not susceptibility
    err = cbind.data.frame(err,dat[,2*i+1])
  }
  
  names(pts) = c("time",series_labs) # Add labels for legend
  pts = melt(pts,id.vars="time") # reshape for ggplot
  names(err) = c("time",series_labs) # Add labels for legend
  err = melt(err,id.vars="time") # reshape for ggplot
  
  pts = cbind.data.frame(pts,err$value) # combine into a single df
  names(pts)[4] = "Err"
  
  series_shape = shapePalette[1:length(unique(pts$variable))]
  
  graph = ggplot(data=pts,aes(x=time,y=value,colour=variable,shape=variable))+ 
    geom_errorbar(aes(ymin=value-Err,ymax=value+Err),width=7,size=0.75,alpha=0.5)+
    geom_line(data=model,aes(x=time,y=value,colour=variable),linetype="dashed", show.legend = FALSE, size=1,alpha=0.5)+
    geom_line(aes(x=time,y=value,colour=variable),size=0.75)+
    geom_point(aes(x=time,y=value,colour=variable,shape=variable),size=2.5)+ # graph on top of everything else
    labs(x="Time (h)",y="Fraction Resistant to San23")+
    theme_bw(9)+
    scale_colour_manual(name="",labels=c( series_labs,"Single phage, no antibiotic (model)   ","Single phage, with antibiotic (model)   ","Two phages, no antibiotic (model)   " ),values=series_cols)+ # "" Removes legend title TODO: THIS IS A HACK for the labels
    scale_shape_manual(name="",labels=c( series_labs,"Single phage, no antibiotic (model)   ","Single phage, with antibiotic (model)   ","Two phages, no antibiotic (model)   " ),values=c(series_shape,rep(NA,3)))+ # "" Removes legend title
    guides(fill = guide_legend(override.aes = list(linetype = 0, shape='')), colour = guide_legend(override.aes = list(linetype=c(rep(1,length(series_labs)),rep(2,3)))))
    theme(legend.position = "right",legend.text=element_text(size=7))+
    theme(legend.margin=margin(t=-0.2, r=0, b=0, l=0, unit="cm")) # move legend up
  
  return(graph)
}

dat6 = read_csv("dat/PhageResData20162.csv")
dat7 = read_csv("dat/PhageResData20171.csv")

param = getParameters()
transfer_times = c(24,48,72,120,144,168,192) # Transfer times in hours. Add one extra transfer time at end as an endpoint
solsTransfers = runTransferExperiment(param,transfer_times)[c(1,2,3)] # do not plot two phages and antibiotic, also change order for legend

frac1 = solsTransfers[[1]]$time
for(sol in solsTransfers) {
  frac1 = cbind.data.frame(frac1,sol$fracRes1) # Get only data used in this graph
}

frac1 = frac1[seq(1, nrow(frac1), 10),] # get every 10th row to render dashed lines correctly

names(frac1) = c("time","Single phage, no antibiotic (model)   ","Single phage, with antibiotic (model)   ","Two phages, no antibiotic (model)   ") # Add labels for legend
mod = melt(frac1,id.vars="time") # reshape for ggplot

pal6 = c(cbPalette[1:4],cbPalette[2:5]) # color palette to be used, including colors for model time series
pal7 = c(cbPalette[1:6],cbPalette[2:5]) # color palette to be used, including colors for model time series
plt6 = plotExp(dat6,c("No Phage   ","San23   ","San23 + Chloramphenicol","San23 + San15   "),pal6,mod)
plt7 = plotExp(dat7,c("No Phage   ","San23   ","San23 + Chloramphenicol","San23 + San15   ","San23 + San24   ","San23 + San15 + San24   "),pal7,mod)

g6 = ggplotGrob(plt6)
g7 = ggplotGrob(plt7)

lay = rbind(c(1),c(2)) # layout

gExp = arrangeGrob(g6,g7, layout_matrix = lay) # combine both plots
ggsave(filename="Experiment.png",path="plt/",plot=gExp,width = 8, height = 6, units = c("in"), dpi = 300)
