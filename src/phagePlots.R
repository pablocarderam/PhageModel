
### GRAPHS

# Imports
library(ggplot2) # needed for pretty plots
library(gtable) # needed to align plots, use rbind for grobs
library(grid) # needed to align plots
library(gridExtra) # needed to align plots
library(reshape2) # allows easy ggplot data input format
source("src/dep/multiplot.R") # needed for plots
source("src/dep/fancyScientific.R") # needed for plots

## Colorblind-friendly palette with black, 
# taken from Winston Chang's http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
# developed by Okabe M. and Ito K. http://jfly.iam.u-tokyo.ac.jp/color/
cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Single Phage
singlePhageGraphs = function(sol,title,legend,y_axis) {
  # Phage Compartment
  phaPltDF = sol[,c("time","V1")] # Get only data used in this graph
  names(phaPltDF) = c("time","Phage") # Add labels for legend
  phaLegend = c(cbPalette[4]) # select colors
  phaPltDF = melt(phaPltDF,id.vars="time") # reshape for ggplot
  
  phaPlt = ggplot(phaPltDF)+ # Plot
    geom_line(aes(x=time,y=value,colour=variable),size=0.75)+
    scale_colour_manual("",values=phaLegend)+ # "" Removes legend title
    scale_y_continuous(labels=fancy_scientific)+
    labs(x="Time (h)",y="Phage (PFU/mL)")+
    theme_bw(7) + theme(legend.key.height=unit(0.5,"line"))+ # 24 is label text size
    ggtitle(title)+
    theme(axis.title.x=element_blank())
  
  # Bacterial Compartments
  compPltDF = sol[,c("time","S","R1","R12","I1")] # Get only data used in this graph
  compPltDF$R = compPltDF$R1 + compPltDF$R12 # add the two compartments containing resistants to phage 1
  compPltDF = compPltDF[,c("time","S","R","I1")] # Get only data used in this graph
  names(compPltDF) = c("time","Susceptible","Resistant","Infected") # Add labels for legend
  compLegend = c(cbPalette[3],cbPalette[7],cbPalette[2]) # select colors
  
  compPltDF = melt(compPltDF,id.vars="time") # reshape for ggplot
  
  bacPlt = ggplot(compPltDF)+ # Plot
    geom_line(aes(x=time,y=value,colour=variable),size=0.75)+
    scale_colour_manual("",values=compLegend)+ # "" Removes legend title
    scale_y_continuous(labels=fancy_scientific)+
    labs(x="Time (h)",y="Bacteria (CFU/mL)")+
    theme_bw(7) + theme(legend.key.height=unit(0.5,"line")) # 24 is label text size
  
  if(!legend) { # if no legend is to be shown, get rid of it
    phaPlt = phaPlt + theme(legend.position = "none")
    bacPlt = bacPlt + theme(legend.position = "none")
  }
  
  if(!y_axis) { # if y axis is not to be shown, get rid of y axis labels
    phaPlt = phaPlt + theme(axis.title.y=element_blank())
    bacPlt = bacPlt + theme(axis.title.y=element_blank())
  }
  
  return(list(phaPlt,bacPlt))
}

## Dual Phage
dualPhageGraphs = function(sol,title,legend,y_axis) {
  # Phage Compartments
  phaPltDF = sol[,c("time","V1","V2")] # Get only data used in this graph
  names(phaPltDF) = c("time","Phage 1","Phage 2") # Add labels for legend
  phaLegend = c(cbPalette[4],cbPalette[6]) # select colors
  
  phaPltDF = melt(phaPltDF,id.vars="time") # reshape for ggplot
  
  phaPlt = ggplot(phaPltDF)+ # Plot
    geom_line(aes(x=time,y=value,colour=variable),size=0.75)+
    scale_colour_manual("",values=phaLegend)+ # "" Removes legend title
    scale_y_continuous(labels=fancy_scientific)+
    labs(x="Time (h)",y="Phage (PFU/mL)")+
    theme_bw(7) + theme(legend.key.height=unit(0.5,"line"))+ # 18 is label text size
    ggtitle(title)+
    theme(axis.title.x=element_blank())
  
  # Large Bacterial Compartments
  bigBactPltDF = sol[,c("time","S","R1","R12","I1","I2")] # Get only data used in this graph
  names(bigBactPltDF) = c("time","Susceptible","Resistant to Phage 1 only","Resistant to both phages","Infected with Phage 1","Infected with Phage 2") # Add labels for legend
  bigBactLegend = c(cbPalette[3],cbPalette[7],cbPalette[1],cbPalette[2],cbPalette[5]) # select colors
  
  bigBactPltDF = melt(bigBactPltDF,id.vars="time") # reshape for ggplot
  
  bigBactPlt = ggplot(bigBactPltDF)+ # Plot
    geom_line(aes(x=time,y=value,colour=variable),size=0.75)+
    scale_colour_manual("",values=bigBactLegend)+ # "" Removes legend title
    scale_y_continuous(labels=fancy_scientific)+
    labs(x="Time (h)",y="Bacteria (CFU/mL)")+
    theme_bw(7) + theme(legend.key.height=unit(0.5,"line"))+ # 18 is label text size
    theme(axis.title.x=element_blank())
  
  # Small Bacterial Compartments
  smlBactPltDF = sol[,c("time","R2")] # Get only data used in this graph
  names(smlBactPltDF) = c("time","Resistant to Phage 2 only") # Add labels for legend
  smlBactLegend = c(cbPalette[8]) # select colors
  
  smlBactPltDF = melt(smlBactPltDF,id.vars="time") # reshape for ggplot
  
  smlBactPlt = ggplot(smlBactPltDF)+ # Plot
    geom_line(aes(x=time,y=value,colour=variable),size=0.75)+
    scale_colour_manual("",values=smlBactLegend)+ # "" Removes legend title
    scale_y_continuous(labels=fancy_scientific)+
    labs(x="Time (h)",y="Bacteria (CFU/mL)")+
    theme_bw(7) + theme(legend.key.height=unit(0.5,"line")) # 18 is label text size
  
  if(!legend) { # if no legend is to be shown, get rid of it
    phaPlt = phaPlt + theme(legend.position = "none")
    bigBactPlt = bigBactPlt + theme(legend.position = "none")
    smlBactPlt = smlBactPlt + theme(legend.position = "none")
  }
  
  if(!y_axis) { # if y axis is not to be shown, get rid of y axis labels
    phaPlt = phaPlt + theme(axis.title.y=element_blank())
    bigBactPlt = bigBactPlt + theme(axis.title.y=element_blank())
    smlBactPlt = smlBactPlt + theme(axis.title.y=element_blank())
  }
  
  return(list(phaPlt,bigBactPlt,smlBactPlt))
}

getFullLegend = function(df,labels,color_key) {
  names(df) = labels # Add labels for legend
  dfLegend = color_key # select colors
  
  df = melt(df,id.vars="time") # reshape for ggplot
  
  dummyPlt = ggplot(df)+ # Plot
    geom_line(aes(x=time,y=value,colour=variable),size=0.75)+
    scale_colour_manual("",values=dfLegend)+ # "" Removes legend title
    theme_bw(7)+
    theme(legend.key.height=unit(0.5,"line"))+
    theme(legend.text=element_text(size=8))+
    theme(legend.position = "bottom",legend.text=element_text(size=8))+
    guides(colour=guide_legend(nrow=2,byrow=F))
  
  
  tmp = ggplot_gtable(ggplot_build(dummyPlt))
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  
  return(legend)
}

# Combine compartment graphs
plotCompGraphs = function(sols) {
  g1P0A = singlePhageGraphs(sols[[1]],"Single Phage, No Antibiotic",legend=F,y_axis=T) # Single Phage, No Antibiotic TODO:fix singlephagegraphs()
  g1P1A = singlePhageGraphs(sols[[2]],"Single Phage, Antibiotic",legend=F,y_axis=F) # Single Phage, Antibiotic
  g2P0A = dualPhageGraphs(sols[[3]],"Two Phages, No Antibiotic",legend=F,y_axis=T) # Dual Phage, No Antibiotic
  g2P1A = dualPhageGraphs(sols[[4]],"Two Phages, Antibiotic",legend=F,y_axis=F) # Dual Phage, Antibiotic
  
  g0A = lapply(c(g1P0A,g2P0A), function(g) {
    g = ggplotGrob(g)
    g = gtable_add_cols(g, unit(0,"mm")) # add a column for missing legend
    g = gtable_add_cols(g, unit(0,"mm")) # add a column for missing legend
  })
  
  g1A = lapply(c(g1P1A,g2P1A), function(g) {
    g = ggplotGrob(g)
  })
  
  g0A = rbind(g0A[[1]],g0A[[2]],g0A[[3]],g0A[[4]],g0A[[5]],size="max") # stack the plots
  g1A = rbind(g1A[[1]],g1A[[2]],g1A[[3]],g1A[[4]],g1A[[5]],size="max") # stack the plots
  g0A$widths = unit.pmax(g0A$widths,g1A$widths) # use the largest widths
  g1A$widths = unit.pmax(g0A$widths,g1A$widths) # use the largest widths
  
  legendDF = sols[[3]][,c("time","S","R1","R2","R12","I1","I2","V1","V2")] # Get only data used in this graph
  
  leg = getFullLegend(legendDF,c("time","Susceptible","Resistant to Phage 1 only   ","Resistant to Phage 2 only   ","Resistant to both phages  ","Infected with Phage 1      ","Infected with Phage 2      ","Phage 1","Phage 2"),c(cbPalette[3],cbPalette[7],cbPalette[8],cbPalette[1],cbPalette[2],cbPalette[5],cbPalette[4],cbPalette[6]) )
  
  lay = rbind(c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(3,3)) # layout with legend at right, small
  
  gCompartments = arrangeGrob(g0A,g1A,leg, layout_matrix = lay) # combine all compartment plots
  
  return(gCompartments)
}

## Different Treatment Comparison
plotTotalGraphs = function(sols) {
  # Total Population
  totPopPltDF = sols[[1]]$time
  for(sol in sols) {
    totPopPltDF = cbind.data.frame(totPopPltDF,sol$N) # Get only data used in this graph
  }
  
  names(totPopPltDF) = c("time","Single phage, no antibiotic","Single phage, antibiotic","Two phages, no antibiotic","Two phages, antibiotic") # Add labels for legend
  totPop1PhaLegend = c(cbPalette[2],cbPalette[3]) # select colors
  totPop2PhaLegend = c(cbPalette[4],cbPalette[5]) # select colors
  
  totPopPltDF = melt(totPopPltDF,id.vars="time") # reshape for ggplot
  
  totPop1PhaPltDF = subset(totPopPltDF,variable=="Single phage, no antibiotic" | variable=="Single phage, antibiotic")
  totPop2PhaPltDF = subset(totPopPltDF,variable=="Two phages, no antibiotic" | variable=="Two phages, antibiotic")
  
  pop1PhaPlt = ggplot(data=totPop1PhaPltDF)+
    geom_line(aes(x=time,y=value,colour=variable),size=0.75)+
    scale_colour_manual("",values=totPop1PhaLegend)+ # "" Removes legend title
    scale_y_continuous(labels=fancy_scientific)+
    labs(x="",y="Total Bacteria \n(CFU/mL)")+
    theme_bw(7) + theme(legend.key.height=unit(0.5,"line"))+ # 18 is label text size
    theme(axis.title.x=element_blank())+
    theme(legend.position = "none")
  
  pop2PhaPlt = ggplot(data=totPop2PhaPltDF)+
    geom_line(aes(x=time,y=value,colour=variable),size=0.75)+
    scale_colour_manual("",values=totPop2PhaLegend)+ # "" Removes legend title
    scale_y_continuous(labels=fancy_scientific)+
    labs(x="",y="Total Bacterial \n(CFU/mL)")+
    theme_bw(7) + theme(legend.key.height=unit(0.5,"line"))+ # 18 is label text size
    theme(axis.title.x=element_blank())+
    theme(legend.position = "none")
  
  # Resistance to Phage 1 Fractions
  fra1PopPltDF = sols[[1]]$time
  for(sol in sols) {
    fra1PopPltDF = cbind.data.frame(fra1PopPltDF,sol$fracRes1) # Get only data used in this graph
  }
  
  names(fra1PopPltDF) = c("time","Single phage, no antibiotic","Single phage, antibiotic","Two phages, no antibiotic","Two phages, antibiotic") # Add labels for legend
  fra1PopLegend = c(cbPalette[2],cbPalette[3],cbPalette[4],cbPalette[5]) # select colors
  fra1Pop1PhaLegend = c(cbPalette[2],cbPalette[3]) # select colors
  fra1Pop2PhaLegend = c(cbPalette[4],cbPalette[5]) # select colors
  
  fra1PopPltDF = melt(fra1PopPltDF,id.vars="time") # reshape for ggplot
  
  fra1Pop1PhaPltDF = subset(fra1PopPltDF,variable=="Single phage, no antibiotic" | variable=="Single phage, antibiotic")
  fra1Pop2PhaPltDF = subset(fra1PopPltDF,variable=="Two phages, no antibiotic" | variable=="Two phages, antibiotic")
  
  fra11PhaPlt = ggplot(data=fra1Pop1PhaPltDF)+
    geom_line(aes(x=time,y=value,colour=variable),size=0.75)+
    scale_colour_manual("",values=fra1Pop1PhaLegend)+ # "" Removes legend title
    labs(x="",y="Fraction Resistant \nto Phage 1")+
    theme_bw(7) + theme(legend.key.height=unit(0.5,"line"))+ # 18 is label text size
    theme(axis.title.x=element_blank())+
    theme(legend.position = "none")
  
  fra12PhaPlt = ggplot(data=fra1Pop2PhaPltDF)+
    geom_line(aes(x=time,y=value,colour=variable),size=0.75)+
    scale_colour_manual("",values=fra1Pop2PhaLegend)+ # "" Removes legend title
    labs(x="",y="Fraction Resistant \nto Phage 1")+
    theme_bw(7) + theme(legend.key.height=unit(0.5,"line"))+ # 18 is label text size
    theme(axis.title.x=element_blank())+
    theme(legend.position = "none")
  
  # Resistance to Treatment Fractions
  fraTPopPltDF = sols[[1]]$time
  for(sol in sols) {
    fraTPopPltDF = cbind.data.frame(fraTPopPltDF,sol$fracResT) # Get only data used in this graph
  }
  
  names(fraTPopPltDF) = c("time","Single phage, no antibiotic   ","Single phage, antibiotic   ","Two phages, no antibiotic   ","Two phages, antibiotic   ") # Add labels for legend
  fraTPopLegend = c(cbPalette[2],cbPalette[3],cbPalette[4],cbPalette[5]) # select colors
  
  fraTPopPltDF = melt(fraTPopPltDF,id.vars="time") # reshape for ggplot
  
  fraTPlt = ggplot(data=fraTPopPltDF)+
    geom_line(aes(x=time,y=value,colour=variable),size=0.75)+
    scale_colour_manual("",values=fraTPopLegend)+ # "" Removes legend title
    labs(x="Time (h)",y="Fraction Resistant \nto Phages 1 and 2")+
    theme_bw(7) + theme(legend.key.height=unit(0.5,"line"))+ # 18 is label text size
    theme(legend.position = "bottom",legend.text=element_text(size=8))+
    guides(colour=guide_legend(nrow=2,byrow=TRUE))
  
  # Align and plot graphs as per https://github.com/tidyverse/ggplot2/wiki/Align-two-plots-on-a-page combining two plots
  gPop1Pha = ggplotGrob(pop1PhaPlt)
  gPop2Pha = ggplotGrob(pop2PhaPlt)
  gFra11Pha = ggplotGrob(fra11PhaPlt)
  gFra12Pha = ggplotGrob(fra12PhaPlt)
  gFraT = ggplotGrob(fraTPlt)
  
  gComparison = rbind(gPop1Pha,gPop2Pha,gFra11Pha,gFra12Pha,gFraT,size="first") # stack the two plots
  gComparison$widths = unit.pmax(gPop1Pha$widths,gPop2Pha$widths,gFra11Pha$widths,gFra12Pha$widths,gFraT$widths) # use the largest widths
  
  return(gComparison)
}

# Plot all:
gCompartments = plotCompGraphs(sols)
gTotals = plotTotalGraphs(sols)

ggsave(filename="compartments.png",path="plt/",plot=gCompartments,width = 7, height = 6, units = c("in"), dpi = 300)
ggsave(filename="totals.png",path="plt/",plot=gTotals,width = 4, height = 6, units = c("in"), dpi = 300)
