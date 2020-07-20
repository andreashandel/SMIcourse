###################################################
##a script that uses results produced by the accompanying R fitting scripts and produces figures 
##written by Andreas Handel (ahandel@uga.edu), last change 11/1/17
###################################################
rm(list=ls()) #this clears the workspace to make sure no leftover variables are floating around. Not strictly needed
library(dplyr)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#load data from hayden
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=read.csv('hayden96data.csv')

nobs = nrow(data) #number of observations
nmax=3; #number of doses to consider - all being fit at the same time

U0=4e8; #number of initial target cells  
Vmin=0.5; #level of detection for virus in log10 units according to paper

#various plotting settings for all plots
cexaxis=1;  cexlab=1; cexval=1;  lw1=2;  lw2=2; symsize=1.2;

upos=1; ipos=2; vpos=3;  #position/index for the different variables

plotfigvec=c(1) #specify which plots to produce

load("haydenfit-timeseries-e1.Rdata")
e1model = timeseries.list
load("haydenfit-timeseries-e2.Rdata")
e2model = timeseries.list

txstatusvec=c('earlytx','latetx','notx')
txtimevec=c(29/24,50/24,500/24) #time of treatment start - last value means never


###################################################
#plotting data and best fit of model 
###################################################
if (1 %in% plotfigvec) #plotting best fit figure
{ 
  
  graphics.off(); #close all graphics windows 
  ww=5; wh=1/1*ww; #size is 17.8cm, needs to be in inches
  windows(width=ww, height=wh) #for windows: opens window of the specified dimensions
  
  
  #par(oma=c(0.2,0.2,0.2,0.2)) #outer margins 
  #par(mar=c(0, 0, 0, 0)) #bottom, left, top, right margins
  #par(mai=c(0.5, 0.5, 0.1, 0.1)) #bottom, left, top, right margins
  
  tmaxplot=8;
  
  mycol=c('green','red','blue','black','orange')
  mylty=c(1,2,3,4,5)
  mypch=c(21,22,23,24,25)
  
  ymin=1e-1; ymax=1e4;
  
  #do it for 1st model
  for (n in 1:nmax)
  {
    oderes=e1model[[n]]
    timemodel=oderes[,'time']
    resmodel=oderes[,'V']
    if (n==1) #set up plot for 1st entry
    {
      plot(timemodel,resmodel,lwd=1,type='l',lty=mylty[1],col=mycol[n],xlim=c(0,tmaxplot),ylim=c(ymin,ymax),xlab='',ylab='',log="y",cex=cexval,cex.lab=cexlab,cex.axis=cexaxis)
    }
    lines(timemodel,resmodel,lwd=lw1,lty=mylty[1],col=mycol[n]) 
  
    txstatus = txstatusvec[n]
    txtime = txtimevec[n]
    timedata <- data %>% filter(Condition==txstatus) %>% select(DaysPI) %>% unlist()  %>% as.numeric()
    virusdata <- data %>% filter(Condition==txstatus) %>% select(LogVirusLoad) %>% unlist()  %>% as.numeric()
    
    points(timedata,10^virusdata,pch=mypch[n],col=mycol[n],bg=mycol[n],cex=symsize)
    
  }         
  #do again for 2nd model
  for (n in 1:nmax)
  {
    oderes=e2model[[n]]
    timemodel=oderes[,'time']
    resmodel=oderes[,'V']
    lines(timemodel,resmodel,lwd=lw1,lty=mylty[2],col=mycol[n]) 
  }         
  
  
  #lines(timemodel,rep(10^Vmin,length(timemodel)),lty=2,lwd=1,col='black') #limit of detection
  mtext("Days post infection",side=1,line=2.5) #x-axis - change line to move in/out
  mtext("Virus load",side=2,line=2) #y-axis
  
  legend('topleft',legend=c('early tx','late tx','no tx'),col=mycol,pch=mypch,lty=mylty, lwd=2, pt.bg = mycol)
  
  #save in various formats
  #dev.print(device=pdf,width=ww,height=wh, paper="special",file="iav-fit.pdf"); 
  #dev.print(device=tiff,filename ="fig2.tif",width=ww, height=wh, units="in",pointsize = 12, compression =c("lzw"), res=300) 
  dev.print(device=png,width=ww,height=wh,units="in",res=600,file="../figures/fludrug.png")
  
}

