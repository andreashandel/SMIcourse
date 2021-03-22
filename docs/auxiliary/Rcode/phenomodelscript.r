############################################################
#script produces a simple figure showing hypothetical drug effect
#for within-host intro chapter
##written by Andreas Handel (ahandel@uga.edu), last change 10/3/17
############################################################
rm(list=ls()) #this clears the workspace to make sure no leftover variables are floating around. Not strictly needed

###################################################################
#main program
###################################################################

nobs=20
concvec=10^seq(-4,4,length=nobs)
chalf=1;
Vmax=1e6;
VAUC=Vmax*(1-concvec/(concvec+chalf))


graphics.off(); #close all graphics windows
ww=4; wh=1/1*ww; #plotting window width and height in inches
windows(width=ww, height=wh) #for windows: opens window of the specified dimensions

par(mfrow=c(1,1)) #margin at bottom, left, top, right

par(oma=c(0.2,0.2,0.2,0.2)) #outer margins
par(mar=c(0, 0, 0, 0)) #bottom, left, top, right margins
par(mai=c(0.7, 0.7, 0.1, 0.1)) #bottom, left, top, right margins


  #plot results
  #first column contains time vector, the following columns contain variables 1 (bacteria) and 2 (immune response)
  plot(concvec,VAUC+0.5*runif(nobs,min=-VAUC ,max=VAUC),type="p",xlab="",ylab="",col="blue",lwd=2,log="xy",xlim=c(min(concvec),max(concvec)),ylim=c(1e2,1e6))
  lines(concvec,VAUC,col="blue",lwd=2)
  
  #legend("topleft", c("Bacteria","Immune Response"),col = c("blue","red"),lwd=2,lty=c(1,2))
  
  mtext("Drug concentration",side=1,line=2.5) #x-axis - change line to move in/out
  mtext("Total viral load",side=2,line=2) #y-axis
  
  dev.print(device=png,width=ww,height=wh,units="in",res=600,file="../figures/phenomodel.png")
  
  
  
###################################################################
#end main program
###################################################################                           
