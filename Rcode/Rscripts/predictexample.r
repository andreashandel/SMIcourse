############################################################
##a simple model for the impact of drug on infection dynamics
#script produces figures for within-host intro chapter
##written by Andreas Handel (ahandel@uga.edu), last change 9/18/17
############################################################
rm(list=ls()) #this clears the workspace to make sure no leftover variables are floating around. Not strictly needed
graphics.off(); #close all graphics windows
library(deSolve)  #loads ODE solver package
library(caTools) #to do integral with trapz

#functions come first, main program below

###################################################################
#function that specificies the ode model called by lsoda (the ode solver) 
###################################################################
odeequations=function(t,y,parameters) 
{ 
  	Utc=y[1]; Itc=y[2]; Vir=y[3];  #uninfected target cells, infected target cells, virus
		b=parameters[1]; delta=parameters[2]; p=parameters[3]; clear=parameters[4]; #four model parameters, passed into function by main program
        n=parameters[5]; m=parameters[6]; bp=parameters[7]; e=parameters[8];
			 
	  #these are the 3 differential equations
		dUtcdt=-b*Vir*Utc-m*Utc+n;
		dItcdt=b*Vir*Utc-delta*Itc-m*Itc;
		dVirdt=p*(1-e)*Itc-clear*Vir-bp*Vir*Utc;
 
		return(list(c(dUtcdt,dItcdt,dVirdt))); 

} #end function specifying the ODEs

###################################################################
#main program
###################################################################
                
  Utc0=1e8; #initial number of uninfected cells 
	Itc0=0; #initial number of infected cells
	Vir0=10; #initial condition for free virus V
  Y0=c(Utc0, Itc0, Vir0);  #combine initial conditions into a vector 
  
  #values for model parameters, units are assumed to be 1/days
	b=1e-8;
	delta=2;
	p=1e2;
	clear=10; #Note: we don't call the parameter "c" since c is a reserved function in R that creates vectors (see next line for a use of it)
    n=Utc0;
    m=1;
    bp=b;
    nobs=51
    evec = seq(0,1,length=nobs)
    VAUC=rep(0,nobs)
    timevec=seq(0,30,0.1); #vector of times for which integration is evaluated (from 0 to 10 days in steps of 0.1)
    for (n in 1:nobs)
    {
        parameters=c(b,delta,p,clear,0,0,bp,e=evec[n]); #vector of parameters which is sent to the ODE function
        odeoutput=lsoda(Y0,timevec,odeequations,parameters);
        VAUC[n]=trapz(odeoutput[,1],odeoutput[,4])   
        #browser()
    }




graphics.off(); #close all graphics windows
ww=17.8/2.54; wh=1/2*ww; #size is 17.8cm, needs to be in inches
windows(width=ww, height=wh) #for windows: opens window of the specified dimensions

par(mfrow=c(1,1)) #margin at bottom, left, top, right

par(oma=c(0.2,0.2,0.2,0.2)) #outer margins
par(mar=c(0, 0, 0, 0)) #bottom, left, top, right margins
par(mai=c(0.7, 0.7, 0.1, 0.1)) #bottom, left, top, right margins


  #plot results
  plot(evec,VAUC,type="p",xlab="",ylab="",col="black",lwd=2,log="y",xlim=c(0,1),ylim=c(1,1e9))
 
  mtext("Drug efficacy, e",side=1,line=2.5) #x-axis - change line to move in/out
  mtext("Total virus load, VAUC",side=2,line=2) #y-axis
 
  dev.print(device=png,width=ww,height=wh,units="in",res=600,file="../figures/predictexample.png")
  
  
###################################################################
#end main program
###################################################################                           
