############################################################
##a simple model for an acute or chronic virus infection
#script produces time-series figures for within-host intro chapter
##written by Andreas Handel (ahandel@uga.edu), last change 9/18/17
############################################################
rm(list=ls()) #this clears the workspace to make sure no leftover variables are floating around. Not strictly needed
graphics.off(); #close all graphics windows
library(deSolve)  #loads ODE solver package

#functions come first, main program below

###################################################################
#function that specificies the ode model called by lsoda (the ode solver) 
###################################################################
odeequations=function(t,y,parameters) 
{ 
  	Utc=y[1]; Itc=y[2]; Vir=y[3];  #uninfected target cells, infected target cells, virus
		b=parameters[1]; delta=parameters[2]; p=parameters[3]; clear=parameters[4]; #four model parameters, passed into function by main program
        n=parameters[5]; m=parameters[6]; bp=parameters[7];
			 
	  #these are the 3 differential equations
		dUtcdt=-b*Vir*Utc-m*Utc+n;
		dItcdt=b*Vir*Utc-delta*Itc-m*Itc;
		dVirdt=p*Itc-clear*Vir-bp*Vir*Utc;
 
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
	#always make sure that your variable names do not clash with built-in definitions, or you might get unpredictable results
    parameters1=c(b,delta,p,clear,0,0,bp); #vector of parameters which is sent to the ODE function
    parameters2=c(b,delta,p,clear,n,m,bp); #vector of parameters which is sent to the ODE function

	timevec1=seq(0,30,0.1); #vector of times for which integration is evaluated (from 0 to 10 days in steps of 0.1)
    timevec2=seq(0,30,0.1); #vector of times for which integration is evaluated (from 0 to 10 days in steps of 0.1)

	#call ode-solver to integrate ODEs 
  odeoutput=lsoda(Y0,timevec1,odeequations,parameters1);
odeoutput2=lsoda(Y0,timevec2,odeequations,parameters2);



graphics.off(); #close all graphics windows
ww=17.8/2.54; wh=1/2*ww; #size is 17.8cm, needs to be in inches
windows(width=ww, height=wh) #for windows: opens window of the specified dimensions

par(mfrow=c(1,2)) #margin at bottom, left, top, right

par(oma=c(0.2,0.2,0.2,0.2)) #outer margins
par(mar=c(0, 0, 0, 0)) #bottom, left, top, right margins
par(mai=c(0.7, 0.7, 0.1, 0.1)) #bottom, left, top, right margins


  #plot results
  plot(odeoutput[,1],odeoutput[,2],type="l",xlab="",ylab="",col="green",lwd=2,log="y",xlim=c(0,max(timevec1)),ylim=c(1,1e9))
  lines(odeoutput[,1],odeoutput[,3],type="l",col="red",lwd=2,lty=2)
  lines(odeoutput[,1],odeoutput[,4],type="l",col="blue",lwd=2,lty=3)
  legend("topright", c("Uninfected","Infected","Virus"),col = c("green","red","blue"),lwd=2,lty=c(1,2,3))

  mtext("Days post infection",side=1,line=2.5) #x-axis - change line to move in/out
  mtext("Virus and cell numbers",side=2,line=2) #y-axis
  
  
  plot(odeoutput2[,1],odeoutput2[,2],type="l",xlab="",ylab="",col="green",lwd=2,log="y",xlim=c(0,max(timevec2)),ylim=c(1,1e9))
  lines(odeoutput2[,1],odeoutput2[,3],type="l",col="red",lwd=2,lty=2)
  lines(odeoutput2[,1],odeoutput2[,4],type="l",col="blue",lwd=2,lty=3)
  mtext("Days post infection",side=1,line=2.5) #x-axis - change line to move in/out
  
    
  dev.print(device=png,width=ww,height=wh,units="in",res=600,file="../figures/virusmodels.png")
  
  
###################################################################
#end main program
###################################################################                           
