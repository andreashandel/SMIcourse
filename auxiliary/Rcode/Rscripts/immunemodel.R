##################################################################################
#uses the best fit model and runs simulations of the model for a range of inoculum doses
#also do a single run to get time-series for each inoculum for plotting
##################################################################################
rm(list=ls()) #this clears the workspace to make sure no leftover variables are floating around. Not strictly needed
graphics.off(); #close all graphics windows
library(deSolve)  #loads ODE solver package

modelequations=function(t,y,parms)
{
  with(
    as.list(c(y,parms)), #lets us access variables and parameters stored in y and pars by name
    {
      dF = 24;
      pF = dF;
      
      dUdt = - b*V*U;
      dIdt = b*V*U - dI*I - kT*T*I;
      dVdt = p/(1 + sF*F) * I - dV*V - kAp*A*V - bp*U*V;
      dFdt = pF - dF*F + V / (V + hV) * gF * (Fmax - F); #assume that F is at 1 at steady state in absence of virus   
      dBdt = F * V / (F * V + hF) * gB * B ; 
      dAdt = rA*B  - dA*A - kA*A*V;
      dTdt = a*V + gT*T
      
      list(c(dUdt, dIdt, dVdt, dFdt, dBdt, dAdt, dTdt))
    }
  ) #close with statement
} #end function withinode


eps=1e-12;
atolv=1e-12; rtolv=1e-12; #tolerances for ODE solver
U0=1e7; #number of cells

upos=1; ipos=2; vpos=3; fpos=4; bpos=5; apos=6;  #position/index for the different variables
tpos=7;

parvec = c(b = 1e-6, dI = 1, kT = 1e-8, p = 1e4, sF = 1e3, dV = 10, kAp = 1e-8, bp = 1e-6, hV = 1e4, gF = 1e-2, Fmax = 1e4, hF = 1e3, gB = 0.1, rA = 1e-3, dA = 1, kA = 1e-8, a = 1e-4, gT = 0.05)

tmax=21*24; #end of simulation in hours. 
t.int.vec=seq(0,tmax,by=0.1); #vector for integration times, 

Y0=c(U = U0, I = 0, V = 1, F = 1, B = 1, A = 0, T = 0)
odesol = ode(Y0,t.int.vec,modelequations,parms=parvec,atol=atolv,rtol=rtolv,method="vode"); #

graphics.off(); #close all graphics windows 
ww=17.8/2.54; wh=2/3*ww; #size is 17.8cm, needs to be in inches
windows(width=ww, height=wh) #for windows: opens window of the specified dimensions

#par(oma=c(0.2,0.2,0.2,0.2)) #outer margins 
par(mar=c(4, 3, 0.5, 0.2)) #bottom, left, top, right margins

plot(odesol[,1],odesol[,vpos+1],log="y",ylim=c(1,1e12),type="l",xlab='',ylab='')
lines(odesol[,1],odesol[,ipos+1],col='green')
lines(odesol[,1],odesol[,fpos+1],col='orange')
lines(odesol[,1],odesol[,bpos+1],col='blue')
lines(odesol[,1],odesol[,apos+1],col='red')
lines(odesol[,1],odesol[,upos+1],col='cyan') 
lines(odesol[,1],odesol[,tpos+1],col='gray') 
legend("topleft",legend=c('U','I','V','F','B','A','T'), col=c('cyan','green','black','orange','blue','red','gray'),lwd=2)

mtext("Hours post infection",side=1,line=2.5) #x-axis - change line to move in/out
#mtext("B cells",side=2,line=2) #y-axis

#save in various formats
#dev.print(device=pdf,width=ww,height=wh, paper="special",file="iav-fit.pdf"); 
#dev.print(device=tiff,filename ="fig2.tif",width=ww, height=wh, units="in",pointsize = 12, compression =c("lzw"), res=300) 
dev.print(device=png,width=ww,height=wh,units="in",res=600,file="../figures/immunemodel.png")

       
