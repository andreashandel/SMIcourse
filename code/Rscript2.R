# R Code that shows some additional ways to interact with the underlying DSAIRM simulation functions
# here we look at the basic bacteria model exploration app/function 
# this code also shows a for loop that can be used to scan over parameters

#load the package
library('DSAIRM')

#call the function which runs the basic bacteria model
#and does so for several settings of a chosen parameter
#specifically, you set the samples you want, the lower and upper value, and the parameter to sample
#see the documentation inside the DSAIRM graphical interface for more details
res <- simulate_basicbacteria_modelexploration( B = 100, I = 10,
                                                g = 2,  Bmax = 1e+05,  dB = 1,  k = 1e-04,  r = 1e-04,
                                                dI = 2,  tstart = 0,  tfinal = 300,  dt = 0.1,
                                                samples = 10,  parmin = 2,  parmax = 10,  samplepar = "g", pardist = "lin"
                                                )

#looking at return structure
#this is not a time-series anymore
#instead, the simulator calls the basic bacteria ODE model for each parameter
#then takes the returned result and extracts peak and steady state
#(similar to the code shown in the previous R script)
str(res)

#plotting one of the outcomes against the varied parameter (here g)
#this plot looks similar to the one you see in the graphical interface
plot(res$dat[,"xvals"],res$dat[,"Bpeak"],type = "p")

#if you interact with the simulator functions through code, you have more flexibility
#here is an example where we scan over 2 parameters by writing our own for-loop

#we want to scan over the parameter dI, in addition to another one
#so we create values that we want to run the model at
dIvec = seq(1,10,length=10) #looping over dI

#number of samples for the other parameter
#defined here so we can easier adjust
psamples = 20 #number of samples for other parameter

#set up empty arrays that will contain our results of interest
#here we look at peak for B and I as function of all combinations of the 
#2 parameters of interest
Bpeakall = array(0, dim = c(length(dIvec),psamples) ) #we assume we are interested in Bpeak
Ipeakall = array(0, dim = c(length(dIvec),psamples) ) #we assume we are interested in Bpeak


# this for-loop loops over the dI parameter
# then for each value of that parameter, we let the simulator function loop over the other

for (n in 1:length(dIvec))
{
  #run the simulator for different values stored in dIvec
  #here, we tell the function to loop over g as well
  #you could change this by changing the entry for samplepar
  res <- simulate_basicbacteria_modelexploration( B = 100, I = 10,
                                                  g = 2,  Bmax = 1e+05,  dB = 1,  k = 1e-04,  r = 1e-04,
                                                  dI = dIvec[n],  tstart = 0,  tfinal = 300,  dt = 0.1,
                                                  samples = psamples,  parmin = 2,  parmax = 20,  samplepar = "g", pardist = "lin"
  )
  #for each value of dI (outer loop), store all values of Bpeak and Ipeak 
  #for the other parameter, which is done by the function we just called
  Bpeakall[n,] = res$dat$Bpeak
  Ipeakall[n,] = res$dat$Ipeak

}

#plot results as a simple 2d plot
image(dIvec,seq(2,10,length=psamples),Bpeakall)
image(dIvec,seq(2,10,length=psamples),Ipeakall)

#############################
#2 loops by hand
#############################
#note that above we used the model exploration function to do one of the loops for us
#that's convenient, but if we wanted to loop over other models that don't have an exploration function, this wouldn't work
#here is a quick example of doing 2 loops to scan over parameters for an underlying function that isn't already 
#pre-coded with an model exploration function, namely the virus and immune response model.

#parameters we want to sample
#choosing only a few samples so code doesn't run too long
dIvec=seq(0.5,1.5,length=5)
dVvec=seq(1,5,length=10)

#for each run, we'll record the maximum of the virus
#of course it's completely up to you what outcome you want to track
Vpeakall = array(0, dim = c(length(dIvec),length(dVvec)) ) 

#outer loop over dI, inner loop over dV
#the order doesn't matter
for (n1 in 1:length(dIvec))
{
  #for each value of dI, we run over all dV values
  for (n2 in 1:length(dVvec))
  {
    #run the model with default settings, but different values
    #of infected cell and virus death rates each time
    res <- simulate_virusandir_ode(dI = dIvec[n1], dV = dVvec[n2])
    
    #store peak of virus for each dI/dV combination
    Vpeakall[n1,n2] = max(res$ts$V)    

  } #closes inner loop

} #closes outer loop

#shows heat map of virus peak as function of dI and dV
image(dIvec,dVvec,Vpeakall)
