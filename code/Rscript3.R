# This script shows code for interacting with the
# uncertainty and sensitivity analysis app/model

#load pacakge
library('DSAIRM')

#run the uncertainty/sensitivity analysis simulation
sim_result <- simulate_usanalysis(Bmin = 100, Bmax = 200, Imin = 1, Imax = 2,
                                  Bmaxmin=1e5, Bmaxmax=2e5,
                                  dBmin=0.5, dBmax = 1, kmin=1e-4, kmax=2e-4,
                                  rmin=1e-4, rmax=2e-4, dImin=1, dImax=2,
                                  gmean=2, gvar=0.5,
                                  samples = 50, rngseed = 100,
                                  tstart = 0, tfinal = 300, dt = 0.05)

# look at output
# this output again does not contain time-series
# as for the model exploration function, the US analysis function
# runs the underlying model, then extracts quantities of interest
# and only returns those
str(sim_result)

#make boxplots using base plot commands
#those are similar to the ones shown in the graphical interface
par(mfrow=c(1,3)) #tells R to divide plot window into 3 columns
boxplot(sim_result$dat$Isteady,xlab = "Isteady")
boxplot(sim_result$dat$Bsteady,xlab = "Bsteady")
boxplot(sim_result$dat$Bpeak,xlab = "Bpeak")

#make a scatterplot for several of the parameters and outcomes
#there are scatterplots to look at for every outcome-parameter combination
#for this example, 3 outcomes and 8 sampled parameters = 24 scatter plots
par(mfrow=c(3,2))
plot(x = sim_result$dat$g, y = sim_result$dat$Isteady, type = "p")
plot(x = sim_result$dat$B, y = sim_result$dat$Isteady, type = "p")
plot(x = sim_result$dat$g, y = sim_result$dat$Bsteady, type = "p")
plot(x = sim_result$dat$B, y = sim_result$dat$Bsteady, type = "p")
plot(x = sim_result$dat$g, y = sim_result$dat$Bpeak, type = "p")
plot(x = sim_result$dat$B, y = sim_result$dat$Bpeak, type = "p")


# the following bit of code illustrates the impact of random number seeds
# if you run the simulation with a fixed number in myseed (first line)
# the resulting plots will look the same every time you run the code.
# this is because random numbers are not truly random, they are reproducible.
# If you run the simulation using Sys.time(), you get a different random seed
# each time. The result is that every run looks different.
# try it out, run the model and plot results twice with the same fixed seed
# then repeat with the varying seed.
myseed = 100 #fixed random seed
#myseed = as.numeric(Sys.time()) #changing random seed

#run again with more samples
sim_result <- simulate_usanalysis(Bmin = 100, Bmax = 200, Imin = 1, Imax = 2,
                                  Bmaxmin=1e5, Bmaxmax=2e5,
                                  dBmin=0.5, dBmax = 1, kmin=1e-4, kmax=2e-4,
                                  rmin=1e-4, rmax=2e-4, dImin=1, dImax=2,
                                  gmean=2, gvar=0.5,
                                  samples = 50, rngseed = myseed,
                                  tstart = 0, tfinal = 300, dt = 0.05)

#make a scatterplot for one of the parameters and outcomes
#this will not change if the model is re-run with the same random seed
#if the random seed is changed, the result will change
par(mfrow=c(1,1)) #one large plot window
plot(x = sim_result$dat$g, y = sim_result$dat$Isteady, type = "p")


#################################################################################
# this example shows how one can combine a systematic scanning over some parameter with sampling for others
# we run a loop over some parameter we are most interested in (here, we pick dI)
# for each value of this parameter, we sample the others
# now instead of getting just one value per main parameter (as we did with the model exploration functions)
# we get distributions of outcomes for each value of the main parameter.
#################################################################################

dIvec = seq(1,5,length=10) #looping over dI
psamples = 20 #number of random samples for all parameters

#we assume we are interested in Bpeak
#you could pick any outcome you want
Bpeakall = array(0, dim = c(length(dIvec),psamples) ) 

#we loop over dI
for (n in 1:length(dIvec))
{

  #for each value of dI, we sample all other parameters
  #we prevent dI from being sampled, and instead make sure it's the fixed value we want
  #by setting both dImin and dImax to the same value.
  #that means dI is sampled, but each sample is exactly the same value, the one we want.
  res <- simulate_usanalysis(Bmin = 100, Bmax = 200, Imin = 1, Imax = 2,
                                    Bmaxmin=1e5, Bmaxmax=2e5,
                                    dBmin=0.5, dBmax = 1, kmin=1e-4, kmax=2e-4,
                                    rmin=1e-4, rmax=2e-4, dImin=dIvec[n], dImax=dIvec[n],
                                    gmean=2, gvar=0.5,
                                    samples = psamples, rngseed = myseed,
                                    tstart = 0, tfinal = 300, dt = 0.05)

  #for each value of dI, we store results for Bpeak
  #based on sampling all the other parameters
  Bpeakall[n,] = res$dat$Bpeak

}

#this series of boxplots corresponds to the following scientific question:
#How does some outcome of interest (here peak of bacteria) vary with some parameter of interest (here dI)
#and how uncertain are the outcomes given that there is uncertainty in the values of the other model parameters
boxplot(t(Bpeakall))




