################################################
#Antibody Dependent Enhancement (ADE)
################################################

# This script calls and explores a simple toy model that is meant to 
# mimic antibody dependent enhancement.
# the uderlying function that runs the model is
# simulate_virusandir_new.R
# that function is a modification of the DSAIRM model
# simulate_virusandir_ode.R
# see the function for some more explanation of what's new



# load the function for our simulator model
# the function needs to be in the same folder as this script
# you also need to make sure that your working directory is the directory
# of the current script.
# to do so in Rstudio, go to "session" -> Set working directory -> Source file
source('simulate_virusandir_new.R')

#note that we don't need to load the DSAIRM package
#we are now using our own function and don't rely on 
#functionality from the DSAIRM package anymore


#remove earlier plots
graphics.off() 

#this parameter controls the strength of ADE
#we vary it here to explore how ADE might impact outcomes
b2vec = c(0,10^seq(-14,-6,length=30)) #strength of ADE

#we assume we are interested in Vpeak and dead cells
#those are thought of as proxies of infection severity
#we expect they might increase as the ADE effect in the model
#(parameter b2) increases

#make empty vectors to store outcomes
Vpeakall = rep(0, length = length(b2vec) ) 
Dall = rep(0, length = length(b2vec) ) 


#run model for different levels of ADE
for (nn in 1:length(b2vec))
{

  #turn on B-cells/Ab and turn off T-cells
  #also lower rate of infection
  #all other model inputs are kept at their defaults
  res <- simulate_virusandir_new(B = 1, A = 10, 
                                 b = 1e-7,
                                 gT=0,rT=0, 
                                 b2 = b2vec[nn])

  #for a better understanding of what's going on as we loop over b2
  #we plot the time-series each time
  
  plot(res$ts$time,res$ts$V,log='y',type = 'l',ylim=c(1e-2,1e6))
  lines(res$ts$time,res$ts$I,col="red")
  lines(res$ts$time,res$ts$F,col="blue")
  lines(res$ts$time,res$ts$A,col="green",lty=2)
  lines(res$ts$time,res$ts$B,col="magenta",lty=2)
  
  #The browser() command is a useful R command. It makes the code run to this point and then stop
  #you can then manually restart the code by pressinc 'c'
  #that means you can run the code for one value of b2, get the plot and stare at it.
  #once you understood the plot, you press c to go to the next value of b2
  #this is often better than having all plots zip by
  #uncomment this line to pause after each plot
  #browser()
  
  #this stores the final value for the dead cells (which is the total, since they accumulate in the D variable)
  #it also stores the peak of the virus
  Dall[nn] = tail(res$ts$D,1)
  Vpeakall[nn] = max(res$ts$V)  
  
}

# here we plot our outcomes of interest as function of b2 (i.e. ADE strength)
# divide plotting area into 2 parts
# we see that ADE leads to increased virus peak and more dead cells
# things jump from a regime of no infection to one where ADE enables infection to saturation
# of course these results are contingent on the values of all other parameters
# it might be worth exploring how changes in other model parameters impact these results/figures
par(mfrow=c(1,2))
plot(b2vec,Vpeakall,log='xy')
plot(b2vec,Dall,log='xy')



