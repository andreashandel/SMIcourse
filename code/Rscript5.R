################################################
#Pre-existing antibody levels and clearance
################################################
# This uses the same function we built for ADE (see Rscript4)
# Here, we want to address a different question
# We are asking how different levels of initial antibodies and B-cells
# impact the peak and duration of the infection
# We won't use ADE, so that parameter is 0

# load the function for our simulator model
source('simulate_virusandir_new.R')

#remove earlier plots
graphics.off() 

#vectors for different values of initial antibodies and B-cells
#if you want to explore how different values in only one quantity impact outcomes
#you can simply make a vector of values for the other quantity that's contains
#all the same values, as shown 
parsamples = 30 #number of samples/values to loop over

Avec = 10^seq(5,9,length=parsamples) #levels of antibodies
Bvec = 10^seq(1,3,length=parsamples) #levels of B cells

#if you want to not have varying antibodies or B-cells just set start and end
#in the seq() command to the same value and uncomment this line
#Avec = c(0,10^seq(3,3,length=parsamples)) #levels of antibodies - same value of 10^3 repeated
#Bvec = c(0,10^seq(1,1,length=parsamples)) #levels of B cells


#we assume we are interested in Vpeak and time of clearance
#this sets up empty vectors to store it
Vpeakall = rep(0, length = parsamples ) 
Duration = rep(0, length = parsamples )


#run model for different levels of antibodies and B-cells
for (nn in 1:parsamples)
{

  #turn on B-cells/Ab and turn off T-cells
  #also lower rate of infection
  #for each loop (value of nn), pull out a value for B and A from Bvec and Avec
  #and run the model
  res <- simulate_virusandir_new(B = Bvec[nn], A = Avec[nn], 
                                 gB = 2,
                                 gT=0,rT=0)

  #for a better understanding of what's going on as we loop over b2
  #we plot the time-series each time
  
  plot(res$ts$time,res$ts$V,log='y',type = 'l',ylim=c(1e-2,1e6))
  lines(res$ts$time,res$ts$I,col="red")
  lines(res$ts$time,res$ts$F,col="blue")
  lines(res$ts$time,res$ts$A,col="green",lty=2)
  lines(res$ts$time,res$ts$B,col="magenta",lty=2)
  
  #you could add the browser() statement here, see previous R script.
  
  #getting peak of virus for each value of B and A
  Vpeakall[nn] = max(res$ts$V)  
  
  #find time at which virus drops below 1
  ind = min(which(res$ts$V<1)) #index for which it is below 1 for the first time
  Duration[nn] = res$ts[ind,"time"] #time at which that drop below 1 happens
}


#plot the results
#this answers our question: 
#How does increases in antibodies and/or B-cells impact the peak and duration
#of the infection
#The duration of infection shows an interesting pattern of down, then up then down to zero
#see the recording for a discussion of a potential reason
#or explore yourself :)

# divide plotting area into 2 parts
par(mfrow=c(1,2))
plot(Avec,Vpeakall,log='xy')
plot(Avec,Duration,log='x')
