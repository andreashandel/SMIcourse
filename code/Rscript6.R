################################################
#Antibody Dependent Enhancement (ADE) - version 2
################################################

# The previous script and underlying function didn't discriminate between 
# B-cells/antibodies for different serotypes
# To account for that, we built a model with two types of B-cells and antibodies
# that model is called simulate_virusandir_new2.R 
# and is an extension of simulate_virusandir_new.R
# see that function for some details on the extension.
# This script calls the new fuction to explore ADE, now taking 
# account of 2 different B-cell and antibody populations


# load the function for our simulator model
source('simulate_virusandir_new2.R')

#remove earlier plots
graphics.off() 

# this is similar to before, but now ADE, implemented by parameter b2
# is induced by antibodies of type 2, while type 1 antibodies/B-cells
# are specific for the new infection and don't cause ADE
b2vec = c(0,10^seq(-14,-6,length=30)) #strength of ADE

#we assume that we are again interested in Vpeak and dead cells
Vpeakall = rep(0, length = length(b2vec) ) 
Dall = rep(0, length = length(b2vec) ) 

#run model for different levels of ADE
for (nn in 1:length(b2vec))
{
  #turn on B-cells/Ab and turn off T-cells
  #notice that we have starting values for both antibodies and B-cells from the previous infection
  #(B2 and A2) and the new infection (A and B). The varying ADE (parameter b2)
  #is mediated by A2, not A.
  res <- simulate_virusandir_new2(B = 1, A = 10, B2 = 10, A2 = 100, 
                                 b = 1e-6, gB2 = 0,
                                 gT=0,rT=0, 
                                 b2 = b2vec[nn])

  # plot time-series after each run
  # now also showing new B-cell/Ab compartments, A2 and B2
  
  plot(res$ts$time,res$ts$V,log='y',type = 'l',ylim=c(1e-2,1e6))
  lines(res$ts$time,res$ts$I,col="red")
  lines(res$ts$time,res$ts$F,col="blue")
  lines(res$ts$time,res$ts$A,col="green",lty=2)
  lines(res$ts$time,res$ts$B,col="magenta",lty=2)
  lines(res$ts$time,res$ts$A2,col="green",lty=3)
  lines(res$ts$time,res$ts$B2,col="magenta",lty=3)
  
  Dall[nn] = tail(res$ts$D,1)
  Vpeakall[nn] = max(res$ts$V)  
  
}

#looking at results, i.e. virus peak and total dead cells as function of ADE
#it seems the extended model leads to similar results
#of course, further careful exploration is needed
#and if one wanted to apply this to a specific pathogen (e.g. Dengue) "for real"
#one would need to carefully match parameters and decide if the components included in the model
#and the various processes are good approximations of what is known about Dengue within-host infections and ADE
par(mfrow=c(1,2)) # divide plotting area into 2 parts
plot(b2vec,Vpeakall,log='xy')
plot(b2vec,Dall,log='xy')
