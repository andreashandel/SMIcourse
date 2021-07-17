# This is short R script which illustrates how one can
# use code to call the DSAIRM simulator functions and process the results,
# without going through the graphical interface.

# load the package
library(DSAIRM)

# run the basic bacteria model with default settings
result <- simulate_basicbacteria_ode()

#inspect structure of object returned from the simulation function
str(result)

#look at the first few rows of the time-series (ts) data frame
head(result$ts)

#assign the time-series data frame that's inside the result list
#to a new variable for easier typing
ts = result$ts

#make a plot for B and I
#using the numeric values returned from the simulator
plot(ts[,"time"],ts[,"B"],type = "l")
lines(ts[,"time"],ts[,"I"],col='red')

#this uses a built-in DSAIRM function to produce a plot
#this plot looks like the one seen in the graphical interface
generate_ggplot(list(result))

# run model with some custom settings
result <- simulate_basicbacteria_ode(B = 100, I = 1, g = 2, 
                                     Bmax = 1e+05, dB = 0.5, 
                                     k = 1e-4, r = 1e-4, dI = 1, 
                                     tstart = 0, tfinal = 100, dt = 0.05)

# look at plot for custom settings model run
generate_ggplot(list(result))

#get max and final values of B
max(ts$B)
tail(ts$B,1)

#get time at which maximum happens
ind = which.max(ts$B)
ts[ind,"time"]

# run the discrete time version of the basic bacteria model
# same settings as ODE model
result2 <- simulate_basicbacteria_discrete(B = 100, I = 1, g = 2, 
                                           Bmax = 1e+05, dB = 0.5, 
                                           k = 1e-4, r = 1e-4, dI = 1, 
                                           tstart = 0, tfinal = 100, dt = 0.05)

#plot ODE and discrete time results to compare
# ODe model is solid, discrete time is dashed
#for the time step chosen here, they are not the same
plot(result$ts[,"time"],result$ts[,"B"],type = "l")
lines(result$ts[,"time"],result$ts[,"I"],col='red')
lines(result2$ts[,"time"],result2$ts[,"B"],col='black',lty=2)
lines(result2$ts[,"time"],result2$ts[,"I"],col='red',lty=2)


# this shows how one can apply the same code and call another one of the 
# underlying simulators, in this case the basic virus model
vir_res <- simulate_basicvirus_ode(U = 1e+05, I = 0, V = 10,
  n = 0, dU = 0, dI = 1, dV = 4, b = 1e-06, p = 100,  g = 1,
  tstart = 0, tfinal = 50, dt = 0.1)

#look at return structure
#this looks similar, but now has 4 columns in ts, 
#1 for time and 3 for the 3 model variables
str(vir_res)

# assign data frame to a new variable
ts_vir = vir_res$ts 

# make a plot
plot(ts_vir[,"time"],ts_vir[,"U"],type = "l",log = "y", ylim = c(1,1e7) )
lines(ts_vir[,"time"],ts_vir[,"I"],col='red')
lines(ts_vir[,"time"],ts_vir[,"V"],col='blue')

#this shows a quick example of using the output and making a ggplot
#load some other packages for data processing and plotting
library(tidyr)
library(ggplot2)

#ggplot needs the data in a certain format
#this line reformats the data to get it into a shape that ggplot wants
#note that the PACKAGE::FUNCTION notation is not required
#after loading the package, one could just write pivot_longer without the tidyr::
#I'm adding it here so it's clear that the pivot_longer function 
#comes from the tidyr package
df <- tidyr::pivot_longer(ts_vir[1:4], 2:4, names_to = "variable", values_to = "value")

# look at structure of new data frame
head(df)

#make a basic ggplot
pl <- ggplot2::ggplot(data = df, aes(x = time, y = value, group = variable, col = variable)) + geom_line() + scale_y_log10()
plot(pl)





