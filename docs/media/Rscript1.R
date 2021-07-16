# load the package
library(DSAIRM)

# run model with default settings
result <- simulate_basicbacteria_ode()

#inspect structure of object returned from the simulation function
str(result)

#look at the first few rows
head(result$ts)

#assign the time-series data frame inside results
#to a new variable for easier typing
ts = result$ts

#make a plot for B and I
#using the numeric values returned from the simulator
plot(ts[,"time"],ts[,"B"],type = "l")
lines(ts[,"time"],ts[,"I"],col='red')

#this uses a built-in DSAIRM function to produce a plot
#this plot looks like the one seen in the UI
generate_ggplot(list(result))

# run model with some custom settings
result <- simulate_basicbacteria_ode(B = 100, I = 1, g = 2, 
                                     Bmax = 1e+05, dB = 0.5, 
                                     k = 1e-4, r = 1e-4, dI = 1, 
                                     tstart = 0, tfinal = 100, dt = 0.01)

# look at plot for custom settings model run
generate_ggplot(list(result))

#get max and final values of B
max(ts$B)
tail(ts$B,1)

#get time at which maximum happens
ind = which.max(ts$B)
ts[ind,"time"]


# run the discrete time version of the basic bacteria model
result <- simulate_basicbacteria_discrete(  B = 10,  I = 1,  g = 1,
  Bmax = 1e+06,  dB = 0.1,  k = 1e-07,  r = 0.001,  dI = 1,  tstart = 0,  tfinal = 30,  dt = 0.01)

# you can do anything to the result you want (e.g. plot or pull out specific values)


# run the basic virus model
vir_res <- simulate_basicvirus_ode(U = 1e+05, I = 0, V = 10,
  n = 0, dU = 0, dI = 1, dV = 4, b = 1e-06, p = 100,  g = 1,
  tstart = 0, tfinal = 50, dt = 0.1)

#look at return structure
str(vir_res)

# assign data frame to a new variable
ts = vir_res$ts 

# make a plot
plot(ts[,"time"],ts[,"U"],type = "l",log = "y", ylim = c(1,1e7) )
lines(ts[,"time"],ts[,"I"],col='red')
lines(ts[,"time"],ts[,"V"],col='blue')

#load some other packages for data processing and plotting
library(tidyr)
library(ggplot2)

#reformat the data to get it into a shape that ggplot wants
df <- tidyr::pivot_longer(ts[1:4], 2:4, names_to = "variable", values_to = "value")

# look at structure of new data frame
head(df)

#make a ggplot
pl <- ggplot(data = df, aes(x = time, y = value, group = variable, col = variable)) + geom_line() + scale_y_log10()
plot(pl)





