#code to create figures for course/book


#figure to distinguish between phenomenological and mechanistic model
library(DSAIDE)
set.seed(101)

Dvec = seq(3,10,by=0.1) #values of recovery rate, g, for which to run the simulation 
gvec = 1/Dvec
allinf = rep(0,length(gvec)) #this will record infected for each g
for (n in 1:length(gvec))
{
    #call the simulator function with different values of g each time
    result <- simulate_sir_ode(S = 500, I = 5, tfinal = 200, g = gvec[n],  b = 1/2500)
    allinf[n] <- tail(result$ts[,"R"],1) #record total number of infected for each value of g
}

infdata = allinf + rnorm(length(gvec),mean = 0, sd = allinf/5) 
fit=lm(allinf/500 ~ Dvec + I(Dvec^2) + I(Dvec^3))
predicted.intervals <- predict(fit,data.frame(x=Dvec),interval='confidence', level=0.99)


ww=17.8/2.54; wh=2/3*ww; #size is 17.8cm, needs to be in inches
windows(width=ww, height=wh) #for windows: opens window of the specified dimensions

plot(Dvec,infdata/500,type='p',xlab='Duration of infection, D',ylab='Fraction infected, F',ylim=c(0,1))
lines(Dvec,predicted.intervals[,1],col='red',lwd=3)
dev.print(device=png,width=ww,height=wh,units="in",res=600,file="phenoexample.png")


lines(Dvec, allinf/500,type='l',col = 'green',lwd=3)
legend("topleft",c("Data","Cubic","Model"), col=c("black","red","green"), lwd=3)

#save in various formats
#dev.print(device=pdf,width=ww,height=wh, paper="special",file="figinoculum.pdf");
#dev.print(device=tiff,filename ="../resultfiles/fig7.tif",width=ww, height=wh, units="in",pointsize = 12, compression =c("lzw"), res=300)
dev.print(device=png,width=ww,height=wh,units="in",res=600,file="SIRpredictexample.png")
