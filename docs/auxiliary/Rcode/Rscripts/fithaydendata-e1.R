##################################################################################
#fitting data from IAV infections at different inoculum doses
#written by Andreas Handel (ahandel@uga.edu)
#last changed on 10/28/2017 
##################################################################################
rm(list=ls()) #this clears the workspace to make sure no leftover variables are floating around. Not strictly needed
graphics.off(); #close all graphics windows
library(deSolve)  #loads ODE solver package
library(nloptr) #for fitting
library(dplyr) #for data prep
library(snow) #for parallel computing


odeequations=function(t,y,parms,txtime)
{
    with(
        as.list(c(y,parms)), #lets us access variables and parameters stored in y and pars by name
        {
             
            
            if (t<txtime) {e1=0}
            
            dUdt = - b*(1-e1)*V*U;
            dIdt = b*(1-e1)*V*U - dI*I;
            dVdt = p * I - dV*V - bp*(1-e1)*U*V;
              
            #browser()
            
            list(c(dUdt, dIdt, dVdt))
        }
    ) #close with statement
} #end function withinode


###################################################################
#function that fits the ODE model to data
###################################################################
fitfunction <- function(pars.all,parnames)
{
    names(pars.all)=parnames; #nloptr strips names from parameters, need to re-assign
    pars.ode=10^(pars.all); #transform parameters back to exponential since we fit in log space

    SSR_virus=rep(0,nmax);
    allvirus_data = NULL;
    allvirus_model = NULL;

    txstatusvec=c('earlytx','latetx','notx')
    txtimevec=c(29/24,50/24,500/24) #time of treatment start - last value means never
    
    for (nn in 1:nmax) #loop over treatment conditions
    {
       
        txstatus = txstatusvec[nn]
        txtime = txtimevec[nn]
        timedata <- mydata %>% filter(Condition==txstatus) %>% select(DaysPI) %>% unlist()  %>% as.numeric()
        virusdata <- mydata %>% filter(Condition==txstatus) %>% select(LogVirusLoad) %>% unlist()  %>% as.numeric()

        allvirus_data = c(allvirus_data, virusdata)        

        t.int.vec=seq(0,max(timedata),by=0.01); #vector for integration times - units of hours
        Y0=c(U=U0, I=0, V=as.numeric(pars.ode[1])) #starting conditions
        
        odestack=try(ode(Y0,t.int.vec,odeequations,parms=pars.ode,atol=atolv,rtol=rtolv,method=odemethod,txtime=txtime)); #runs the ODE equations
        #odestack=try(lsoda(Y0,t.int.vec,iavequations,parms=pars.ode,atol=atolv,rtol=rtolv)); #runs the ODE equations
        if (length(odestack)==1) {cat('!!unresolvable integrator error - triggering early return from optimizer!!'); return(1e10) } #catching errors that might happen during fitting
        vir.model.lin=odestack[match(timedata,odestack[,1]),vpos+1]; #extract values for virus load at time points corresponding to experimental measurements
        vir.model=log10(pmax(eps,vir.model.lin)); #convert to log scale. prevent potential negative entries by setting them to a small number (eps).
        #for censored data, set difference to zero if model prediction is lower than censored value/level of detection
        vir.model[(virusdata<=Vmin & (virusdata-vir.model)>0)]=Vmin

        #save results from model for virus load and cell damage/death for each inoculum dose to compute SSR below
        allvirus_model = c(allvirus_model, vir.model)        

        if (ploton==1) #plotting, for diagnostic/debugging purposes only
        {
            if (nn==1)
            {
                plot(timedata,10^virusdata,type='p',lty=1,lwd=2,ylab="",xlab="hours p.i.",log="y",xlim=c(0,10),ylim=c(1e-1,1e9),pch=nn,col=nn);
            }
            points(timedata,10^virusdata,lty=1,pch=19,col=nn);
            lines(odestack[,1],odestack[,vpos+1],col=nn,lty=1,lwd=2)  #virus
            lines(odestack[,1],odestack[,upos+1],col=nn,lty=2,lwd=2)  #uninfected
            lines(odestack[,1],odestack[,ipos+1],col=nn,lty=3,lwd=2)  #infected
            
              browser()
        } #end plotting part
    } #end loop over all inoculum doses

    #compute SSR for virus and damage. 
    #Extra scaling for damage since dead cells in model 
    #might not be a 1:1 mapping to experimentally reported lung damage 
    SSR.virus = sum( (allvirus_data-allvirus_model)^2 )
    
    #objective function that is optimized
    Fobject=SSR.virus;
    if (is.na(Fobject)) {Fobject=1e10}; #if something went wrong, e.g. because ODE solver messed up and returned nonsense leading to NA, we set Fobject to a really large number, otherwise optimizer might stop
    return(sum(Fobject))
} #end of fitting function

#############################################################
#wrapper function that runs parallel over different methods
#############################################################
outerfitfc <- function(solvertype) 
{

    #check if we already have previous best fit estimates
    scind=nrow(previous.res); 
    
    if (length(scind)>0) #if previous results exist, use them, otherwise create random starting guess
    { 
        minind = which.min(previous.res[,"SSR"])
        p0 = log10(previous.res[minind, match(parnames,names(previous.res[minind,]))]);
    }
    else 
    {
        p0 = (logub+loglb)/2
    }  

    #p0=log10(c(V1=1e-1, b=1e-8,bp=1e-8,dI=2,dV=10,p=5e1,e1=0.5))
      
    #browser()
    names(p0)=parnames; #assign names to parameters

    #error check
    if (min(logub-p0)<0 | min(p0-loglb)<0) { print(sprintf('initial condition out of bound')); p0=pmax(pmin(logub,p0),loglb);  }  

    #nloptr global solvers - setting solver to specific numbers below lets one run different solvers
    if (solvertype==11)  {alg.name="NLOPT_GN_DIRECT";    }  #solver generally performs poorly
    if (solvertype==12) {alg.name="NLOPT_GN_CRS2_LM";     }
    if (solvertype==13) {  alg.name="NLOPT_GN_ISRES";     }
    if (solvertype==14) { alg.name="NLOPT_GN_MLSL_LDS";   }
   
    #nloptr local gradient-free solvers
    if (solvertype==21)      {         alg.name="NLOPT_LN_SBPLX"       }
    if (solvertype==22)         {        alg.name="NLOPT_LN_COBYLA"     }
    if (solvertype==23)        {        alg.name="NLOPT_LN_BOBYQA"      }     #constructs a quadratic approximation of the objective
    if (solvertype==24)      {         alg.name="NLOPT_LN_PRAXIS"        }
    if (solvertype==25)      {        alg.name="NLOPT_LN_NELDERMEAD";      }

    if (solvertype!=14) {fres <- nloptr(x0=p0,eval_f=fitfunction,lb=loglb,ub=logub,opts=list("algorithm"=alg.name,xtol_rel=xerr.general,maxeval=max.stps,maxtime=max.time,print_level=plevel),parnames=parnames)}
    if (solvertype==14) {fres <- nloptr(x0=p0,eval_f=fitfunction,lb=loglb,ub=logub,opts=list("algorithm"=alg.name,xtol_rel=xerr.global, maxeval=max.stps,maxtime=max.time,print_level=plevel,"local_opts"=list("algorithm"="NLOPT_LN_SBPLX",xtol_rel=xerr.local,maxeval=max.stps.local)),parnames=parnames) }

    fpar=10^(fres$solution); Obj.final=fres$objective; solverstatus=fres$status;
    names(fpar) <- parnames
    
    if (length(scind)>0) #if previous fit exists, compute improvement in fit
    {
        improvement=(previous.res[minind,2]-Obj.final)/previous.res[minind,2]*100
    }
    else 
    {
        improvement=0
    }
    
    print(sprintf('Finished  using solver %s. SSR %f', alg.name, Obj.final));
    print(sprintf('Steps taken %d/%d. Solver status is %d (between 1-4 is good, see nlopt webpage for details)',fres$iterations,max.stps,fres$status));
    print(sprintf('Obj.Fct=%e; Percent improvement in fit: %f',Obj.final,improvement))
    
    resvector=c(solvertype,Obj.final,improvement,solverstatus,fpar)
    names(resvector) <-  matnames
    
    return(resvector) #return results from optimizer to main function
    
} #finish function for each scenario

#################################
#main program
#################################
tstart=proc.time(); #capture current time to measure duration of process
eps=1; #cut-off such that anything below 1 virion is scored as 0 on a log scale

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#data from Hayden 1996 JAMA
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mydata=read.csv('hayden96data.csv') #data file needs to be in same directory as this script

nmax = 3; #number of conditions/simultaneous fits

nobs = nrow(mydata) #number of data points

U0=4e8; #number of cells - estimate/guess  
#Vmin=0.5; #level of detection for virus in log10 units according to paper
Vmin=0; #use this level since data are averages 

#done with data part
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#bounds on initial conditions
#bounds for parameters improves solver convergence
V1low=1e-1; V1high=10^2; #for V0

#bounds on parameters, all is in units of  HOURS
blow=1e-14; bhigh=1e1; #cell infection rate
bplow=1e-16; bphigh=1e1; #cell infection rate
plow=1e-10; phigh=1e10; #virus production rate

#rate of death/clearance for infected cells, innate IR, virus, B-cells, Ab
dIlow= 1e-3; dIhigh= 1e2;  #rate of death of infected cells 
dVlow=1e-3; dVhigh=1e4; #rate of virus clearance         
e1low=1e-15; e1high=1;

#mapping of experimental cell damage to dead cells
#dmaplow=1e-2; dmaphigh=1e2;


#vector of parameter names and lower/upper bounds
parnames=c('V1', 'b','bp',  'dI', 'dV', 'p',  'e1' );
lb=c(V1low, blow,bplow, dIlow,dVlow, plow, e1low);
ub=c(V1high, bhigh,bphigh,dIhigh,dVhigh,phigh,  e1high);

loglb=log10(lb); logub=log10(ub); #fit in log space


#load previous best fit and save new best fit
filename.load=paste("hayden-fit-results-e1.Rdata",sep="");
filename.save=paste("hayden-fit-results-e1.Rdata",sep="");


resmatrix=NULL;
load(filename.load); #loads previous results (saved in variable resmatrix) and uses them as starting values. If not available, resmatrix remains NULL
previous.res=resmatrix;

upos=1; ipos=2; vpos=3; #position/index for the different variables

atolv=1e-12; rtolv=1e-12; #tolerances for ODE solver
odemethod="vode"

#solver types to run
#can choose multiple solvers (will be run in parallel) or a single one
#solvertype = c(21,22,23,24,25);
#solvertype = c(12,14,21,25);
solvertype = c(12,13,14,21,22,23,24,25);
#solvertype = c(21);

#matrix containing all results
matnames=c('Solver','SSR','Improvement','Solverstatus',parnames)
resmatrix=matrix(0,nrow=length(solvertype),ncol=length(matnames))
colnames(resmatrix) <- matnames;

ploton=0; #could turn plotting on for diagnostic purposes
plevel=0; #how much diagnostics to print to the screen. Should be 0 for real runs. Nothing shown for parallel fitting
max.stps = -1000; #max steps to take per optimizer iteration  - set to negative if we don't want to use this criterion. Used for both local and global solvers
max.time= 20 * 60*60 # runtime per scenario (in seconds, so 1*60*60 is 1 hour) - only one scenario/outer loop here6
max.stps.local = 10000; #max steps taken by local solver if we do the combined global/local (i.e. solver 14)
xerr.global=1e-8; xerr.local=1e-5; #used for multi-start solver MLS
xerr.general=1e-6; #used for all other solvers

real.time.start=date(); #get current time to measure length of optimization for each strain
print(sprintf('Optimization started at %s ',real.time.start))

parallel.comp=1; #turn on or off parallel computing -
node.num=length(solvertype); #number of sockets/nodes to use for parallel computing
node.type=1; #choose socket/node type. 1 for SOCK (can be run locally), 2 for MPI


if (parallel.comp==0) #standard non-parallel run using 1 core
{
    node.num=1; #just for printing purposes below
    node.type=1;
    reslist <- lapply(solvertype, outerfitfc)   #reslist contains best fits for all strains and all models
}
if (parallel.comp==1) #using snow package to do parallel computing
{
    if (node.type==1) {clust <- makeCluster(node.num, type = "SOCK")} #for local machines
    if (node.type==2) {clust <- makeCluster(node.num, type = "MPI")} #if run on a cluster
    clusterExport(clust,ls()) #make global variables available on each node/slave
    clusterEvalQ(clust, library(nloptr)) #load packages on each node
    clusterEvalQ(clust, library(deSolve))
    clusterEvalQ(clust, library(dplyr))
    reslist <- clusterApplyLB(clust, solvertype, outerfitfc)
    stopCluster(clust)
}


ct=1;
for (nn in 1:length(solvertype))
{
    
    ind=match(names(reslist[[ct]]),names(resmatrix[nn,]));
    resmatrix[nn,ind]=reslist[[ct]]
    ct=ct+1;
}
#save best fit parameter values and fitting diagnostics
save(resmatrix, file = filename.save);


real.time.stop=date();
tend=proc.time(); #capture current time
tdiff=tend-tstart;
runtime.minutes=tdiff[[3]]/60;

#print(resmatrix)
if (node.type==2) {mpi.quit() }

print(sprintf('Optimization ended at %s and took %f minutes using %d sockets with solver type %d',real.time.stop,tdiff[[3]]/60,node.num,solvertype));


#################################################
#do final runs to get time-series for plotting

txstatusvec=c('earlytx','latetx','notx')
txtimevec=c(29/24,50/24,500/24) #time of treatment start - last value means never

#take fit with smallest SSR (best fit) and process
minind = which.min(resmatrix[,"SSR"])
resvec = resmatrix[minind,];

#save time-series for each dose
timeseries.list=list()

timeseriesname = "haydenfit-timeseries-e1.Rdata"


for (nn in 1:nmax) #loop over treatment conditions
{
  
  txstatus = txstatusvec[nn]
  txtime = txtimevec[nn]
  
  t.int.vec=seq(0,10,by=0.01); #vector for integration times - units of hours
  Y0=c(U=U0, I=0, V=as.numeric(resvec["V1"]) ) #starting conditions
  
  odestack=ode(Y0,t.int.vec,odeequations,parms=resvec,atol=atolv,rtol=rtolv,method=odemethod,txtime=txtime); #runs the ODE equations
  timeseries.list[[nn]]=odestack
}

save(timeseries.list, file = timeseriesname);

#compute AICc for final fit
dp=nrow(mydata); #number of data points
xm=length(parnames)+1; #number of parameters - add extra if SSR framework is used 
AIC=2*xm+dp*log(as.numeric(resvec['SSR'])/dp); #for SSR approach
AICc=AIC+2*xm*(xm+1)/(dp-xm-1); #used if all parameters, including nuisance ones are fit. 
if ((dp-xm-1)<1) {AICc=NA;} #we have too many parameters, AICc makes no sense
print(sprintf('AICc is %f',AICc))

###################################################################
#end main program
###################################################################
