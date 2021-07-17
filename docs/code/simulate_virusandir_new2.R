# Simulation of a viral infection model with an immune response
# Modified version of simulate_virusandir_ode from DSAIRM
# This is the second modified function we built
# This function builds on simulate_virusandir_new and adds two more equations
# namely for a second type of B cells and antibodies (B2 and A2)
# the idea is that B2 and A2 are immune memory from a previous infection,
# which means A2 are the antibodies that will be detrimental and could cause ADE (through the modified term b2*A2*V*U) 
# B and A are now the antibodies and B-cells specific for the new infection,
# which means those do not cause ADE
# Rscript6.R calls and uses this function

simulate_virusandir_new2 <- function(U = 1e5, I = 0, V = 10, T = 0, B = 0, A = 0, B2 = 0, A2 = 0, D=0, 
                                    n=0, dU = 0, dI = 1, dV = 4, 
                                    b = 1e-5, p = 1e3, 
                                    sF=1e-2,kA=1e-5,kT=1e-5,pF=1,dF=1,gF=1,Fmax=1e3,
                                    hV=1e-6,hF=1e-5,gB=1,gT=1e-4,rT=0.5,rA=10,dA=0.2, 
                                    gB2=1,rA2=10,dA2=0.2,kA2=1e-5,
                                    b2 = 0,
                                    tstart = 0, tfinal = 30, dt = 0.05)
{

  #function that specificies the ode model
  virusandirode <- function(t, y, parms)
  {
    with(
      as.list(c(y,parms)), #lets us access variables and parameters stored in y and parms by name
      {

        dUdt = n - dU*U - b*V*U - b2*A2*V*U #extra term for ADE
        dIdt = b*V*U + b2*A2*V*U - dI*I - kT*T*I #extra term for ADE
        dVdt = p*I/(1+sF*F) - dV*V - b*V*U - kA*A*V - b2*A2*V*U - kA2*A2*V #extra term for ADE
        dFdt = pF - dF*F + V / (V + hV) * gF * (Fmax - F)
        dTdt = F * V * gT + rT * T
        dBdt = F * V / (F * V + hF) * gB * B
        dAdt = rA*B  - dA*A - kA*A*V

        dB2dt = F * V / (F * V + hF) * gB2 * B2 #previous B cells
        dA2dt = rA2*B2  - dA2*A2 - kA2*A2*V #previous Ab
        
        dDdt = dI*I + kT*T*I #keeping track of all dead cells
        
        list(c(dUdt, dIdt, dVdt, dFdt, dTdt, dBdt, dAdt, dB2dt, dA2dt, dDdt))
      }
    ) #close with statement
  } #end function specifying the ODEs


  #combine initial conditions into a vector
  #some initial conditions are set to fixed values and can't be adjusted in the app
  Y0 = c(U = U, I = I, V = V, F=pF/dF, T=T, B=B, A=A, B2=B2, A2=A2, D=D);
  timevec = seq(tstart, tfinal, by = dt); #vector of times for which solution is returned (not that internal timestep of the integrator is different)

  #combining parameters into a parameter vector
  pars = c(n=n,dU=dU,dI=dI,dV=dV,b=b,p=p,sF=sF,kA=kA,kT=kT,pF=pF,dF=dF,gF=gF,Fmax=Fmax,hV=hV,hF=hF,gB=gB,gT=gT,rT=rT,rA=rA,dA=dA, gB2=gB2,rA2=rA2,dA2=dA2,kA2=kA2,b2 = b2);

  #this line runs the simulation, i.e. integrates the differential equations describing the infection process
  #the result is saved in the odeoutput matrix, with the 1st column the time, all other column the model variables
  #in the order they are passed into Y0 (which needs to agree with the order in virusode)
  odeoutput = deSolve::ode(y = Y0, times = timevec, func = virusandirode, parms=pars, atol=1e-12, rtol=1e-12);

  #return result as list, with element ts containing the time-series
  result = list()
  result$ts = as.data.frame(odeoutput)
  return(result)
}
