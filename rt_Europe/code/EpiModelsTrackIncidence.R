# Script to run the humidity forcing SIRS models for influenza outbreak simulation
# Run continuously if by default
# The realdata switch controls whether to track incidence (newI), or just S and I
# Note: newI is cumulative incidence (not weekly incidence)
##  More information in 'Comparison of filtering methods for the modeling and retrospective forecasting of influenza epidemics' (PLoS Compute Biol)
##  by Wan Yang, Alicia Karspeck, and Jeffrey Shaman, 2014

## [july 27, 2015] fixed an error in numerical integration

propagateSIR<-function(tm_strt, tm_end, tm_step, S0, I0, N, D, L, beta, realdata=FALSE){
  # function to integrate to the next time step
  # use SIR model, integrate dicretely with Poisson distributions
  # input: tm_strt: starting time; tm_end: ending time; tm_step: time step
  #         S0, I0: initial states; N: population size
  #         D: infection period, day; L: immune period, day; 
  #         alpha: rate from exposed to infectious; beta: transmission matrix
  # output: S, I for all time steps
  cnt=1;
  tm_strt=tm_strt-tm.range[1]+1; # adjust the index to match beta
  tm_end=tm_end-tm.range[1]+1;
  tm_vec=seq(tm_strt,tm_end,by=tm_step)
  tm_sz=length(tm_vec)+1
  S=I=newI=rep(0,tm_sz)
  S[1]=S0; I[1]=I0; # R[1]=N-S0-I0;
  newI[1]=0;
  if(! exists("discrete")) discrete=FALSE; 
  if (discrete){
    S[1]=round(S0,0); I[1]=round(I0,0); # R[1]=N-S0-I0;
    for (t in tm_vec){
      cnt=cnt+1;
      
      Eimmloss=tm_step*(1/L*(N-S[cnt-1]-I[cnt-1]))
      Einf=tm_step*(beta[t,]%*%I[cnt-1]*S[cnt-1]/N)
      Erecov=tm_step*(1/D*I[cnt-1])
      Eimmloss[Eimmloss<0]=0   # adjust, set <0 to 0
      Einf[Einf<0]=0
      Erecov[Erecov<0]=0
      smcl=rpois(1,Eimmloss)
      smci=rpois(1,Einf)
      smcr=rpois(1,Erecov)
      
      sk1=smcl-smci
      ik1=smci-smcr
      ik1a=smci;
      Ts1=S[cnt-1]+round(sk1/2,0) # [july 27, 2015] fixed
      Ti1=I[cnt-1]+round(ik1/2,0)
      
      Eimmloss=tm_step*(1/L*(N-Ts1-Ti1))
      Einf=tm_step*(beta[t,]%*%Ti1*Ts1/N)
      Erecov=tm_step*(1/D*Ti1)
      Eimmloss[Eimmloss<0]=0   # adjust, set <0 to 0
      Einf[Einf<0]=0
      Erecov[Erecov<0]=0
      smcl=rpois(1,Eimmloss)
      smci=rpois(1,Einf)
      smcr=rpois(1,Erecov)
      sk2=smcl-smci
      ik2=smci-smcr
      ik2a=smci;
      Ts2=S[cnt-1]+round(sk2/2,0)
      Ti2=I[cnt-1]+round(ik2/2,0)
      
      Eimmloss=tm_step*(1/L*(N-Ts2-Ti2))
      Einf=tm_step*(beta[t,]%*%Ti2*Ts2/N)
      Erecov=tm_step*(1/D*Ti2)
      Eimmloss[Eimmloss<0]=0   # adjust, set <0 to 0
      Einf[Einf<0]=0
      Erecov[Erecov<0]=0                   
      smcl=rpois(1,Eimmloss)
      smci=rpois(1,Einf)
      smcr=rpois(1,Erecov)
      sk3=smcl-smci
      ik3=smci-smcr
      ik3a=smci;
      Ts3=S[cnt-1]+round(sk3,0)
      Ti3=I[cnt-1]+round(ik3,0)
      
      Eimmloss=tm_step*(1/L*(N-Ts3-Ti3))
      Einf=tm_step*(beta[t,]%*%Ti3*Ts3/N)
      Erecov=tm_step*(1/D*Ti3)
      Eimmloss[Eimmloss<0]=0   # adjust, set <0 to 0
      Einf[Einf<0]=0
      Erecov[Erecov<0]=0
      smcl=rpois(1,Eimmloss)
      smci=rpois(1,Einf)
      smcr=rpois(1,Erecov)
      sk4=smcl-smci
      ik4=smci-smcr
      ik4a=smci;
      
      seed=rpois(1,.1)
      S[cnt]=S[cnt-1]+round(sk1/6+sk2/3+sk3/3+sk4/6,0)-seed
      I[cnt]=I[cnt-1]+round(ik1/6+ik2/3+ik3/3+ik4/6,0)+seed
      newI[cnt]=round(newI[cnt-1]+ik1a/6+ik2a/3+ik3a/3+ik4a/6+seed,0);  
      # R[cnt]=N-S[cnt]-I[cnt]
    }
  } else {
    # run continuously
    for (t in tm_vec){
      cnt=cnt+1;
      
      Eimmloss=tm_step*(1/L*(N-S[cnt-1]-I[cnt-1]))
      Einf=tm_step*(beta[t,]%*%I[cnt-1]*S[cnt-1]/N)
      Erecov=tm_step*(1/D*I[cnt-1])
      Eimmloss[Eimmloss<0]=0   # adjust, set <0 to 0
      Einf[Einf<0]=0
      Erecov[Erecov<0]=0
      smcl=Eimmloss
      smci=Einf
      smcr=Erecov
      
      sk1=smcl-smci
      ik1=smci-smcr
      ik1a=smci;
      Ts1=S[cnt-1]+round(sk1/2,0)
      Ti1=I[cnt-1]+round(ik1/2,0)
      
      Eimmloss=tm_step*(1/L*(N-Ts1-Ti1))
      Einf=tm_step*(beta[t,]%*%Ti1*Ts1/N)
      Erecov=tm_step*(1/D*Ti1)
      Eimmloss[Eimmloss<0]=0   # adjust, set <0 to 0
      Einf[Einf<0]=0
      Erecov[Erecov<0]=0
      smcl=Eimmloss
      smci=Einf
      smcr=Erecov
      sk2=smcl-smci
      ik2=smci-smcr
      ik2a=smci;
      Ts2=S[cnt-1]+round(sk2/2,0)
      Ti2=I[cnt-1]+round(ik2/2,0)
      
      Eimmloss=tm_step*(1/L*(N-Ts2-Ti2))
      Einf=tm_step*(beta[t,]%*%Ti2*Ts2/N)
      Erecov=tm_step*(1/D*Ti2)
      Eimmloss[Eimmloss<0]=0   # adjust, set <0 to 0
      Einf[Einf<0]=0
      Erecov[Erecov<0]=0                   
      smcl=Eimmloss
      smci=Einf
      smcr=Erecov
      sk3=smcl-smci
      ik3=smci-smcr
      ik3a=smci;
      Ts3=S[cnt-1]+sk3
      Ti3=I[cnt-1]+ik3
      
      Eimmloss=tm_step*(1/L*(N-Ts3-Ti3))
      Einf=tm_step*(beta[t,]%*%Ti3*Ts3/N)
      Erecov=tm_step*(1/D*Ti3)
      Eimmloss[Eimmloss<0]=0   # adjust, set <0 to 0
      Einf[Einf<0]=0
      Erecov[Erecov<0]=0
      smcl=Eimmloss
      smci=Einf
      smcr=Erecov
      sk4=smcl-smci
      ik4=smci-smcr
      ik4a=smci;
      
      # seed=rpois(1,.1)
      seed=.1;
      S[cnt]=S[cnt-1]+round(sk1/6+sk2/3+sk3/3+sk4/6,0)-seed
      I[cnt]=I[cnt-1]+round(ik1/6+ik2/3+ik3/3+ik4/6,0)+seed
      newI[cnt]=newI[cnt-1]+ik1a/6+ik2a/3+ik3a/3+ik4a/6+seed;
      # R[cnt]=N-S[cnt]-I[cnt]
    }
  }
  
  # S=round(S,0); I=round(I,0);
  if (realdata==FALSE){
    rec=cbind(S,I,deparse.level=0); 
  } else {
    rec=cbind(S,I,newI,deparse.level=0); 
  }
  
  rec;
}



# function to integrate forward in time with SIRS model with vector opporation for multiple simulations. 
propagateParSIR<-function(tm_strt, tm_end, tm_step, S0, I0, N, D, L, beta, realdata=FALSE){
  # function to integrate to the next time step
  # use SIR model, integrate dicretely with Poisson distributions
  # input: tm_strt: starting time; tm_end: ending time; tm_step: time step
  #         S0, I0: initial states; N: population size
  #         D: infection period, day; L: immune period, day; 
  #         alpha: rate from exposed to infectious; beta: transmission matrix
  # output: S, I for all time steps
  cnt=1;
  # beta stores only data during the time used for the truth
  tm_strt=tm_strt-tm.range[1]+1; # adjust the index to match beta
  tm_end=tm_end-tm.range[1]+1;
  tm_vec=seq(tm_strt,tm_end,by=tm_step)
  tm_sz=length(tm_vec)+1; # including the initial conditions and results of length(tm_vec) integration steps
  Np=length(S0); # number of particles
  # if integrating parallelly, S0, and I0 should be passed in as a vector (Np particles)
  # but if so, can't do multiple age groups
  S=I=newI=matrix(0,Np,tm_sz)
  S[,1]=S0; I[,1]=I0; # R[,1]=N-S0-I0;
  newI[,1]=0;
  if(! exists("discrete")) discrete=FALSE; 
  
  if (discrete){
    S[,1]=round(S0,0); I[,1]=round(I0,0); # R[,1]=N-S0-I0;
    for (t in tm_vec){
      cnt=cnt+1;
      
      Eimmloss=tm_step*(1/L*(N-S[,cnt-1]-I[,cnt-1]))
      Einf=tm_step*(beta[t,]*I[,cnt-1]*S[,cnt-1]/N)
      Erecov=tm_step*(1/D*I[,cnt-1])
      Eimmloss[Eimmloss<0]=0   # adjust, set <0 to 0
      Einf[Einf<0]=0
      Erecov[Erecov<0]=0
      smcl=rpois(Np,Eimmloss)
      smci=rpois(Np,Einf)
      smcr=rpois(Np,Erecov)
      
      sk1=smcl-smci
      ik1=smci-smcr
      ik1a=smci;
      Ts1=S[,cnt-1]+round(sk1/2,0) # [july 27, 2015] fixed
      Ti1=I[,cnt-1]+round(ik1/2,0)
      
      Eimmloss=tm_step*(1/L*(N-Ts1-Ti1))
      Einf=tm_step*(beta[t,]*Ti1*Ts1/N)
      Erecov=tm_step*(1/D*Ti1)
      Eimmloss[Eimmloss<0]=0   # adjust, set <0 to 0
      Einf[Einf<0]=0
      Erecov[Erecov<0]=0
      smcl=rpois(Np,Eimmloss)
      smci=rpois(Np,Einf)
      smcr=rpois(Np,Erecov)
      sk2=smcl-smci
      ik2=smci-smcr
      ik2a=smci;
      Ts2=S[,cnt-1]+round(sk2/2,0)
      Ti2=I[,cnt-1]+round(ik2/2,0)
      
      Eimmloss=tm_step*(1/L*(N-Ts2-Ti2))
      Einf=tm_step*(beta[t,]*Ti2*Ts2/N)
      Erecov=tm_step*(1/D*Ti2)
      Eimmloss[Eimmloss<0]=0   # adjust, set <0 to 0
      Einf[Einf<0]=0
      Erecov[Erecov<0]=0                   
      smcl=rpois(Np,Eimmloss)
      smci=rpois(Np,Einf)
      smcr=rpois(Np,Erecov)
      sk3=smcl-smci
      ik3=smci-smcr
      ik3a=smci;
      Ts3=S[,cnt-1]+round(sk3,0)
      Ti3=I[,cnt-1]+round(ik3,0)
      
      Eimmloss=tm_step*(1/L*(N-Ts3-Ti3))
      Einf=tm_step*(beta[t,]*Ti3*Ts3/N)
      Erecov=tm_step*(1/D*Ti3)
      Eimmloss[Eimmloss<0]=0   # adjust, set <0 to 0
      Einf[Einf<0]=0
      Erecov[Erecov<0]=0
      smcl=rpois(Np,Eimmloss)
      smci=rpois(Np,Einf)
      smcr=rpois(Np,Erecov)
      sk4=smcl-smci
      ik4=smci-smcr
      ik4a=smci;
      
      seed=rpois(Np,.1)
      # seed=.1;
      S[,cnt]=S[,cnt-1]+round(sk1/6+sk2/3+sk3/3+sk4/6,0)-seed
      I[,cnt]=I[,cnt-1]+round(ik1/6+ik2/3+ik3/3+ik4/6,0)+seed
      newI[,cnt]=round(newI[,cnt-1]+ik1a/6+ik2a/3+ik3a/3+ik4a/6+seed,0);
      # R[,cnt]=N-S[,cnt]-I[,cnt]
    }
  } else {
    # run continuously
    for (t in tm_vec){
      cnt=cnt+1;
      
      Eimmloss=tm_step*(1/L*(N-S[,cnt-1]-I[,cnt-1]))
      Einf=tm_step*(beta[t,]*I[,cnt-1]*S[,cnt-1]/N)
      Erecov=tm_step*(1/D*I[,cnt-1])
      Eimmloss[Eimmloss<0]=0   # adjust, set <0 to 0
      Einf[Einf<0]=0
      Erecov[Erecov<0]=0
      smcl=Eimmloss
      smci=Einf
      smcr=Erecov
      
      sk1=smcl-smci
      ik1=smci-smcr
      ik1a=smci;
      Ts1=S[,cnt-1]+sk1/2
      Ti1=I[,cnt-1]+ik1/2
      
      Eimmloss=tm_step*(1/L*(N-Ts1-Ti1))
      Einf=tm_step*(beta[t,]*Ti1*Ts1/N)
      Erecov=tm_step*(1/D*Ti1)
      Eimmloss[Eimmloss<0]=0   # adjust, set <0 to 0
      Einf[Einf<0]=0
      Erecov[Erecov<0]=0
      smcl=Eimmloss
      smci=Einf
      smcr=Erecov
      sk2=smcl-smci
      ik2=smci-smcr
      ik2a=smci;
      Ts2=S[,cnt-1]+sk2/2
      Ti2=I[,cnt-1]+ik2/2
      
      Eimmloss=tm_step*(1/L*(N-Ts2-Ti2))
      Einf=tm_step*(beta[t,]*Ti2*Ts2/N)
      Erecov=tm_step*(1/D*Ti2)
      Eimmloss[Eimmloss<0]=0   # adjust, set <0 to 0
      Einf[Einf<0]=0
      Erecov[Erecov<0]=0                   
      smcl=Eimmloss
      smci=Einf
      smcr=Erecov
      sk3=smcl-smci
      ik3=smci-smcr
      ik3a=smci;
      Ts3=S[,cnt-1]+sk3
      Ti3=I[,cnt-1]+ik3
      
      Eimmloss=tm_step*(1/L*(N-Ts3-Ti3))
      Einf=tm_step*(beta[t,]*Ti3*Ts3/N)
      Erecov=tm_step*(1/D*Ti3)
      Eimmloss[Eimmloss<0]=0   # adjust, set <0 to 0
      Einf[Einf<0]=0
      Erecov[Erecov<0]=0
      smcl=Eimmloss
      smci=Einf
      smcr=Erecov
      sk4=smcl-smci
      ik4=smci-smcr
      ik4a=smci;
      
      # seed=rpois(Np,.1)
      seed=.1;
      S[,cnt]=S[,cnt-1]+sk1/6+sk2/3+sk3/3+sk4/6-seed
      I[,cnt]=I[,cnt-1]+ik1/6+ik2/3+ik3/3+ik4/6+seed
      newI[,cnt]=newI[,cnt-1]+ik1a/6+ik2a/3+ik3a/3+ik4a/6+seed;
      # R[,cnt]=N-S[,cnt]-I[,cnt]
    }
  }
  
  # S=t(round(S,0)); I=t(round(I,0));
  S=t(S); I=t(I); newI=t(newI);
  if (realdata==FALSE){
    rec=list(S=S,I=I); 
  } else {
    rec=list(S=S,I=I,newI=newI); 
  }
  rec;
}


