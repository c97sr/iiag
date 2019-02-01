## Function the check aphysicality of the variables/parameters
##  More information in 'Comparison of filtering methods for the modeling and retrospective forecasting of influenza epidemics' (PLoS Compute Biol)
##  by Wan Yang, Alicia Karspeck, and Jeffrey Shaman, 2014

Fn_checkxnobounds<-function(xnew){
    n.ens=dim(xnew)[2]; # get the number of particles
    n.var=dim(xnew)[1]; # get the number of variables

  ug=max(xnew[1,]);    #  Corrects if S>N
  if (ug>N){  
    for (jj in 1:n.ens){
      if (xnew[1,jj]>N){
        xnew[1,jj]=N-1;
      }
    }
  }
  ug=max(xnew[2,]);    #  Corrects if I>N
  if (ug>N){
    for (jj in 1:n.ens){
      if (xnew[2,jj]>N){
        xnew[2,jj]=median(xnew[2,]);
      }
    }
  }
  ug=min(xnew);     #  Corrects if any state or parameter nudges negative
  if (ug<=0){
      for (jj in 1:n.ens){
          for (ii in 1:n.var){
              if (xnew[ii,jj]<=0){
                  xnew[ii,jj]=max(mean(xnew[ii,]),1);
              }
          }
      }
  }
  ug=min(xnew[3,]);     #  Corrects if L < 200 days
  if (ug<200){
    for (jj in 1:n.ens){
      if (xnew[3,jj]<200){
        xnew[3,jj]=max(median(xnew[3,]),200);
      }
    }
  }
  ug=min(xnew[4,]);     #  Corrects if D < .5
  if (ug<=1){
    for (jj in 1:n.ens){
      if (xnew[4,jj]<=.5){
        xnew[4,jj]=max(median(xnew[4,]),.5);
      }
    }
  }
  ug=min(xnew[5,]-xnew[6,]);  #  Runs if R0mx <= R0mn
  if (ug<=0){
    for (jj in 1:n.ens){
      if (xnew[5,jj]<=xnew[6,jj]){
        xnew[5,jj]=xnew[6,jj]+0.01;
      }
    }
  }
  
  xnew;
}

Fn_checkDA<-function(xnew,bound.low,bound.up){
  b.low=bound.low;
  b.up=bound.up;
  n.var=nrow(xnew); n.ens=ncol(xnew);
  for(vi in 1:n.var){
    #  Corrects if <b.low
    ug=min(xnew[vi,]);
    if (ug<b.low[vi]){  
      for (jj in 1:n.ens){
        if (xnew[vi,jj]<b.low[vi]){
          xnew[vi,jj]=b.low[vi];
        }
      }
    }
    ug=max(xnew[vi,]);
    if (ug>b.up[vi]){  
      for (jj in 1:n.ens){
        if (xnew[vi,jj]>b.up[vi]){
          xnew[vi,jj]=b.up[vi];
        }
      }
    }
  }
  xnew;
}