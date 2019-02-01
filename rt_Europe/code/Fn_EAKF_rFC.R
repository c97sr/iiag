## R Function to run Ensemble Adjustment Kalman Filter with inflation and then forecast
## Based on Anderson JL (2001) Mon Weather Rev 129: 2884-2903.
##          and Karspeck AR, Anderson JL (2007) J Climate 20: 4638-4658.
##          and Shaman J, Karspeck A (2012) Proc Natl Acad Sci USA 109: 20425-20430.
## nsn: number of weeks in a season, b/c it's run retrospectively, we have the whole time series
## ntrn: number of weeks for the training process, beyond which we do forecast
##  More information in 'Comparison of filtering methods for the modeling and retrospective forecasting of influenza epidemics' (PLoS Compute Biol)
##  by Wan Yang, Alicia Karspeck, and Jeffrey Shaman, 2014

EAKF_rFC<-function(num_ens, tmstep, param.bound, obs_i=obs_i,ntrn=1, 
                    obs_vars,tm.ini=273, tm.range=273:500){
  ## EAKF to run retrospecitve forecast of influenza
  ## Inputs: 
  ##        num_ens: number of ensemble members
  ##        tmstep: time step (e.g., 7 days if weekly observations are used)
  ##        param.bound: bounds for the parameters 
  ##          (a 2-column matrix prescribed the lower (col#1) and upper (col#2 of the bounds))
  ##        obs_i: observations (the time series);
  ##        ntrn: number of observations for training the filter
  ##        obs_vars: variances of the obseravations
  ##        tm.ini: the initial time (the number of day in the year)
  ##        tm.range: the sequence of day from the first to the last day of simulation 
  ##         from the current year to the next year of a season (used to match the time of climatological humidity)
  ## Outputs:
  ##        train: modeled time series of the variables/parameters (posteriors) for the training period
  ##        trainsd: standard deviation of the variables/parameters during the training
  ##        fcast: forecasted time series of the variables (i.e., S, I, and newI)
  ##        fcsd: standard deviation of the variables of the forecast
  ##        metrics: root mean squared error (rms), correlation with the observations, etc.
  ##        trainprior: prior estimates of the variables/parameters during the training
  ##        Note: if only wish to save the metrics, set metricsonly=TRUE
  
  library("truncnorm"); library("tgp"); # for lhs
  library("MASS"); # for multivariate normal distribution
  source(paste(dir.code,"EpiModelsTrackIncidence.R",sep="")); # SIRS model
  source(paste(dir.code,"Fn_checkxnobounds.R",sep="")); # function to check data aphysicality
  ## read list of variables (generated from a latin hypercube for the prior of variables/parameter)
  if (! exists("phi")){
    phi=read.csv(paste(dir.data,'phi.csv',sep=""),header=F);
    params=phi[,1:4];
    susceps=infects=NULL; 
    for (i in 5:31){
      susceps=append(susceps,phi[,i])
      infects=append(infects,phi[,27+i])
    }
  }

  num_times=floor(length(tm.range)/tmstep);
  nsn=nsn # length(obs_i); 
  nfc=nsn-ntrn; # number of weeks for forecasting;
  
  theta_low=param.bound[,1];
  theta_up=param.bound[,2];
  
  So=matrix(0,6,num_ens);
  xprior=array(0,c(7,num_ens,ntrn+1));
  xpost=array(0,c(7,num_ens,ntrn));
  fcast=array(0,c(3,num_ens,nfc)); # fcast: S, I, newI;
  
  rnd=cbind(ceiling(27*1e4*runif(num_ens)), ceiling(27*1e4*runif(num_ens)), 
            ceiling(1e4*runif(num_ens)));
  So[1,]=susceps[rnd[,1]]/5; # per 100,000 population
  So[2,]=infects[rnd[,2]]/5;
  So[3,]=params[rnd[,3],1]*365;
  So[4,]=params[rnd[,3],2]*365;
  So[5,]=params[rnd[,3],3];
  So[6,]=params[rnd[,3],4];

  ## Calculate the reproductive number at time t BT1 and the transmission rate 
  ## according to Shaman J, Karspeck A (2012) Proc Natl Acad Sci USA 109: 20425-20430.
  beta.range=tm.range[1]:(tail(tm.range,1)+2*tmstep); 
  ## make sure to include enough beta, b/c the prior integrate 1 time step more than needed
  AHpt=AH[beta.range,]; # only use AH during the time of interest, to minimize space needed
  ones1=matrix(1,length(AHpt),1); # matrix of 1 to make the dimension match
  b=log(So[5,]-So[6,]); b=ones1%*%b; # expand b to a matrix with each col for each ensemble member
  a=-180;
  ones2=matrix(1,1,num_ens); AHpt=as.matrix(AHpt,length(AHpt),1)
  AHpt=AHpt%*%ones2
  BT1=exp(a*AHpt+b)+ones1%*%So[6,];
  beta=BT1/(ones1%*%So[4,]);
  tcurrent=tm.ini;
  ## integrate forward 1 step foreward
  Sr_tmp=propagateParSIR(tcurrent+dt,tcurrent+tmstep,dt,S0=So[1,], 
                         I0=So[2,], N, D=So[4,], L=So[3,], beta, realdata=T)
  xprior[1,,1]=tail(Sr_tmp$S,1) 
  xprior[2,,1]=tail(Sr_tmp$I,1)
  xprior[3:6,,1]=So[3:6,];
  xprior[7,,1]=tail(Sr_tmp$newI,1)
  
  ####  Define the mapping operator H
  H=7;   #  observing the 7th variable (new infections)
  
  #### Begin looping through observations
  #### Training process
  for (tt in 1:ntrn){
    if(!is.na(obs_i[tt])) { # Update state variables and parameters, then integrate forward
      # inflate all states and parameters
      inflat=diag(x=c(lambda,lambda,lambda,lambda,lambda,lambda,lambda),7,7);
      xmn=rowMeans(xprior[,,tt]);
      xprior[,,tt]=inflat%*%(xprior[,,tt]-xmn%*%matrix(1,1,num_ens))+xmn%*%matrix(1,1,num_ens)
      
      ####  Get the variance of the ensemble
      obs_var = obs_vars[tt];
      prior_var = var(xprior[H,,tt]);
      post_var = prior_var*obs_var/(prior_var+obs_var);
      
      prior_mean = mean(xprior[H,,tt]);
      post_mean = post_var*(prior_mean/prior_var + obs_i[tt]/obs_var);
      
      #### Compute alpha and adjust distribution to conform to posterior moments
      alp = sqrt(obs_var/(obs_var+prior_var));
      
      dy = post_mean + alp*((xprior[H,,tt])-prior_mean)-xprior[H,,tt];
      ###  Getting the covariance of the prior state space and
      ###  observations  (which could be part of state space, e.g. infections)
      rr=NULL;
      for (j in 1:dim(xprior[,,tt])[1]){
        C=cov(xprior[j,,tt],xprior[H,,tt])/prior_var;  # covariance/variance of x.obs
        rr=append(rr,C);
      }
      dx=rr%*%t(dy);
      
      ###  Get the new ensemble and save prior and posterior
      xnew = xprior[,,tt] + dx;
      
      #  Corrections to data aphysicalities
      xnew[1:6,]=Fn_checkxnobounds(xnew[1:6,]);
      
      xpost[,,tt]=xnew;
      
      #  Integrate forward one time step
      b=log(xpost[5,,tt]-xpost[6,,tt]); b=ones1%*%b; # expand b to a matrix with each col for each particle
      a=-180;
      BT1=exp(a*AHpt+b)+ones1%*%xpost[6,,tt];
      beta=BT1/(ones1%*%xpost[4,,tt]); 
      
      tcurrent = tm.ini+tmstep*tt;
      Sr_tmp=propagateParSIR(tcurrent+dt,tcurrent+tmstep,dt,xpost[1,,tt],xpost[2,,tt],N,
                             D=xpost[4,,tt],L=xpost[3,,tt],beta,realdata=T)
      xprior[1,,tt+1]=tail(Sr_tmp$S,1);
      xprior[2,,tt+1]=tail(Sr_tmp$I,1);
      xprior[3:6,,tt+1]=xpost[3:6,,tt];
      xprior[7,,tt+1]=tail(Sr_tmp$newI,1);
      
    } else { ## Simply integrate forward without updating anything
      ##  Don't use xpost for propagating - use xprior that resulted from previous integration
      ##  Integrate forward one time step
      b=log(xprior[5,,tt]-xprior[6,,tt]); b=ones1%*%b; # expand b to a matrix with each col for each particle
      a=-180;
      BT1=exp(a*AHpt+b)+ones1%*%xprior[6,,tt];
      beta=BT1/(ones1%*%xprior[4,,tt]); 
      tcurrent = tm.ini+tmstep*tt;
      Sr_tmp=propagateParSIR(tcurrent+dt,tcurrent+tmstep,dt,xprior[1,,tt],xprior[2,,tt],N,
                             D=xprior[4,,tt],L=xprior[3,,tt],beta,realdata=T)
      xprior[1,,tt+1]=tail(Sr_tmp$S,1);
      xprior[2,,tt+1]=tail(Sr_tmp$I,1);
      xprior[3:6,,tt+1]=xprior[3:6,,tt];
      xprior[7,,tt+1]=tail(Sr_tmp$newI,1);
    }
    
  } # end of training
  
  #### Forecast
  b=log(xpost[5,,ntrn]-xpost[6,,ntrn]); b=ones1%*%b; # expand b to a matrix with each col for each particle
  a=-180;
  BT1=exp(a*AHpt+b)+ones1%*%xpost[6,,ntrn];
  beta=BT1/(ones1%*%xpost[4,,ntrn]); 
  tcurrent = tm.ini+tmstep*ntrn;
  Sr_tmp=propagateParSIR(tcurrent+dt,tcurrent+tmstep*nfc,dt,xpost[1,,ntrn],xpost[2,,ntrn],N,
                         D=xpost[4,,ntrn],L=xpost[3,,ntrn],beta,realdata=T)
  fcast[1,,]=t(Sr_tmp$S[tmstep*(1:nfc)+1,]); # num_ens, time 
  fcast[2,,]=t(Sr_tmp$I[tmstep*(1:nfc)+1,]);
  fcast[3,,]=t(Sr_tmp$newI[tmstep*(1:nfc)+1,]-Sr_tmp$newI[tmstep*(0:(nfc-1))+1,]);
  ## get weekly incidence. (Note: newI in the SIRS function is cummulative new cases)
  ## calculate the mean of ensemble
  xprior_mean=xpost_mean=xsd=matrix(0,7,ntrn)
  for (tt in 1:ntrn){
    xprior_mean[,tt]=apply(xprior[,,tt],1,mean, na.rm=T)
    xpost_mean[,tt]=apply(xpost[,,tt],1,mean, na.rm=T)
    xsd[,tt]=apply(xpost[,,tt],1,sd, na.rm=T)
  }
  
  obs.na <- which(is.na(obs_i))[which(is.na(obs_i)) %in% 1:(dim(xpost_mean)[2])]
  xpost_mean[,obs.na] <- NA
  
  fcast_mean=fcast_sd=matrix(0,3,nfc)
  for (tt in 1:nfc){
    fcast_mean[,tt]=apply(fcast[,,tt],1,mean, na.rm=T);
    fcast_sd[,tt]=apply(fcast[,,tt],1,sd, na.rm=T);
  }
  
  #### Metrics for comparison
  Y=c(xpost_mean[7,],fcast_mean[3,]);  # newI
  pkwks=rep(0,num_ens);
  for (i in 1:num_ens){
    pkwks[i]=which.max(c(xpost[7,i,1:ntrn],fcast[3,i,1:nfc]));
  }
  pkwk_mode=MODE(pkwks)[1];
  pkwk_mode_perc=MODE(pkwks)[2]/num_ens;
  leadpkwk_mode=pkwk_mode-ntrn;
  pkwk_var=var(pkwks);
  
  corr=cor(Y[1:length(obs_i)], obs_i, use='pairwise.complete.obs') # correlation for fitting only
  rms=sqrt(mean((Y[1:length(obs_i)] - obs_i)^2,na.rm=T)) # RMSE for fitting only
  
  ## check the mean as well
  pkwk_mean=which.max(Y);
  leadpkwk_mean=pkwk_mean-ntrn;
  
  ## intensity
  peak_intensity <- max(Y,na.rm=T)
  ili_peaks = rep(0,num_ens);
  for (i in 1:num_ens){
    ili_peaks[i] = max(c(xpost[7,i,1:ntrn], fcast[3,i,1:nfc]));
  }
  peak_intensity_var <- var(ili_peaks)
  
  ## AR
  tot_attack <- sum(Y[!is.na(obs_i)])
  
  ## next 4 weeks
  ## these are manually added to the metrics collection at the end of this file
  
  ## onset, end, duration
  #onset5 <- findOnset(Y, 500)$onset
  end5 <- findOnset(Y, 500)$end
  dur5 <- findOnset(Y, 500)$duration
  
  ## Dists
  onsets5 = ends = durations = rep(0,num_ens)
  peakWeeks = peakIntensities = rep(NA, num_ens)
  nextILI = matrix(NA, nrow=4, ncol=num_ens)
  for (i in 1:num_ens){
    yy=c(xpost[7,i,1:ntrn], fcast[3,i,1:nfc]); # trajectory for that particle
    
    onsets5[i] = findOnset(yy, 500)$onset
    ends[i] = findOnset(yy, 500)$end
    durations[i] = findOnset(yy, 500)$duration
    
    ## find peak week as predicted by this ensemble member
    peakWeeks[i]  = which.max(yy);
    if(!is.na(peakWeeks[i]))
      peakWeeks[i] = peakWeeks[i] + wk_start - 1
    
    ## find peak intensity as predicted by this ensemble member
    peakIntensities[i] = max(yy, na.rm=T)
    
    ## find forecasts of this ensemble member for the next 4 weeks; row_i = ith week ahead
    for(j in 1:min(nfc, 4)){
      nextILI[j, i] = fcast[3, i, j]
    }
    
  } ## end of ensemble loop
  
  ## calculate prob. distribution for onsets
  onsets5DistNA = matrix(c(-1, round(length(onsets5[is.na(onsets5)])/length(onsets5), 4)), nrow=1, ncol=2)
  onsets5 = onsets5[!is.na(onsets5)];
  onsets5Dist = matrix(NA, nrow=length(unique(onsets5)), ncol=2)
  row=1
  for(i in sort(unique(onsets5))){
    onsets5Dist[row, 1] = i;
    onsets5Dist[row, 2] = round(length(onsets5[onsets5==i])/300, 4)
    row=row+1
  }
  if (onsets5DistNA[, 2] > max(onsets5Dist[, 2])) {
    onset5 <- NA
  } else {
    onset5 <- findOnset(Y, 500)$onset
  }
  onsets5Dist = rbind(onsets5Dist, onsets5DistNA)
  
  endsDistNA = matrix(c(-1, round(length(ends[is.na(ends)])/length(ends), 4)), nrow=1, ncol=2)
  ends = ends[!is.na(ends)]
  endsDist = matrix(NA, nrow=length(unique(ends)), ncol=2)
  row=1
  for(i in sort(unique(ends))) {
    endsDist[row, 1] = i
    endsDist[row, 2] = round(length(ends[ends==i])/length(ends), 4)
    row=row+1
  }
  endsDist = rbind(endsDist, endsDistNA)
  
  dursDistNA = matrix(c(-1, round(length(durations[is.na(durations)])/length(durations), 4)), nrow=1, ncol=2)
  durations = durations[!is.na(durations)]
  dursDist = matrix(NA, nrow=length(unique(durations)), ncol=2)
  row=1
  for(i in sort(unique(durations))) {
    dursDist[row, 1] = i
    dursDist[row, 2] = round(length(durations[durations==i])/length(durations), 4)
    row=row+1
  }
  dursDist = rbind(dursDist, dursDistNA)
  
  ## calculate prob. distribution for peak weeks
  peakWeeksDistNA = matrix(c(-1, round(length(peakWeeks[is.na(peakWeeks)])/length(peakWeeks), 4)), nrow=1, ncol=2)
  
  peakWeeks = peakWeeks[!is.na(peakWeeks)]
  peakWeeksDist = matrix(NA, nrow=length(unique(peakWeeks)), ncol=2)
  row=1
  for(i in sort(unique(peakWeeks))){
    peakWeeksDist[row, 1] = i;
    peakWeeksDist[row, 2] = round(length(peakWeeks[peakWeeks==i])/length(peakWeeks), 4)
    row=row+1
  } 
  peakWeeksDist = rbind(peakWeeksDist, peakWeeksDistNA)
  
  ## calculate prob. distribution for peak intensities
  reqLimits = seq(0, 1e4, by=1e3)
  peakIntensities=peakIntensities[!is.na(peakIntensities)]
  peakIntensitiesDist = matrix(NA, nrow=length(reqLimits), ncol=2)
  row=1
  for(i in 2:length(reqLimits)){
    peakIntensitiesDist[row, 1] = reqLimits[i]
    peakIntensitiesDist[row, 2] = round(length(peakIntensities[peakIntensities >= reqLimits[i-1] & peakIntensities < reqLimits[i]])/length(peakIntensities), 4)
    row=row+1    
  }
  ## > 1e4
  peakIntensitiesDist[row, 1] = 1e5
  peakIntensitiesDist[row, 2] = round(length(peakIntensities[peakIntensities >= max(reqLimits)])/length(peakIntensities), 4)
  
  ## calculate prob. distribution for next 4 week forecasts
  ## we bin in increments of 1e3 upto 1e4. An extra bin for >1e4
  nextILIDist = matrix(NA, nrow=length(reqLimits), ncol=nrow(nextILI)+1)
  row=1
  for(i in 2:length(reqLimits)){
    nextILIDist[row, 1] = reqLimits[i]
    for(j in 1:nrow(nextILI)){
      values = nextILI[j,]
      values = values[!is.na(values)]
      nextILIDist[row, j+1] = round(length(values[values >= reqLimits[i-1] & values < reqLimits[i]]) / length(values), 4)
    }
    row=row+1
  }
  
  ## > 1e4
  nextILIDist[row, 1] = 1e5
  for(j in 1:nrow(nextILI)){
    values = nextILI[j,]
    values = values[!is.na(values)]
    nextILIDist[row, j+1] = round(length(values[values  >= max(reqLimits)]) / length(values), 4)
  }
  
  ## Calculate onset variances
  onset5_var=var(onsets5,na.rm=T)
  end_var=var(ends,na.rm=T)
  dur_var=var(durations,na.rm=T)
  
  ## output data
  ## output the prediction from the last iteration
  tstep=seq(tm.ini+tmstep,nsn*tmstep+tm.ini,by=tmstep); # time
  fc_start=ntrn+wk_start-1; # the time of forecast
  
  out1=cbind(rep(fc_start,length((1:ntrn)[!is.na(obs_i[1:ntrn])])),
             tstep[(1:ntrn)[!is.na(obs_i[1:ntrn])]],
             (1:ntrn)[!is.na(obs_i[1:ntrn])]+39,
             t(xpost_mean)[(1:ntrn)[!is.na(obs_i[1:ntrn])],],
             t(xsd)[(1:ntrn)[!is.na(obs_i[1:ntrn])],]);
  # out2=cbind(rep(fc_start,ntrn),tstep[1:ntrn],t(xsd));
  out3=cbind(rep(fc_start,nfc),tstep[(ntrn+1):nsn],((ntrn+1):nsn)+39,t(fcast_mean),t(fcast_sd));
  # out4=cbind(rep(fc_start,nfc),tstep[(ntrn+1):nsn],t(fcast_sd));
  out5=t(c(fc_start,pkwk_mode+wk_start-1,pkwk_mean+wk_start-1,leadpkwk_mode,leadpkwk_mean,sqrt(pkwk_var),
           peak_intensity,sqrt(peak_intensity_var),tot_attack,fcast_mean[3,1],fcast_mean[3,2],fcast_mean[3,3],
           fcast_mean[3,4],rms,corr,onset5,sqrt(onset5_var),end5,sqrt(end_var),dur5,sqrt(dur_var)))
#   out5=t(c(fc_start,rms,
#            corr,delta_sum_newI,delta_pkwk,leadpkwk_mode,pkwk_var,pkwk_mode_perc,
#            pkwk_mean,leadpkwk_mean,delta_pkwk_mean)); # original
  out6=cbind(rep(fc_start,ntrn),tstep[1:ntrn],t(xprior_mean)); # This probably still has values where data is NA, but we're not using this so far, so it probably doesn't matter; also note that we'd have to remove one FORWARD from missing data
  out7 = formatDistENDS(onsets5Dist, endsDist, dursDist, peakWeeksDist, peakIntensitiesDist, nextILIDist)
  out7 = cbind(rep(fc_start, nrow(out7)), out7)
  colnames(out1)=c("fc_start","time","week","S","I","L","D","R0max","R0min","Est","S_sd","I_sd",
                   "L_sd","D_sd","R0max_sd","R0min_sd","Est_sd");
  colnames(out3)=c("fc_start","time","week","S","I","Est","S_sd","I_sd","Est_sd");
  colnames(out5)=c('fc_start','pkwk_mode','pkwk_mean','leadpkwk_mode','leadpkwk_mean','pkwk_sd',
                   'peak_intensity','peak_intensity_sd','tot_attack','fcast_1week','fcast_2week',
                   'fcast_3week','fcast_4week','rms','corr','onset5','onset5_sd','end','end_sd',
                   'duration','dur_sd')
  out8 = cbind(c('pi', '1week', '2week', '3week', '4week'), rep(fc_start, 5), rbind(peakIntensities, nextILI))
  rownames(out8) <- NULL
  colnames(out8) <- c('metrics', 'fc_start', 1:300)
  
  if(metricsonly==F){
    out=list(train=out1,fcast=out3,metrics=out5,trainprior=out6,dist=out7,ensembles=out8);
  }else{
    out=list(metrics=out5);
  }
}
