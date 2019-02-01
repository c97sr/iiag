## Code to read in all functions and set all parameters necessary for forecast generation

### Read in model functions
source(paste0(dir.code,"EpiModelsTrackIncidence.R")) # SIRS model
source(paste0(dir.code,"Fn_checkxnobounds.R")) # function to check DA aphysicality
source(paste0(dir.code,'Util.R')) #buildOutput and ConvertDaysToDate functions; simple utilities
source(paste0(dir.code,'Fn_initializations.R')) # source of some funtions for initializations
source(paste0(dir.code,'Fn_EAKF_rFC.R'))

### Specify headers for output file
metrics_header <- c('country', 'season', 'run', 'scaling', 'fc_start', 'pkwk_mode',
                    'pkwk_mean', 'leadpkwk_mode', 'leadpkwk_mean', 'pkwk_sd', 'peak_intensity',
                    'peak_intensity_sd', 'tot_attack', 'fcast_1week', 'fcast_2week',
                    'fcast_3week', 'fcast_4week', 'rms', 'corr', 'onset', 'onset_sd',
                    'end', 'end_sd', 'duration', 'dur_sd')
output_header <- c('country', 'season', 'run', 'scaling', 'result', 'fc_start', 'time', 'time_date',
                   'week', 'estimate', 'est_sd', 'S', 'S_sd', 'I', 'I_sd', 'L', 'L_sd', 'D', 'D_sd',
                   'R0max', 'R0max_sd', 'R0min', 'R0min_sd', 'observed')
dist_header <- c('country', 'season', 'run', 'scaling', 'fc_start', 'metric', 'estimate', 'prob')
ens_header <- c('country', 'season', 'run', 'scaling', 'oev', 'lambda', 'metric')

### Read in list of variables
if(!exists("phi")) {
  phi <- read.csv(paste(dir.data, 'phi.csv', sep=''), header=F)
  params <- phi[,1:4]
  susceps = infects = NULL
  for(i in 5:31) {
    susceps <- append(susceps, phi[,i])
    infects <- append(infects, phi[,27+i])
  }
}

### Global variables
N <- 1e5 # population
dt <- 1 # time step for SIRS integration
tmstep <- 7 #data is weekly

### Parameter boundaries
D_low <- 1.5; L_low <- 1*365; Rmx_low <- 1.3; Rmn_low <- 0.8;
D_up <- 7; L_up <- 10*365; Rmx_up <- 4; Rmn_up <- 1.2;
theta_low <- c(L_low, D_low, Rmx_low, Rmn_low)
theta_up <- c(L_up, D_up, Rmx_up, Rmn_up)
param.bound <- cbind(theta_low, theta_up)

### Parameters for the filters
discrete <- FALSE # run the SIRS model continuously
metricsonly <- FALSE # save all outputs
lambda <- 1.03 # inflation factor for the ensemble filters
oev_denom <- 1 #10 #c(1,10,100) # denominator for observation error variance
num_ens <- 300 # use 300 for ensemble filters, 10000 for particle filters
num_runs <- 5
wk_start <- 40

### Read in humidity data
ah <- read.csv(paste(dir.data, 'ah_realtime.csv', sep=''), header=T)
colnames(ah)[c(28, 30, 39)] <- c('Republic of Moldova', 'Russian Federation',
                                 'United Kingdom of Great Britain and Northern Ireland')

### Get dates/times for the simulation
tmp <- Fn_dates(season)
weeks <- 1:ntrn
start_date <- tmp$start_date
end_date <- tmp$end_date
nsn <- tmp$nsn

clim_start <- as.numeric(start_date - as.Date(paste('20',
                                                    substr(season, gregexpr('-', season)[[1]][1]-2,
                                                           gregexpr('-', season)[[1]][1]-1),
                                                    '-01-01', sep=''))) + 1 - 6
clim_end <- as.numeric(end_date - as.Date(paste('20',
                                                substr(season, gregexpr('-', season)[[1]][1]-2,
                                                       gregexpr('-', season)[[1]][1]-1),
                                                '-01-01', sep=''))) + 1
tm.ini <- clim_start - 1 # time to initiate model
tm.range <- clim_start:clim_end # time span of model run
