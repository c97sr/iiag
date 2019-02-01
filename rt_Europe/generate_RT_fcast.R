
### Code to read and format data, generate real-time forecasts, and output formatted results
library(reshape2)
library(stringr)
require(plyr)

# Before running forecasts, you will have to download FluNet and FluID data for the countries
# of your choosing, then edit the code below to point to where the data are stored. Supported
# countries are:
countries.supported <- c('Austria', 'Belarus', 'Belgium', 'Bulgaria', 'Croatia', 'Czechia', 'Denmark', 'Estonia',
                         'Finland', 'France', 'Georgia', 'Germany', 'Greece', 'Hungary', 'Iceland', 'Ireland',
                         'Israel', 'Italy', 'Kazakhstan', 'Kyrgyzstan', 'Latvia', 'Lithuania', 'Luxembourg',
                         'Netherlands', 'Norway', 'Poland', 'Portugal', 'Republic of Moldova', 'Romania',
                         'Russian Federation', 'Serbia', 'Slovakia', 'Slovenia', 'Spain', 'Sweden', 'Switzerland',
                         'Turkey', 'Ukraine', 'United Kingdom of Great Britain and Northern Ireland', 'Uzbekistan',
                         'Canada', 'Mexico', 'Morocco')

# ********************************************************************************************************************************** #
### EDIT THIS SECTION ###
dir.WHO.data <- '/WHOdata/' # directory where influenza data are stored
dir.data <- '/data/' # directory where other data are stored
dir.code <- '/code/' # directory where code is stored
dir.results <- '/results/' # directory to store results

### Season and week of forecast
season <- '2018-19'
fcast.week <- 55
# Can also edit file "fcast_setup.R" to change parameter ranges, number of runs/ensemble members, etc.

# ********************************************************************************************************************************** #
### Read in and format data ###
source(paste0(dir.code, 'kramer_funcs.R'))

# Initial formatting:
vir.dat <- read.csv(paste0(dir.WHO.data, 'FluNetInteractiveReport.csv'), header = F)
vir.dat <- format.virologic.data(vir.dat)

syn.dat.list <- list.files(dir.WHO.data, pattern = 'FLUID')
syn.dat <- format.syndromic.data(syn.dat.list, season); rm(syn.dat.list)

# Get vector of countries being forecasted:
countries <- unique(syn.dat$COUNTRY)
countries <- countries[which(countries %in% countries.supported)]

# Keep only relevant data types:
data.types <- read.csv(paste0(dir.data, 'data_types.csv'))
syn.dat <- select.data.type.syndromic(syn.dat, data.types)

# Keep only countries with data:
countries <- countries[countries %in% unique(vir.dat$Country)]
vir.dat <- vir.dat[vir.dat$Country %in% countries, ]; vir.dat$Country <- factor(vir.dat$Country)
syn.dat <- syn.dat[syn.dat$COUNTRY %in% countries, ]; syn.dat$COUNTRY <- factor(syn.dat$COUNTRY)

# Check for duplicate data:
if (check.duplicate.data(syn.dat, vir.dat)) {
  stop('Data contain duplicate rows. Please fix this issue before continuing.')
}

# Multiply syndromic and virologic data:
comb.dat <- multiply.syn.vir(syn.dat, vir.dat)

# Scale data:
scalings <- read.csv(paste0(dir.data, 'scaling_values.csv'))
comb.dat <- scale.data(comb.dat, scalings)

# Remove country/year/week combinations w/ no data:
comb.dat <- comb.dat[!is.na(comb.dat$data.plus), ]
comb.dat$COUNTRY <- factor(comb.dat$COUNTRY)

# Format data for forecasting:
comb.dat <- flip.data(comb.dat)

# ********************************************************************************************************************************** #
### Perform the real-time forecasts ###
ntrn <- fcast.week - 40 + 1 # number of observations for the training period; could be 1 to 51 weeks
source(paste0(dir.code, 'fcast_setup.R'))

# Initialize output data frames:
outputMetrics <- NULL # contains measures such as peak timing, peak intensity, etc.
outputOP <- NULL # contains values for fitted model states and parameters
outputDist <- NULL # contains probability distributions for metrics
outputEns <- NULL # contains intensity predictions (peak and 1-4 week) for each ensemble member

for (i in 1:length(countries)) {
  country <- countries[i]
  
  # Get country-specific humidity data:
  AH <- matrix(c(ah[, country], ah[, country]), 365*2, 1)
  
  # Get country-specific influenza data:
  obs_i <- as.matrix(comb.dat[, country], dim(comb.dat)[1], dim(comb.dat)[2])
  obs_i[obs_i == -1] <- NA
  
  # Main loop:
  if(!all(is.na(obs_i)) & any(obs_i > 1, na.rm=T) & !is.na(obs_i[ntrn]) & obs_i[ntrn] > 0) {
    print(paste0('Running forecast for ', country, '...'))
    
    obs_i <- replace.leading.na(obs_i) # convert any leading NAs to 0
    obs_vars <- calculate.oev(obs_i, oev_denom) # calculate estimated "observed" error variance
    scale.temp <- scalings$scaling[scalings$countries == country] # record country's scaling factor
    
    for (run in 1:num_runs) {
      res <- EAKF_rFC(num_ens, tmstep, param.bound, obs_i = obs_i, ntrn, obs_vars,
                      tm.ini, tm.range)
      
      outputMetrics <- rbind(outputMetrics, cbind(country, season, run, scale.temp, res$metrics))
      outputDist <- rbind(outputDist, cbind(country, season, run, scale.temp, res$dist))
      outputEns <- rbind(outputEns, cbind(country, season, run, scale.temp, oev_denom, lambda, res$ensembles))
      
      for (j in 1:2) {
        temp_res <- buildOutput(res[[j]], 'long') # reorder columns
        temp_res <- convertDaysToDate(temp_res, start_date) #convert 'time' to date
        obs_output <- NULL
        
        if (names(res)[j] == "fcast"){
          if (length(obs_i) == nsn){
            obs_output <- obs_i[(ntrn + 1):nsn]
          }
          
          else {
            if (ntrn + 1 <= length(obs_i)) {
              obs_output <- obs_i[(ntrn + 1):length(obs_i)]
            }
            
            obs_output <- c(obs_output, rep(-1, nsn - length(obs_i)))
          }
        } else {
          obs_output <- obs_i[1:ntrn][!is.na(obs_i[1:ntrn])]
        }
        
        outputOP <- rbind(outputOP, cbind(country, season, run, scale.temp, names(res)[i], temp_res, obs_output))
      }
      
    }
  
  }
}

# Rename output file columns
dimnames(outputMetrics) <- list(c(), as.list(metrics_header))
dimnames(outputOP)[[2]] <- as.list(output_header)
dimnames(outputDist) <- list(c(), as.list(dist_header))
dimnames(outputEns)[[2]][1:7] <- as.list(ens_header)

# Reformat dist file to match format on website
# https://cpid.iri.columbia.edu/
# Intensity reported as percent of historical maximum, not as scaled or unscaled cases
hist.max <- read.csv(paste0(dir.data, 'historical_maximums.csv'))
outputDistPerc <- reformat.ens(outputEns)

outputDist <- as.data.frame(outputDist)
outputDist <- outputDist[outputDist$metric != 'pi', ]
for (i in c(3:5, 7:8)) {
  outputDist[, i] <- as.numeric(as.character(outputDist[, i]))
}
outputDist <- rbind(outputDist, outputDistPerc)

# Save results files
write.csv(outputMetrics, paste0(dir.results, 'outputMet_Week', fcast.week, '.csv'), row.names = FALSE)
write.csv(outputOP, paste0(dir.results, 'outputOP_Week', fcast.week, '.csv'), row.names = FALSE)
write.csv(outputDist, paste0(dir.results, 'outputDist_Week', fcast.week, '.csv'), row.names = FALSE)
write.csv(outputEns, paste0(dir.results, 'outputEns_Week', fcast.week, '.csv'), row.names = FALSE) # note: large file
rm(list=ls())
