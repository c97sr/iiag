
# ********************************************************************************************************************************** #
### Functions for data formatting:

format.virologic.data <- function(d1) {
  # Calculates the percentage of tests positive for any influenza for each country and week
  
  colnames(d1) <- as.character(unlist(d1[3, ]))
  d1 <- d1[-c(1:3), ]
  rownames(d1) <- as.numeric(rownames(d1)) - 3
  
  d1[, 8] <- as.numeric(as.character(d1[, 8])) # number of specimens received
  d1[, 9] <- as.numeric(as.character(d1[, 9])) # number of specimens processed
  d1[, 20] <- as.numeric(as.character(d1[, 20])) # number of influenza positive results
  
  d1[, 9][is.na(d1[, 9]) & !is.na(d1[, 8])] <- d1[, 8][is.na(d1[, 9]) & !is.na(d1[, 8])]
  
  d1[, 9][d1[, 9] == 0 & !is.na(d1[, 9])] <- NA
  d1$posprop <- d1[, 20] / d1[, 9]
  d1$posprop[d1$posprop > 1 & !is.na(d1$posprop)] <- NA
  
  return(d1)
}

format.syndromic.data <- function(dnames, season) {
  # Reads in and compiles syndromic data, and removes age- and subtype-specific data
  #
  # Args:
  #   dnames: A vector containing the names of all syndromic data files
  #   season: The season for which influena is being forecasted
  #
  # Returns:
  #   d1: A dataframe including syndromic data from all countries
  
  data.list <- list()
  
  d1 <- NULL
  for (i in 1:length(dnames)) {
    filename <- dnames[i]
    d.temp <- read.csv(paste0(dir.WHO.data, filename), na.strings = '')
    
    d.temp$MEASURE_CODE <- as.character(d.temp$MEASURE_CODE)
    d.temp$AGEGROUP_CODE <- as.character(d.temp$AGEGROUP_CODE)
    d.temp$DataValue <- as.character(d.temp$DataValue)
    d.temp$COUNTRY <- as.character(d.temp$COUNTRY)
    
    d.temp <- d.temp[!is.na(d.temp$DataValue),]
    
    d1 <- rbind(d1, d.temp)
  }; rm(i, d.temp)
  
  d1$MEASURE_CODE <- factor(d1$MEASURE_CODE)
  d1$AGEGROUP_CODE <- factor(d1$AGEGROUP_CODE)
  
  d1 <- d1[d1$AGEGROUP_CODE == 'ALL' & !is.na(d1$AGEGROUP_CODE), ] # remove age-specific data
  d1$AGEGROUP_CODE <- factor(d1$AGEGROUP_CODE)
  
  d1 <- d1[is.na(d1$VIRUSTYPE_CODE) & is.na(d1$VIRUS_SUB_CODE) & is.na(d1$INF_CONF_CODE) & is.na(d1$SAMPLED_CODE), ]
  # remove subtype/positivity-specific data
  
  d1 <- d1[d1$YEAR_CODE >= substr(season, 1, 4) & ((d1$YEAR_CODE == substr(season, 1, 4) & d1$WEEK_CODE >= 40) |
                                                     d1$YEAR_CODE == (as.numeric(substr(season, 1, 4)) + 1)), ]
  # remove any data from before the season of interest (season beginning in week 40)
  
  d1$DataValue <- as.numeric(d1$DataValue)
  
  return(d1)
}

select.data.type.syndromic <- function(d1, types) {
  # Keeps only supported data types for each country
  #
  # Args:
  #   d1: A formatted dataframe containing syndromic data
  #   types: A dataframe of supported data types by country
  #
  # Returns:
  #   d1: A dataframe including only supported data types for each country
  
  rows.to.keep <- c()
  for (country in countries) {
    if (types$data.type[types$countries == country] == 'ILI') {
      rows.to.keep <- c(rows.to.keep, rownames(d1[d1$COUNTRY == country & d1$MEASURE_CODE == 'ILI CASES', ]))
    } else if (types$data.type[types$countries == country] == 'ARI') {
      rows.to.keep <- c(rows.to.keep, rownames(d1[d1$COUNTRY == country & d1$MEASURE_CODE == 'ARI CASES', ]))
    } else {
      print('ERROR: Data type not supported.')
    }
  }
  d1 <- d1[rows.to.keep, ]
  
  return(d1)
}

check.duplicate.data <- function(d1, d2) {
  # Checks syndromic and virologic data for duplicate rows; if duplicates
  # are identified, run is halted
  #
  # Args:
  #   d1: A formatted dataframe containing syndromic data
  #   d2: A formatted dataframe containing virologic data
  
  any.err <- FALSE
  for (country in countries) {
    for (year in unique(d1$YEAR_CODE)) {
      for (week in unique(d1$WEEK_CODE)) {
        
        syn.temp <- d1[d1$COUNTRY == country & d1$YEAR_CODE == year & d1$WEEK_CODE == week, ]
        vir.temp <- d2[d2$Country == country & d2$Year == year & d2$Week == week, ]
        
        if (length(syn.temp$COUNTRY) > 1) {
          print('ERROR: Duplicate syndromic data')
          any.err <- TRUE
        }
        if (length(vir.temp$Country) > 1) {
          print('ERROR: Duplicate virologic data')
          any.err <- TRUE
        }
        
      }
    }
  }
  
  return(any.err)
}

multiply.syn.vir <- function(d1, d2) {
  # Multiplies syndromic data by proportion of tests positive for influenza to
  # produce "syndromic+" data
  #
  # Args:
  #   d1: A formatted dataframe containing syndromic data
  #   d2: A formatted dataframe containing virologic data
  #
  # Returns:
  #   d1: A dataframe including only country, year, week, and syndromic+ data
  
  d1$data.plus <- NA
  
  for (country in levels(d1$COUNTRY)) {
    for(year in unique(d1$YEAR_CODE)) {
      for (week in unique(d1$WEEK_CODE)) {
        syn.temp <- d1[d1$COUNTRY == country & d1$YEAR_CODE == year & d1$WEEK_CODE == week, 'DataValue']
        vir.temp <- d2[d2$Country == country & d2$Year == year & d2$Week == week, 'posprop']
        
        if (length(syn.temp) == 1 & length(vir.temp) == 1) {
          d1$data.plus[d1$COUNTRY == country & d1$YEAR_CODE == year & d1$WEEK_CODE == week] <-
            syn.temp * vir.temp
        }
        if (length(syn.temp) > 1 | length(vir.temp) > 1) {
          print('ERROR: Duplicate data may exist')
        }
      }
    }
  }
  
  return(d1[, c(1:3, 11)])
}

scale.data <- function(d1, s) {
  # Scales data using country-specific scaling factors
  # For information on how scaling factors were determined, see:
  # Kramer & Shaman (Accepted), PLoS Computational Biology
  #
  # Args:
  #   d1: A formatted dataframe containing syndromic+ data
  #   s: A dataframe of scaling values
  # Returns:
  #   d1: A formatted dataframe containing scaled syndromic+ data
  
  for (country in levels(d1$COUNTRY)) {
    d1$data.plus[d1$COUNTRY == country] <-
      d1$data.plus[d1$COUNTRY == country] * s$scaling[s$countries == country]
  }
  return(d1)
}

flip.data <- function(d1) {
  # Reformats data so that each column represents data over time for a single country
  
  d1$WEEK_CODE <- as.factor(d1$WEEK_CODE)
  levels(d1$WEEK_CODE) <- str_pad(levels(d1$WEEK_CODE), 2, pad = '0')
  d1$time <- paste(d1$YEAR_CODE, d1$WEEK_CODE, sep = '_')
  d1$time <- factor(d1$time)
  
  flip.d1 <- dcast(d1, time ~ COUNTRY, fill = (-1), value.var = 'data.plus')
  flip.d1$time <- factor(flip.d1$time)
  
  return(flip.d1)
}

# ********************************************************************************************************************************** #
### Functions used in model fitting and forecasting:
replace.leading.na <- function(dat) {
  # Replaces any leading NAs with 0
  
  start.index <- 1
  while (dat[start.index] < 1 | is.na(dat[start.index])) {
    start.index <- start.index + 1
  }
  start.index <- start.index - 1
  if (start.index > 0) {
    dat[1:start.index] <- 0
  }
  
  return(dat)
}

calculate.oev <- function(dat, oev.denom) {
  # Calculates observational error variance using country-level data as an input
  
  tmp <- rep(0, length(dat))
  for (i in 4:length(dat)) {
    tmp[i] <- mean(dat[(i - 3):(i - 1)], na.rm = TRUE)
  }
  obs.vars <- (1e5 + (tmp^2)/5)/oev.denom
  
  for(i in which(is.na(obs.vars))) {
    dat.before <- dat[i-4]
    k <- i
    while(is.na(dat.before)) {
      k <- k - 1
      dat.before <- dat[k]
    }
    dat.after <- dat[i]
    j <- i
    while(is.na(dat.after)) {
      j <- j + 1
      if (j > length(dat)) {
        break
      }
      dat.after <- dat[j]
    }
    
    art.mean <- mean(c(dat.before, dat.after))
    obs.vars[i] <- (1e5 + (art.mean^2)/5)/oev.denom
  }
  
  return(obs.vars)
}

# ********************************************************************************************************************************** #
reformat.ens <- function(ens) {
  # Uses the individual ensemble member peak intensities to determine probability distributions
  # for peak timing standardized by historical maximum cases, and reformats these results to be
  # entered into outputDist
  
  ens <- as.data.frame(ens)
  ens <- ens[ens$metric == 'pi', ]
  
  # un-scale:
  for (i in c(4, 9:308)) {
    ens[, i] <- as.numeric(as.character(ens[, i]))
  }
  ens[, c(9:308)] <- (ens[, c(9:308)] / ens[, 4])
  
  # standardize by historical max:
  ens <- merge(ens, hist.max, by = 'country')
  ens[, c(9:308)] <- (ens[, c(9:308)] / ens[, 309])
  
  # reformat to resemble outputDist
  reqLimits <- seq(0, 1, by = 0.1)
  
  new.dists <- lapply(1:length(ens$country), function(ix) {
    temp.dat <- as.vector(ens[ix, 9:308])
    peakIntensitiesDist = matrix(NA, nrow = length(reqLimits), ncol = 8)
    peakIntensitiesDist[, 1] <- ens[ix, 1]
    peakIntensitiesDist[, 2] <- ens[ix, 2]
    peakIntensitiesDist[, 3] <- ens[ix, 3]
    peakIntensitiesDist[, 4] <- ens[ix, 4]
    peakIntensitiesDist[, 5] <- ens[ix, 6]
    peakIntensitiesDist[, 6] <- ens[ix, 5]
    row = 1
    for (i in 2:length(reqLimits)) {
      peakIntensitiesDist[row, 7] <- reqLimits[i]
      peakIntensitiesDist[row, 8] <- round(length(temp.dat[temp.dat >= reqLimits[i - 1] & temp.dat < reqLimits[i]]) / length(temp.dat), 4)
      row <- row + 1
    }
    peakIntensitiesDist[row, 7] <- 1.1
    peakIntensitiesDist[row, 8] <- round(length(temp.dat[temp.dat >= max(reqLimits)]) / length(temp.dat), 4)
    colnames(peakIntensitiesDist) <- c('country', 'season', 'run', 'scaling', 'fc_start', 'metric', 'estimate', 'prob')
    return(peakIntensitiesDist)
  })
  new.dists <- as.data.frame(do.call(rbind, new.dists))
  
  new.dists$country <- levels(ens$country)[new.dists$country]
  new.dists$season <- '2018-19'
  new.dists$metric <- 'pi'
  
  return(new.dists)
}
