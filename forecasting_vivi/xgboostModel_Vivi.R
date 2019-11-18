library('idd') # WHO influenza-likeness-illness data; idd package is not available (for R version 3.6.1)
library('ggplot2') # transform integer numbers to categorical variable
library('xgboost') # boosted regression tree
library('formattable') # formatting on data frames
library('gtools') # permutations
library("lattice") # heat plot
library("RColorBrewer") # colour platter
library("raster") # colour platter
library("rasterVis") # colour platter
# library("purrr") # Partial dependence plots
library("tidyverse")
library("dplyr") # mutate
library("lubridate")
library("mlr")
library("Matrix")
library("e1071")
library("DiagrammeR") # tree plot
library("rsvg")
library("DiagrammeRsvg") # multi tree plot
library("data.table") # convert data frame to data table
library("aweek") # convert discrete time variable into contious variable
library("knitr")

rm(list = ls(all = TRUE))

# Load WHO FluID dataset
#' Becuase idd package is not aviable for R version 3.6.1, I downloaded dataset fluIliCountryData.rda in the 
#' idd and loaded in the Rstudio directly.

load(file = "fluIliCountryData.rda")

# data("fluIliCountryData")

#' Load country list for WHO data.
countryISO <- read.csv("C:/Users/hw3616/Desktop/Imperial/Project1_Forecasting/Project_Coding/iiag/data_old/country_list_ISO.csv")

#' Load WHO FluID data set.
#' From 2010 week 1 to 2018 week 9 
load.iiag.data.fluid <- function(datadir="../iiag_data/data_old/") {
  
  ## Helper function to fix some header names to be used below
  fix_headers <- function(x) {
    curnames <- names(x)
    newnames <- curnames
    if (!is.na(match("ISO_Week",curnames))) {
      newnames[match("ISO_Week",curnames)] <- "ISO_WEEK"
    }
    if (!is.na(match("ï..ISO3",curnames))) {
      newnames[match("ï..ISO3",curnames)] <- "ISO3"
    }
    newnames
  }

  ## Define the strings for all the files needed and read 
  # fid_this <- read.csv(paste(datadir,"/2019-2020_FluIDData.csv",sep=""))
  fid_old_3 <- read.csv(paste0(datadir,"/2017-2018_FluIDData.csv",sep=""))
  fid_old_2 <- read.csv(paste0(datadir,"/2014-2016_FluIDData.csv",sep=""))
  fid_old_1 <- read.csv(paste0(datadir,"/2010-2013_FluIDData.csv",sep=""))
  fid_old_0 <- read.csv(paste0(datadir,"/2000-2009_FluIDData.csv",sep=""))
  # fnet_this <- read.csv(paste(datadir,"/2017-2018_FluNetData_20190110.csv",sep=""))
  # fnet_old_3 <- read.csv(paste(datadir,"/2017-2018_FluNetData.csv",sep=""))
  # fnet_old_2 <- read.csv(paste(datadir,"/2014-2016_FluNetData.csv",sep=""))
  # fnet_old_1 <- read.csv(paste(datadir,"/2010-2013_FluNetData.csv",sep=""))
  # fnet_old_0 <- read.csv(paste(datadir,"/2000-2009_FluNetData.csv",sep=""))
  
  ## Fix names for flu id
  # names(fid_this) <- fix_headers(fid_this)
  names(fid_old_1) <- fix_headers(fid_old_1)
  names(fid_old_2) <- fix_headers(fid_old_2)
  names(fid_old_3) <- fix_headers(fid_old_3)
  # names(fnet_old_1) <- fix_headers(fnet_old_1)
  # names(fnet_old_2) <- fix_headers(fnet_old_2)
  # names(fnet_old_3) <- fix_headers(fnet_old_3)
  # names(fnet_this) <- fix_headers(fnet_this)

  ## Use rbind to make the large tables. Should throw an error if the column
  ## names change in the future.
  # dfId <- rbind(fid_old_0,fid_old_2,fid_old_2,fid_old_3,fid_this)
  # dfNet <- rbind(fnet_old_0,fnet_old_2,fnet_old_2,fnet_old_3,fnet_this)
  dfId <- rbind(fid_old_1,fid_old_2,fid_old_3)
  # dfNet <- rbind(fnet_old_1,fnet_old_2,fnet_old_3,fnet_this)

  ## Sort both dataframe just incase
  dfId <- dfId[order(dfId$ISO2,dfId$ISO_YEAR,dfId$ISO_WEEK),]
  # dfNet <- dfNet[order(dfNet$ISO2,dfNet$ISO_YEAR,dfNet$ISO_WEEK),]
  
  ## Return the two datasets as a list
  # list(lab=dfNet,synd=dfId)
  dfId
}

# WHO dataset
fluWHO <- load.iiag.data.fluid(datadir = "C:/Users/hw3616/Desktop/Imperial/Project1_Forecasting/Project_Coding/iiag/data_old")

## Take the raw WHO country database and extract weekly incidence
extract.incidence.who <- function( dfId,
                                   sel_iso3,
                                   sel_ag,
                                   sel_measure,
                                   minYear,
                                   maxYear) {
  
  ## Setup the week scale in a format consistent with the week format
  ## in the data and cope with 53-week years. Needs the list of 53 week years
  ## extending in both directions.
  ## Perhaps should have a few lines to get rid of data NAs and avoid a warning
  ## at the next line?
  dfId$yrweek <- paste(dfId$ISO_YEAR,sprintf("%02d",as.numeric(dfId$ISO_WEEK)),sep="-")
  min(dfId$ISO_YEAR)
  yrs53Weeks <- 2015
  currentYear <- minYear
  vecWeekScale <- NULL
  while (currentYear <= maxYear) {
    if (currentYear %in% yrs53Weeks) {
      max_week <- 53
    } else {
      max_week <- 52
    }
    vecWeekScale <- c(vecWeekScale,
                      paste(currentYear,sprintf("%02d",1:max_week),sep="-"))
    currentYear <- currentYear +1
  }
  
  ## Define the return matrix for the function
  sel_weeks <- vecWeekScale
  rtnmat <- matrix(data=NA,nrow=length(sel_weeks),ncol=length(sel_iso3))
  colnames(rtnmat) <- sel_iso3
  rownames(rtnmat) <- sel_weeks
  
  ## Start outer loop over the country codes
  for (cur_iso3 in sel_iso3) {
    
    ## Define criteria and subset the data
    crit1 <- (dfId$ISO3 == cur_iso3)
    if(!("AGEGROUP_CODE" %in% colnames(dfId))) {
      crit2 <- TRUE
    } else {
      crit2 <- (dfId$AGEGROUP_CODE %in% sel_ag)
    }
    
    crit3 <- (dfId$MEASURE_CODE %in% sel_measure)
    tmpdf <- dfId[crit1 & crit2 & crit3,]
    tmpdf <- tmpdf[order(tmpdf$yrweek),]
    
    ## Setup the preconditions for the nested while loops
    max_ind_rtn <- dim(rtnmat)[1]
    max_ind_df <- dim(tmpdf)[1]
    cur_ind_rtn <- 1
    cur_ind_df <- 1
    
    ## 2-level while loop with index "pointers" into the rtn matrix
    ## and the dataframe. Scans through the data and the rtn matrix
    ## at the same time and adds any none-na value that meets the
    ## criteria for any given week. This works only because the
    ## date format is correctly ordered by sort even though its not
    ## a numeric and the subsetted dataframe _has_ been sorted.
    ## Could be done with a small number of table commands, but I
    ## (SR) wanted to be able to handle any line-by-line cleaning
    ## in future within this loop if needed.
    while (cur_ind_df <= max_ind_df) {
      while (
        sel_weeks[cur_ind_rtn] != tmpdf$yrweek[cur_ind_df] &&
        cur_ind_rtn <= max_ind_rtn
      ) {
        cur_ind_rtn <- cur_ind_rtn + 1
      }
      if (cur_ind_rtn <= max_ind_rtn) {
        val_rtn <- rtnmat[cur_ind_rtn,cur_iso3]
        val_df <- as.numeric(tmpdf$ValueNumeric[cur_ind_df])
        if (!is.na(val_df)) {
          if (is.na(val_rtn)) {
            rtnmat[cur_ind_rtn,cur_iso3] <- val_df
          } else {
            rtnmat[cur_ind_rtn,cur_iso3] <-
              rtnmat[cur_ind_rtn,cur_iso3] + val_df
          }
        }
      }
      cur_ind_df <- cur_ind_df + 1
    }
    
    ## Close the country-level loop
  }
  
  ## Return the populated incidence matrix as only result of function
  rtnmat
  
}


#' I will only extract countries appearred in my report for now to keep consistency
#' To get the country list, I will run the code of filtering countries.

minprop <- 0.5
fluWHO.incidence <- extract.incidence.who(
  fluWHO,
  sel_iso3 = unique(fluWHO$ISO3),
  sel_ag = c("All"),
  sel_measure = c("ILI_CASES"),
  minYear=2010,
  maxYear = 2018)

sel_iso <- names(which(colSums(is.na(fluWHO.incidence))/dim(fluWHO.incidence)[1]<minprop))

fluWHO.incidence <- extract.incidence.who(
  fluWHO,
  sel_iso3 = sel_iso,
  sel_ag = c("All"),
  sel_measure = c("ILI_CASES"),
  minYear=2010,
  maxYear=2018
)

#' Check if the data of left countries is eligible for xgboost

#' extract_incidence is the function in package idd which can'y be loaded
#' therefore copy the code 
extract_incidence.idd <- function(flu_data,
                                  country_code,
                                  year) {
  flu_data <- as.data.frame(flu_data)
  year_names <- rownames(flu_data)
  # start plotting at week 27 of the current year
  row_name_start <- paste0(year, "-27") 
  # stop plotting at week 26 of the next year
  row_name_end <- paste0(year + 1, "-26")
  # find the corresponding weeks in the data
  row_index_start <- which(rownames(flu_data) == row_name_start)
  row_index_end <- which(rownames(flu_data) == row_name_end)
  # extrac the week number and incidence for those weeks
  incidence <- flu_data[seq(row_index_start, row_index_end), 
                        colnames(flu_data) == country_code]
  time_name_vec <- year_names[seq(row_index_start, row_index_end)]
  
  incidence_data <- data.frame(t = seq_along(time_name_vec), 
                               time_name = time_name_vec, 
                               incidence = incidence)
  return(incidence_data)
}

#' XGBoost uses a matrix of input data instead of a data frame,
#' so the output of gbm_complex should be a matrix of dataset.
gbm_complex <- function(data, country, num_category,nWeek_ahead){
  yr <- seq(2010, 2017, by = 1)
  initial_data <- c()
  for (i in 1:length(yr)){
    tmp <- extract_incidence.idd(data, country_code = country, yr[i])
    initial_data <- rbind(initial_data, tmp)
  }
  initial_data <- as.data.frame(initial_data)
  
  # Divide incidence into 10 categories
  initial_data2 <- cbind(initial_data, cut_interval(initial_data$incidence, n=num_category))
  colnames(initial_data2)[4] <- 'category'
  
  # match the intervals with numeric categories
  levels_index <- levels(initial_data2$category)
  for(i in 1:nrow(initial_data2)){
    x <- initial_data2$category[i]
    if (is.na(x)==FALSE){
      y <- which(levels(initial_data2$category)==x)
      initial_data2$category2[i] <- as.numeric(y)
    }else{
      initial_data2$category2[i] <- NA
    }
  }
  
  # convert into suitable stucture for gbm
  incidence_gbm <- matrix(NA, nrow = (dim(initial_data2)[1]-5), ncol=3)
  
  if(nWeek_ahead == 1){
    for (i in 1:nrow(incidence_gbm)){
      incidence_gbm[i,] <- initial_data2$category2[c((i+5),(i+4),(i+3))]
    }
    
    incidence_gbm <- as.data.frame(incidence_gbm)
    colnames(incidence_gbm) <- c('Y_week0','week_1','week_2')
    rownames(incidence_gbm) <- as.character(initial_data2$time_name[6:nrow(initial_data2)])
  }
  
  if(nWeek_ahead == 2){
    for (i in 1:nrow(incidence_gbm)){
      incidence_gbm[i,] <- initial_data2$category2[c((i+5),(i+3),(i+2))]
    }
    
    incidence_gbm <- as.data.frame(incidence_gbm)
    colnames(incidence_gbm) <- c('Y_week0','week_2','week_3')
    rownames(incidence_gbm) <- as.character(initial_data2$time_name[6:nrow(initial_data2)])
  }
  
  if(nWeek_ahead == 3){
    for (i in 1:nrow(incidence_gbm)){
      incidence_gbm[i,] <- initial_data2$category2[c((i+5),(i+2),(i+1))]
    }
    
    incidence_gbm <- as.data.frame(incidence_gbm)
    colnames(incidence_gbm) <- c('Y_week0','week_3','week_4')
    rownames(incidence_gbm) <- as.character(initial_data2$time_name[6:nrow(initial_data2)])
  }
  
  
  if(nWeek_ahead == 4){
    for (i in 1:nrow(incidence_gbm)){
      incidence_gbm[i,] <- initial_data2$category2[c((i+5),(i+1),i)]
    }
    
    incidence_gbm <- as.data.frame(incidence_gbm)
    colnames(incidence_gbm) <- c('Y_week0','week_4','week_5')
    rownames(incidence_gbm) <- as.character(initial_data2$time_name[6:nrow(initial_data2)])
  }
  
  
  # add potential covariates
  # month: Jan - Dec
  # incidence_gbm$week_time <- as.character(initial_data2$time_name[3:nrow(initial_data2)])
  week_time <- rownames(incidence_gbm)
  incidence_gbm$month <- NA
  
  jan_row <- c()
  jan_week <- paste("-0",1:5, sep = "")
  for (i in 1:length(jan_week)){
    tmp <- grep(jan_week[i], rownames(incidence_gbm))
    jan_row <- append(jan_row, tmp)
  }
  incidence_gbm$month[jan_row] <- "January"
  
  feb_row <- c()
  feb_week <- paste("-0",6:9, sep = "")
  for (i in 1:length(feb_week)){
    tmp <- grep(feb_week[i], rownames(incidence_gbm))
    feb_row <- append(feb_row, tmp)
  }
  incidence_gbm$month[feb_row] <- "Feburary"
  
  mar_row <- c()
  mar_week <- paste("-", 10:14, sep = "")
  for (i in 1:length(mar_week)){
    tmp <- grep(mar_week[i], rownames(incidence_gbm))
    mar_row <- append(mar_row, tmp)
  }
  incidence_gbm$month[mar_row] <- "March"
  
  apr_row <- c()
  apr_week <- paste("-", 15:18, sep = "")
  for (i in 1:length(apr_week)){
    tmp <- grep(apr_week[i], rownames(incidence_gbm))
    apr_row <- append(apr_row, tmp)
  }
  incidence_gbm$month[apr_row] <- "April"
  
  may_row <- c()
  may_week <- paste("-", 19:23, sep = "")
  for (i in 1:length(may_week)){
    tmp <- grep(may_week[i], rownames(incidence_gbm))
    may_row <- append(may_row, tmp)
  }
  incidence_gbm$month[may_row] <- "May"
  
  june_row <- c()
  june_week <- paste("-", 24:27, sep = "")
  for (i in 1:length(june_week)){
    tmp <- grep(june_week[i], rownames(incidence_gbm))
    june_row <- append(june_row, tmp)
  }
  incidence_gbm$month[june_row] <- "June"
  
  july_row <- c()
  july_week <- paste("-", 28:32, sep = "")
  for (i in 1:length(july_week)){
    tmp <- grep(july_week[i], rownames(incidence_gbm))
    july_row <- append(july_row, tmp)
  }
  incidence_gbm$month[july_row] <- "July"
  
  aug_row <- c()
  aug_week <- paste("-", 33:36, sep = "")
  for (i in 1:length(aug_week)){
    tmp <- grep(aug_week[i], rownames(incidence_gbm))
    aug_row <- append(aug_row, tmp)
  }
  incidence_gbm$month[aug_row] <- "August"
  
  sep_row <- c()
  sep_week <- paste("-", 37:40, sep = "")
  for (i in 1:length(sep_week)){
    tmp <- grep(sep_week[i], rownames(incidence_gbm))
    sep_row <- append(sep_row, tmp)
  }
  incidence_gbm$month[sep_row] <- "September"
  
  oct_row <- c()
  oct_week <- paste("-", 41:44, sep = "")
  for (i in 1:length(oct_week)){
    tmp <- grep(oct_week[i], rownames(incidence_gbm))
    oct_row <- append(oct_row, tmp)
  }
  incidence_gbm$month[oct_row] <- "October"
  
  nov_row <- c()
  nov_week <- paste("-", 45:48, sep = "")
  for (i in 1:length(nov_week)){
    tmp <- grep(nov_week[i], rownames(incidence_gbm))
    nov_row <- append(nov_row, tmp)
  }
  incidence_gbm$month[nov_row] <- "November"
  
  dec_row <- c()
  dec_week <- paste("-", 49:52, sep = "")
  for (i in 1:length(dec_week)){
    tmp <- grep(dec_week[i], rownames(incidence_gbm))
    dec_row <- append(dec_row, tmp)
  }
  incidence_gbm$month[dec_row] <- "December"
  
  # season 
  incidence_gbm$season <- NA
  for(i in 1:nrow(incidence_gbm)){
    if(incidence_gbm$month[i] %in% c("September","October","November")==TRUE){
      incidence_gbm$season[i] <-"autumn"
    }
    if(incidence_gbm$month[i] %in% c("June","July","August")==TRUE){
      incidence_gbm$season[i] <- "summer"
    }
    if(incidence_gbm$month[i] %in% c("December","January","February")==TRUE){
      incidence_gbm$season[i] <- "winter"
    } 
    if (incidence_gbm$month[i] %in% c("March","April","May")==TRUE){
      incidence_gbm$season[i] <- "spring"
    }
  }
  
  # deal with NA data
  b <- c()
  for(i in 1:nrow(incidence_gbm)){
    if(sum(is.na(incidence_gbm[i,])==TRUE) >= 1){
      tmp <- i
      b <- append(b,tmp)
    }
    b
  }
  
  incidence_gbm <- incidence_gbm[-b, ]
  incidence_gbm <- as.data.frame(incidence_gbm)
  incidence_gbm
}


#' explore how many years data each countries owns
duration <- function(country){
  year_time <- c(2010:2018)
  
  flu_data_complex <- gbm_complex(fluWHO.incidence, country, 10,1)
  year_start <- min(as.numeric(substr(rownames(flu_data_complex),0,4)))
  year_end <- max(as.numeric(substr(rownames(flu_data_complex),0,4)))
  all_year <- as.numeric(substr(rownames(flu_data_complex),0,4))
  
  country_year <- c()
  
  for (i in 1:length(year_time)){
    if (year_time[i] %in% all_year == TRUE){
      tmp <- "Yes"
    }
    if (year_time[i] %in% all_year == FALSE){
      tmp <- "No"
    }
    country_year <- append(country_year, tmp)
  }
  country_year <- c(country, country_year, year_start, year_end)
  
  country_year
}

country_year <- NULL
for (i in 1:length(sel_iso)){
  tmp <- duration(sel_iso[i])
  country_year <- rbind(country_year, tmp)
}
country_year <- as.data.frame(country_year)
colnames(country_year) <- c("Country","2010","2011","2012","2013","2014","2015",
                            "2016","2017","2018","start_year","end_year")
country_year$end_year <- as.numeric(as.character(country_year$end_year))
country_year$start_year <- as.numeric(as.character(country_year$start_year))
rownames(country_year) <- c(1:nrow(country_year))

countryCode_no1718Or10 <- country_year$Country[which(country_year$end_year < 2017 | country_year$`2010`=="No" )]
countryCode_no10 <- country_year$Country[which(country_year$`2010`=="No")]
#' countries will not be used in xgboost model because of lack of 2017 and 2018 data or 2010 data
#' They are Barbados,Belarus,Bhutan, Honduras, New Zealand, Nigeria,Oman, Pakistan,Singapore,Tajikistan,
#' Thailand, Macedonia, the former Yugoslav Republic of, Northern Mariana Islands. 
country_no1718Or10 <- sel_iso[which(sel_iso %in%countryCode_no1718Or10)] # 20 in total

#' Exlcude the countries do not have data of 2010 or 2017 or 2018
#' 51 countries left
sel_iso_xgb <- sel_iso
for (i in 1:length(country_no1718Or10)){
  index <- which(sel_iso_xgb == country_no1718Or10[i])
  sel_iso_xgb <- sel_iso_xgb[-(index)]
}

#' Check different countries in fluIliCountryData with countries in null_models
country.idd <- colnames(fluIliCountryData)
country.null <- colnames(x)

same.country.idd.null <- c()
not.same.country.idd.null <- c()
for (i in 1:length(country.null)){
  if (country.null[i] %in% country.idd == TRUE){
    tmp <- country.null[i]
    same.country.idd.null <- append(same.country.idd.null, tmp)
  }else{
    tmp <- country.null[i]
    not.same.country.idd.null <- append(not.same.country.idd.null,country.null[i])
  }
}

#' The 68 countries in fluIliCountryData are also in the WHO dataset.
#' There are 3 more countries in WHO dataset, they are Cuba CUB, Iceland ISL and Italy ITA
not.same.country.idd.null


# for (i in 1:length(not.same.country.idd.null)){
#  if(not.same.country.idd.null[i] %in% sel_iso_xgb == TRUE){
#    index <- which(sel_iso_xgb == not.same.country.idd.null[i])
#  }
#  sel_iso_xgb <- sel_iso_xgb[-(index)]
#}


#' Check if countries have less than 5 weeks data in a year because I will do the 4-week ahead foreacast which
#' requires data of 5 consective weeks
#' Exclude countires whose datasets are uneligible to do the 4-week ahead forecast
#' GUM 2016 data only get one obeservation
sel_iso_xgb <- sel_iso_xgb[-which(sel_iso_xgb%in%"GUM")] # 50 countries
# KIR 2016 data only 5 observations,  4-week predict.
sel_iso_xgb <- sel_iso_xgb[-which(sel_iso_xgb%in%"KIR")] # 49 countries
# FSM
sel_iso_xgb <- sel_iso_xgb[-which(sel_iso_xgb%in%"FSM")] # 48 countries
# MHL
sel_iso_xgb <- sel_iso_xgb[-which(sel_iso_xgb%in%"MHL")] # 47 countries
# PLW
sel_iso_xgb <- sel_iso_xgb[-which(sel_iso_xgb%in%"PLW")] # 46 countries
# COK
sel_iso_xgb <- sel_iso_xgb[-which(sel_iso_xgb%in%"COK")] # 45 countries
# NRU
sel_iso_xgb <- sel_iso_xgb[-which(sel_iso_xgb%in%"NRU")] # 44 countries
# TUV
sel_iso_xgb <- sel_iso_xgb[-which(sel_iso_xgb%in%"TUV")] # 43 countries
# ASM
sel_iso_xgb <- sel_iso_xgb[-which(sel_iso_xgb%in%"ASM")] # 42 countries
#ISL only 4 consecutive week in 2010
sel_iso_xgb <- sel_iso_xgb[-which(sel_iso_xgb%in%"ISL")] # 41 countries

# extract incidence of 42 eligible countries.
fluWHO.incidence <- extract.incidence.who(fluWHO,
                                          sel_iso_xgb,
                                          sel_ag = c("All"),
                                          sel_measure = c("ILI_CASES"),
                                          minYear = 2010,
                                          maxYear = 2018
                                          )


#### XGBooost ####

# example of datasets
USA_complex1 <- gbm_complex(fluWHO.incidence, "USA", 10,1) # dataframe for 1-week ahead 
USA_complex2 <- gbm_complex(fluWHO.incidence, "USA", 10,2) # dataframe for 2-week ahead
USA_complex3 <- gbm_complex(fluWHO.incidence, "USA", 10,3) # dataframe for 3-week ahead
USA_complex4 <- gbm_complex(fluWHO.incidence, "USA", 10,4) # dataframe for 4-week ahead


country_xgb <- c()
for (i in 1:length(country_xgboost)){
  index <- grep(country_xgboost[i],countryISO$ISO3)
  tmp <- countryISO[index,]
  country_xgb <- rbind(country_xgb, tmp)
}

#' Convert dataframe into matrix
xgboost_dat <- function(flu_data_complex, start_year, end_year){
  require(dplyr)
  
  flu_data_complex$month <- as.factor(flu_data_complex$month)
  flu_data_complex$season <- as.factor(flu_data_complex$season)
  
  start_year <- as.character(start_year)
  end_year <- as.character(end_year)
  
  start_index <- grep(start_year,rownames(flu_data_complex))[1]
  end <- grep(end_year,rownames(flu_data_complex))
  end_index <- end[length(end)]
  
  # XGBoost requires the classes to be in a numeric format, starting with 0.
  dat_labels <- as.numeric(flu_data_complex$Y_week0)-1
  
  dat <- model.matrix(~.+0, contrasts.arg = lapply(flu_data_complex[,4:5], contrasts, contrasts=FALSE),
                      data = flu_data_complex[,-1])
  
  if (class(dat_labels) != "numeric"){
    dat_labels <- as.numeric(as.factor(dat_labels))
  }
  new_dat <- dat[start_index:end_index,]
  new_dat_labels <- dat_labels[start_index:end_index]
  xgb_dat <- xgb.DMatrix(data = new_dat,label = new_dat_labels)
  xgb_dat
}


#' XGBoost model
#' train_num_start: to calculate start year of trainign set. The origin is 2010, 0 represents 2010, 1 represents 2010+1=2011
#' , and on 
#' train_num_end: to calculate the end year of traingin set. End year = start year + train_num_end 
#' = 2010 + train_num_start +train_num_end
xgboost.model.pred <- function(flu_data, country, num_category,
                               train_num_start, train_num_end, nWeek_ahead){
  # set up dataset for xgboost
  flu_data_complex <- gbm_complex(flu_data, country, 10,nWeek_ahead)

  year_start <- min(as.numeric(substr(rownames(flu_data_complex),0,4)))
  year_end <- max(as.numeric(substr(rownames(flu_data_complex),0,4)))
  start_year_tr <- year_start + train_num_start
  end_year_tr <-  start_year_tr + train_num_end
  start_year_ts <- end_year_tr + 1
  
  if ((start_year_ts == 2017 && year_end == 2018) == TRUE ){
    end_year_ts <- year_end
  }else{
    end_year_ts <- start_year_ts
  }
  
  xgb_tr <- xgboost_dat(flu_data_complex, start_year_tr, end_year_tr)
  xgb_ts <- xgboost_dat(flu_data_complex, start_year_ts, end_year_ts)
  
  # train the xgboost model
  params.train <- list(booster = "gbtree", objective = "multi:softprob", gamma=0, num_class = 10,
                       subsample=1, colsample_bytree=1,eval_metric = "mlogloss")
  watchlist <- list(train = xgb_tr, test = xgb_ts)
  xgb_model <- xgb.train(params = params.train, data = xgb_tr, nrounds = 100, 
                         watchlist = watchlist,verbose = 2, print_every_n = 10,
                         early_stopping_round = 20)
  xgb_pred <- predict(xgb_model, newdata = xgb_ts)
  start_year_ts_index <- grep(start_year_ts,rownames(flu_data_complex))[1]
  end_ts <- grep(end_year_ts,rownames(flu_data_complex))
  end_year_ts_index <- end_ts[length(end_ts)]
  xgb_val_out <- matrix(xgb_pred, nrow = 10, ncol = length(xgb_pred)/10) %>% 
    t() %>%
    data.frame() %>%
    mutate(max = max.col(., ties.method = "last"), 
           category = flu_data_complex$Y_week0[start_year_ts_index:end_year_ts_index])
  
  pred_timeseries <- rownames(flu_data_complex)[start_year_ts_index:end_year_ts_index] %>% 
    cbind(xgb_val_out[,(ncol(xgb_val_out)-1):ncol(xgb_val_out)]) %>%
    data.frame()
  colnames(pred_timeseries) <- c("week_time", "Prediction", "Observation")
  pred_timeseries$Observation <- as.numeric(pred_timeseries$Observation)
  pred_timeseries$Prediction <- as.numeric(pred_timeseries$Prediction)
  for (i in 1:nrow(pred_timeseries)){
    if (pred_timeseries[i,2]==pred_timeseries[i,3]){
      pred_timeseries$Accurate[i] <- 1
    }
    else{
      pred_timeseries$Accurate[i] <- 0
    }
  }
  pred_timeseries
}

#' Function of caculating the accuracy metric
compare_accuracy <- function(country_list,train_num_start, train_num_end,nWeek_ahead){
  pred <- NULL
  for (i in 1:length(country_list)){
    individual_pred <- xgboost.model.pred(fluWHO.incidence,country_list[i],10,
                                          train_num_start, train_num_end,nWeek_ahead)
    individual_pred <- cbind(rep(country_list[i], nrow(individual_pred)),individual_pred)
    pred <- rbind(pred,individual_pred)
    
  }
  
  pred <- as.data.frame(pred)
  colnames(pred) <- c("Country","week_time","Observation","Prediction","Accurate")
  pred
}


# one-week ahead forecast
# 2010-2104 training, 2015 test
compare_pred15 <- compare_accuracy(sel_iso_xgb, 0, 4 ,1)

# 2011-2015 training, 2016 test
compare_pred16 <- compare_accuracy(country_xgb$ISO3, 1, 4, 1)

# 2012-2016 training. 2017 test
compare_pred17 <- compare_accuracy(country_xgb$ISO3, 2, 4, 1)

# 2010-2016 traing, 2017 test
compare1016_pred17 <- compare_accuracy(country_xgb$ISO3, 0, 6, 1)

# two-week ahead forecast
# 2010-2104 training, 2015 test
compare_pred15Two <- compare_accuracy(country_xgb$ISO3, 0, 4 ,2)

# 2011-2015 training, 2016 test
compare_pred16Two <- compare_accuracy(country_xgb$ISO3, 1, 4,2)

# 2012-2016 training. 2017 test
compare_pred17Two <- compare_accuracy(country_xgb$ISO3, 2, 4, 2)

# 2010-2016 traing, 2017 test
compare1016_pred17Two <- compare_accuracy(country_xgb$ISO3, 0, 6, 2)

# three week ahead
# 2010-2104 training, 2015 test
compare_pred15Three <- compare_accuracy(country_xgb$ISO3, 0, 4 ,3)

# 2011-2015 training, 2016 test
compare_pred16Three <- compare_accuracy(country_xgb$ISO3, 1, 4, 3)

# 2012-2016 training. 2017 test
compare_pred17Three <- compare_accuracy(country_idd$ISO3, 2, 4, 3)

# 2010-2016 traing, 2017 test
compare1016_pred17Three <- compare_accuracy(country_xgb$ISO3, 0, 6, 3)

# four week ahead
# 2010-2104 training, 2015 test
compare_pred15Four <- compare_accuracy(country_xgb$ISO3, 0, 4 ,4)

# 2011-2015 training, 2016 test
compare_pred16Four <- compare_accuracy(country_xgb$ISO3, 1, 4, 4)

# 2012-2016 training. 2017 test
compare_pred17Four <- compare_accuracy(country_xgb$ISO3, 2, 4, 4)

# 2010-2016 traing, 2017 test
compare1016_pred17Four <- compare_accuracy(country_xgb$ISO3, 0, 6, 4)


#### Heat plot for xgboost model ####
freq_table <- function(prediction, row_col){
  accuracy <- prediction[,2:3]
  for (i in 1:nrow(accuracy)){
    if (accuracy[i,1]==accuracy[i,2]){
      accuracy$accurate[i] <- 1
    }
    else{
      accuracy$accurate[i] <- 0
    }
  }
  per <- gtools::permutations(row_col,2,repeats.allowed = TRUE)
  freq <- cbind(per, rep(0, nrow(per)))
  freq <- as.data.frame(freq)
  colnames(freq) <- c('observation_category', 'forecast_category', 'frequency')
  forecast_freq <- ftable(accuracy[1:2], row.vars = 2:1)
  forecast_freq <- as.data.frame(forecast_freq)
  colnames(forecast_freq) <- c('observation_category', 'forecast_category', 'frequency')
  # calculate frequency 
  j <- 1
  while(j <= nrow(per)){
    for (i in 1:nrow(forecast_freq)){
      if (forecast_freq$observation_category[i] %in% freq$observation_category[j] && 
          forecast_freq$forecast_category[i] %in% freq$forecast_category[j] == TRUE){
        freq$frequency[j] <- forecast_freq$frequency[i]
      }
    }
    j <- j+1
  }
  freq
  
  squared_freq <- matrix(nrow = row_col, ncol = row_col)
  n <- 1
  while(n <= nrow(freq)){
    for (j in 1:ncol(squared_freq)){
      for (i in 1:nrow(squared_freq)){
        squared_freq[i,j] <- freq$frequency[n]
        n <- n+1
      }
    }
  }
  squared_freq <- as.data.frame(squared_freq)
  colnames(squared_freq) <- c('observed_1','observed_2','observed_3','observed_4','observed_5','observed_6',
                              'observed_7','observed_8','observed_9','observed_10')
  rownames(squared_freq) <- c('predicted_1','predicted_2','predicted_3','predicted_4','predicted_5',
                              'predicted_6','predicted_7','predicted_8','predicted_9','predicted_10')
  squared_freq
}

heat_plot <- function(frequencyTable, countryName){
  frequencyMatrix <- as.matrix(frequencyTable)
  colnames(frequencyMatrix) <- c(1:10)
  rownames(frequencyMatrix) <- c(1:10)
  # my.at will be changed according to the number of data points
  my.at <- c(0,1,5,10,15,20,25,30,40,50)
  my.brks <- seq(0, max(frequencyMatrix, na.rm = TRUE), length.out = length(my.at))
  blues <- brewer.pal(9, "Blues")
  reds <-  brewer.pal(9, "Reds")
  mapTheme <- rasterTheme(region= c(blues,reds))
  myColorkey <- list(at=my.brks, labels=list(at=my.brks, labels=my.at), space="right")
  myPanel <- function(x, y, z,...) {
    panel.levelplot(x,y,z,...)
    for (i in 1:nrow(frequencyMatrix)){
      for (j in 1:ncol(frequencyMatrix)){
        panel.text(x=i, y=j, frequencyMatrix[cbind(j,i)])
      }
    }
  }
  country <- countryName 
  levelplot(t(frequencyMatrix), xlab = 'Obeserved', ylab = 'Forecast',
            panel = myPanel, par.settings=mapTheme, at=my.at, colorkey=myColorkey, margin=F)
}

#' Print out heat chart 
# 2010-2014 training, 2015 test
for (i in 1:length(sel_iso_xgb)){
  #pdf(paste0(sel_iso_xgb[i],".pdf"))
  forecast_result <- xgboost.model.pred(fluWHO.incidence, sel_iso_xgb[i], 10, 0, 4, 1)
  frequency_table <- freq_table(forecast_result, 10)
  print(heat_plot(frequency_table, sel_iso_xgb[i]))
  #dev.off()
}

# 2011-2015 train, 2016 test
for (i in 1:length(sel_iso_xgb)){
  pdf(paste0(sel_iso_xgb[i],".pdf"))
  forecast_result <- xgboost.model.pred(fluWHO.incidence, sel_iso_xgb[i], 10, 1, 4, 1)
  frequency_table <- freq_table(forecast_result, 10)
  print(heat_plot(frequency_table, sel_iso_xgb[i]))
  dev.off()
}

# 2012-2016 train, 2017 test
for (i in 1:length(sel_iso_xgb)){
  pdf(paste0(sel_iso_xgb[i],".pdf"))
  forecast_result <-  xgboost.model.pred(fluWHO.incidence, sel_iso_xgb[i], 10, 2, 4, 1)
  frequency_table <- freq_table(forecast_result, 10)
  print(heat_plot(frequency_table, sel_iso_xgb[i]))
  dev.off()
}


#### historical average model ####
#' Historical avarage model to generate category predictions
hist_dataform <- function(flu_data){
  
  year <- unique(substr(rownames(flu_data),0,4))
  year_min <- min(as.numeric(year))
  year_max <- max(as.numeric(year))
  
  # if(year_max == 2018){
  # index_2018 <- grep("2018", rownames(flu_data))
  # flu_data <- flu_data[-index_2018,]
  # }
  
  yr <- seq(year_min,year_max,1)
  week <- c(seq(27,52,1),seq(1,26,1))
  hist <- matrix(nrow = 52, ncol = length(yr))
  for (i in 1:length(yr)){
    year_week <- paste0(yr[i], "-", week)
    for (j in 1:length(year_week))
      if (year_week[j] %in% rownames(flu_data) == TRUE){
        row_index <- grep(year_week[j],rownames(flu_data))
        hist[j,i] <- flu_data$Y_week0[row_index]
      }else{
        hist[j,i] <- NA
      }
  }
  
  colnames(hist) <- year
  rownames(hist) <- paste0("week", week)
  
  #na_index <- c()
  # for (i in 1:nrow(hist)){
  # if(length(which(is.na(hist[i,]))) >= 3){
  # tmp <- i
  # na_index <- append(na_index, tmp)
  # }
  # }
  # hist <- hist[-na_index,]
  # hist <- as.data.frame(hist)
  hist
}

hist_average <- function(flu_data, country,num_category, numWeek_ahead){
  flu_data_complex <- gbm_complex(flu_data, country, num_category,numWeek_ahead)
  flu_data_complex_four <- gbm_complex(flu_data, country, num_category,4)
  # hist_dataset <- hist_dataform(flu_data_complex)
  flu_data_complex <- cbind(substr(rownames(flu_data_complex), 0,4), 
                            substr(rownames(flu_data_complex),6,7),
                            flu_data_complex[,1:3]) %>% as.data.frame()
  colnames(flu_data_complex) <- c("Year","Week","Y_week0","week_1","week_2")
  flu_data_complex$Week <- as.numeric(flu_data_complex$Week)
  
  flu_data_complex_four <- cbind(substr(rownames(flu_data_complex), 0,4), 
                                 substr(rownames(flu_data_complex),6,7),
                                 flu_data_complex[,1:3]) %>% as.data.frame()
  colnames(flu_data_complex_four) <- c("Year","Week","Y_week0","week_4","week_5")
  flu_data_complex_four$Week <- as.numeric(flu_data_complex_four$Week)

  pred <- matrix(NA,nrow = nrow(flu_data_complex), ncol = numWeek_ahead)
  
  if (numWeek_ahead == 1){
    for (i in 1:nrow(flu_data_complex)){
      yr <- flu_data_complex$Year[i]
      week <- flu_data_complex$Week[i]-1
      obsTem <- flu_data_complex[which(flu_data_complex$Week==week),]
      obs <- obsTem$Y_week0[which(obsTem$Year != yr)]
      pred[i,] <- which.max(tabulate(obs))
    }
    
    pred <- cbind(rownames(flu_data_complex),pred,flu_data_complex$Y_week0) %>% 
      as.data.frame()
    colnames(pred) <- c("Week_time","OneWeek_ahead", "Observation")

    # accuracy
    for (i in 1:nrow(pred)){
      if (is.na(pred[i,2])==TRUE || is.na(pred[i,3])==TRUE){
        pred$accurate[i] <- NA
      }else{
        if (pred[i,2]==pred[i,3]){
          pred$Accurate[i] <- 1
        }else{
          pred$Accurate[i] <- 0
        }
      }
    }
  }
  
  if (numWeek_ahead == 2){
    for (i in 1:nrow(flu_data_complex)){
      yr <- flu_data_complex$Year[i]
      week <- flu_data_complex$Week[i]-2
      obsTem <- flu_data_complex[which(flu_data_complex$Week==week),]
      obs <- obsTem$Y_week0[which(obsTem$Year != yr)]
      pred[i,] <- which.max(tabulate(obs))
    }
    pred <- cbind(rownames(flu_data_complex),pred,flu_data_complex$Y_week0) %>% 
      as.data.frame()
    colnames(pred) <- c("Week_time","OneWeek_ahead", "TwoWeek_ahead","Observation")
    
    # accuracy
    for (i in 1:nrow(pred)){
      if (is.na(pred[i,3])==TRUE || is.na(pred[i,4])==TRUE){
        pred$accurate[i] <- NA
      }else{
        if (pred[i,3]==pred[i,4]){
          pred$Accurate[i] <- 1
        }else{
          pred$Accurate[i] <- 0
        }
      }
    }
  }
  
  if (numWeek_ahead == 3){
    for (i in 1:nrow(flu_data_complex)){
      yr <- flu_data_complex$Year[i]
      week <- flu_data_complex$Week[i]-3
      obsTem <- flu_data_complex[which(flu_data_complex$Week==week),]
      obs <- obsTem$Y_week0[which(obsTem$Year != yr)]
      pred[i,] <- which.max(tabulate(obs))
    }
    pred <- cbind(rownames(flu_data_complex),pred,flu_data_complex$Y_week0) %>% 
      as.data.frame()
    colnames(pred) <- c("Week_time","OneWeek_ahead", "TwoWeek_ahead","ThreeWeek_ahead","Observation")
    
    # accuracy
    for (i in 1:nrow(pred)){
      if (is.na(pred[i,4])==TRUE || is.na(pred[i,5])==TRUE){
        pred$accurate[i] <- NA
      }else{
        if (pred[i,4]==pred[i,5]){
          pred$Accurate[i] <- 1
        }else{
          pred$Accurate[i] <- 0
        }
      }
    }
  }
  
  if (numWeek_ahead == 4){
    for (i in 1:nrow(flu_data_complex)){
      yr <- flu_data_complex$Year[i]
      week <- flu_data_complex$Week[i]-4
      obsTem <- flu_data_complex[which(flu_data_complex$Week==week),]
      obs <- obsTem$Y_week0[which(obsTem$Year != yr)]
      pred[i,] <- which.max(tabulate(obs))
    }
    pred <- cbind(rownames(flu_data_complex),pred,flu_data_complex$Y_week0) %>% 
      as.data.frame()
    colnames(pred) <- c("Week_time","OneWeek_ahead", "TwoWeek_ahead","ThreeWeek_ahead","FourWeek_ahead","Observation")
    
    # accuracy
    for (i in 1:nrow(pred)){
      if (is.na(pred[i,5])==TRUE || is.na(pred[i,6])==TRUE){
        pred$accurate[i] <- NA
      }else{
        if (pred[i,5]==pred[i,6]){
          pred$Accurate[i] <- 1
        }else{
          pred$Accurate[i] <- 0
        }
      }
    }
  }
  pred
}

# example 
USA_hist_pred_one <- hist_average(fluWHO.incidence,"USA",10,1)
USA_hist_pred_two <- hist_average(fluWHO.incidence,"USA",10,2)
USA_hist_pred_three <- hist_average(fluWHO.incidence,"USA",10,3)
USA_hist_pred_four <- hist_average(fluWHO.incidence,"USA",10,4)


# compare 1,2,3,4-week ahead forecast accuracy
compare_accuracy_hist <- function(flu_data,country,num_category,numWeek_ahead){
  pred <- NULL
  
  for (i in 1:length(country)){
    individual_pred <- hist_average(flu_data,country[i],num_category,numWeek_ahead)
    individual_pred <- cbind(rep(country[i], nrow(individual_pred)),individual_pred)
    pred <- rbind(pred,individual_pred)
  }
  pred <- as.data.frame(pred)
  colnames(pred) <- c("Country",colnames(pred)[2:ncol(pred)])
  pred
}

# 8842/13603 = 0.65
oneWeek_ahead_totalAccuracy_hist <- compare_accuracy_hist(fluWHO.incidence,sel_iso_xgb,10,1)
length(which(oneWeek_ahead_totalAccuracy_hist$Accurate==1))/nrow(oneWeek_ahead_totalAccuracy_hist)

# 0.64
twoWeek_ahead_totalAccuracy_hist <- compare_accuracy_hist(fluWHO.incidence,sel_iso_xgb,10,2)
length(which(twoWeek_ahead_totalAccuracy_hist$Accurate==1))/nrow(twoWeek_ahead_totalAccuracy_hist)

# 0.63
threeWeek_ahead_totalAccuracy_hist <- compare_accuracy_hist(fluWHO.incidence,sel_iso_xgb,10,3)
length(which(threeWeek_ahead_totalAccuracy_hist$Accurate==1))/nrow(threeWeek_ahead_totalAccuracy_hist)

# 0.62
fourWeek_ahead_totalAccuracy_hist <- compare_accuracy_hist(fluWHO.incidence,sel_iso_xgb,10,4)
length(which(fourWeek_ahead_totalAccuracy_hist$Accurate==1))/nrow(fourWeek_ahead_totalAccuracy_hist)

# plot the accuracy drop off as number of week-ahead increases
compareAccuracy_total_hist <- cbind(c(1:4),c(0.65,0.64,0.62,0.61)) %>% 
  as.data.frame()
colnames(compareAccuracy_total_hist) <- c("nWeek_ahead","percentage")

#### repeat model ####
repeat_model <- function(flu_data,country, num_category, numWeek_ahead){
  require(dplyr)
  
  flu_data_complex <- gbm_complex(flu_data, country, num_category,numWeek_ahead)                                  
  # prediction of the week is the same as the last week
  if(numWeek_ahead == 1){
    pred <- c()
    for (i in 1:nrow(flu_data_complex)){
      tmp <- flu_data_complex$week_1[i]
      pred <- append(pred,tmp)
    }
    pred <- cbind(rownames(flu_data_complex),pred, flu_data_complex$Y_week0) %>% 
      as.data.frame()
    colnames(pred) <- c("Week_time","Prediction","Observation")
    
    # accuracy
    for (i in 1:nrow(pred)){
      if (is.na(pred[i,2])==TRUE || is.na(pred[i,3])==TRUE){
        pred$accurate[i] <- NA
      }else{
        if (pred[i,2]==pred[i,3]){
          pred$Accurate[i] <- 1
        }else{
          pred$Accurate[i] <- 0
        }
      }
    }
  }
  
  if (numWeek_ahead == 2){
    pred <- matrix(NA, nrow = nrow(flu_data_complex),ncol = 2)
    for (i in 1:nrow(pred)){
      pred[i,] <- flu_data_complex$week_2[i]
    }
    
    pred <- cbind(rownames(flu_data_complex),pred,flu_data_complex$Y_week0) %>% 
      as.data.frame()
    colnames(pred) <- c('Week_time','OneWeek_ahead','TwoWeek_ahead', "Observation")
    
    # accuracy
    for (i in 1:nrow(pred)){
      if (is.na(pred[i,3])==TRUE || is.na(pred[i,4])==TRUE){
        pred$accurate[i] <- NA
      }else{
        if (pred[i,3]==pred[i,4]){
          pred$Accurate[i] <- 1
        }else{
          pred$Accurate[i] <- 0
        }
      }
    }
  }
  
  if (numWeek_ahead == 3){
    pred <- matrix(NA, nrow = nrow(flu_data_complex),ncol = 3)
    for (i in 1:nrow(pred)){
      pred[i,] <- flu_data_complex$week_3[i]
    }
    
    pred <- cbind(rownames(flu_data_complex),pred,flu_data_complex$Y_week0) %>% 
      as.data.frame()
    colnames(pred) <- c("Week_time","OneWeek_ahead","TwoWeek_ahead","ThreeWeek_ahead", "Observation")

    # accuracy
    for (i in 1:nrow(pred)){
      if (is.na(pred[i,4])==TRUE || is.na(pred[i,5])==TRUE){
        pred$accurate[i] <- NA
      }else{
        if (pred[i,4]==pred[i,5]){
          pred$Accurate[i] <- 1
        }else{
          pred$Accurate[i] <- 0
        }
      }
    }
  }
  
  if (numWeek_ahead == 4){
    pred <- matrix(NA, nrow = nrow(flu_data_complex),ncol = 4)
    for (i in 1:nrow(pred)){
      pred[i,] <- flu_data_complex$week_4[i]
    }
    
    pred <- cbind(rownames(flu_data_complex),pred,flu_data_complex$Y_week0) %>% 
      as.data.frame()
    colnames(pred) <- c("Week_time","OneWeek_ahead","TwoWeek_ahead","ThreeWeek_ahead","FourWeek_ahead", "Observation")
    # accuracy
    for (i in 1:nrow(pred)){
      if (is.na(pred[i,5])==TRUE || is.na(pred[i,6])==TRUE){
        pred$accurate[i] <- NA
      }else{
        if (pred[i,5]==pred[i,6]){
          pred$Accurate[i] <- 1
        }else{
          pred$Accurate[i] <- 0
        }
      }
    }
  }
  
  pred
}

# example 
USA_repeat_one <- repeat_model(fluWHO.incidence,'USA',10,1)
USA_repeat_two <- repeat_model(fluWHO.incidence,'USA',10,2)
USA_repeat_three <- repeat_model(fluWHO.incidence,'USA',10,3)
USA_repeat_four <- repeat_model(fluWHO.incidence,'USA',10,4)

intersect(intersect(intersect(USA_repeat_one$Week_time,USA_repeat_two$Week_time),USA_repeat_three$Week_time),
          USA_repeat_four$Week_time)


compare_accuracy_repeat <- function(flu_data,country,num_category,numWeek_ahead){
  pred <- NULL
  
  for (i in 1:length(country)){
    individual_pred <- repeat_model(flu_data,country[i],num_category, numWeek_ahead)
    individual_pred <- cbind(rep(country[i], nrow(individual_pred)),individual_pred)
    pred <- rbind(pred,individual_pred)
  }
  pred <- as.data.frame(pred)
  colnames(pred) <- c("Country",colnames(pred)[2:ncol(pred)])
  pred
  
}

# get the total accuract for 1,2,3,4-weeks ahead
oneWeek_ahead_totalAccuracy <- compare_accuracy_repeat(fluWHO.incidence,sel_iso_xgb,10,1)
length(which(oneWeek_ahead_totalAccuracy$Accurate==1))  # 9641 out of 12081
length(which(oneWeek_ahead_totalAccuracy$Accurate==1))/nrow(oneWeek_ahead_totalAccuracy) # 0.79803

twoWeek_ahead_totalAccuracy <- compare_accuracy_repeat(fluWHO.incidence,sel_iso_xgb,10,2)
length(which(twoWeek_ahead_totalAccuracy$Accurate==1)) # 8896
length(which(twoWeek_ahead_totalAccuracy$Accurate==1))/nrow(twoWeek_ahead_totalAccuracy) # 0.7440616

threeWeek_ahead_totalAccuracy <- compare_accuracy_repeat(fluWHO.incidence,sel_iso_xgb,10,3)
length(which(threeWeek_ahead_totalAccuracy$Accurate==1)) # 8332
length(which(threeWeek_ahead_totalAccuracy$Accurate==1))/nrow(threeWeek_ahead_totalAccuracy) #0.7043111

fourWeek_ahead_totalAccuracy <- compare_accuracy_repeat(fluWHO.incidence,sel_iso_xgb,10,4)
length(which(fourWeek_ahead_totalAccuracy$Accurate==1)) # 7880
length(which(fourWeek_ahead_totalAccuracy$Accurate==1))/nrow(fourWeek_ahead_totalAccuracy) # 0.67

# plot the accuracy drop off as number of week-ahead increases
compareAccuracy_total <- cbind(c(1:4),c(10181, 9401, 8415,7996),c(0.75,0.69,0.62,0.59))
compareAccuracy_total <- as.data.frame(compareAccuracy_total)
colnames(compareAccuracy_total) <- c("nWeek_ahead","accuracy","percentage")
compareAccuracy_total$nWeek_ahead <- as.numeric(compareAccuracy_total$nWeek_ahead)

compareAccuracy_total$nWeek_ahead <- factor(compareAccuracy_total$nWeek_ahead,
                                            levels = c("one-week ahead","two-week ahead",
                                                       "three-week ahead","four-week ahead"))
compareAccuracy_total$accuracy <- as.numeric(as.character(compareAccuracy_total$accuracy))
compareAccuracy_total$percentage <- as.numeric(as.character(compareAccuracy_total$percentage))
compareAccuracy_total$percentage <- factor(compareAccuracy_total$percentage,levels = c(0,0.2,0.4,0.6,0.8,1))