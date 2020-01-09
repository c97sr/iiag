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

#' Function that helps identify correponding country's full name according to ISO3 code
iso3_country <- function(iso3_code){
  country <- countryISO$Country[which(iso3_code == countryISO$ISO3)]
  
  return(country)
}


#' extract_incidence is the function in package idd which can'y be loaded
#' therefore copy the code 
extract.incidence.idd <- function(flu_data,
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

#' check the data availablity in each year
duration <- function(country,maxYear, minYear){
  year_time <- c(minYear:maxYear)
  
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

#' Keep data frames used for 1,2,3,4-week ahead forecast are in the same size
adjust.data.size <- function(flu_data,country,category,numWeek_ahead){
  complex1 <- gbm_complex(flu_data,country,category,1)
  complex2 <- gbm_complex(flu_data,country,category,2)
  complex3 <- gbm_complex(flu_data,country,category,3)
  complex4 <- gbm_complex(flu_data,country,category,4)
  week <- intersect(rownames(complex1),intersect(rownames(complex2),
                                                 intersect(rownames(complex3),rownames(complex4))))
  
  complex <- gbm_complex(flu_data,country,category,numWeek_ahead)
  index <- c()
  for (i in 1:nrow(complex)){
    if(rownames(complex)[i] %in% week == FALSE){
      tmp <- i
      index <- append(index, tmp)
    }
  }
  if (length(index) == 0){
    complex <- complex
  }else{
    complex <- complex[-(index),]
  }
  
  complex
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
#' train_num_start: calculates start year of trainign set. The origin is 2010, 0 represents 2010, 1 represents 2010+1=2011
#' train_num_end: calculates the end year of traingin set. End year = start year + train_num_end 
#' = 2010 + train_num_start +train_num_end
xgboost.model.pred <- function(flu_data, country, num_category,
                               train_num_start, train_num_end, nWeek_ahead){
  # set up dataset for xgboost
  flu_data_complex <- adjust.data.size(flu_data, country, num_category, nWeek_ahead)
  
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

#' Function that gives individual country forecast result and accuracy score
compare_accuracy_indi <- function(individual_country, flu_data, num_category,train_num_start, 
                                  train_num_end,nWeek_ahead){
  country_list <- individual_country
  individual_pred <- xgboost.model.pred(fluWHO.incidence,country_list,num_category,
                                        train_num_start, train_num_end,nWeek_ahead)
  individual_pred <- cbind(rep(individual_country, nrow(individual_pred)),individual_pred)
  individual_pred <- as.data.frame(individual_pred)
  colnames(individual_pred) <- c("Country","week_time","Observation","Prediction","Accurate")
  
  score <- length(which(individual_pred$Accurate == 1))/nrow(individual_pred)
  
  result <- NULL
  result$individual_pred <- individual_pred
  result$score <- score
  
  return(result)
}


#' Function of caculating the accuracy metric of xgboost model
compare_accuracy <- function(country_list,flu_data,num_category, train_num_start, train_num_end,nWeek_ahead){
  pred <- NULL
  for (i in 1:length(country_list)){
    individual_pred <- xgboost.model.pred(flu_data,country_list[i],num_category,
                                          train_num_start, train_num_end,nWeek_ahead)
    individual_pred <- cbind(rep(country_list[i], nrow(individual_pred)),individual_pred)
    pred <- rbind(pred,individual_pred)

  }
  pred <- as.data.frame(pred)
  colnames(pred) <- c("Country","week_time","Observation","Prediction","Accurate")
  
  score <- length(which(pred$Accurate == 1))/nrow(pred)
  
  result <- NULL
  result$pred <- pred
  result$score <- score
  
  return(result)

}


#' shows the output of xgboost model 
xgboost.model.train <- function(flu_data, country, num_category,
                                train_num_start, train_num_end, nWeek_ahead){
  # set up dataset for xgboost
  flu_data_complex <- gbm_complex(flu_data, country, num_category, nWeek_ahead)
  
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
  xgb_model
}

#' shows the prediction results of xgboost model
xgboost.model.pred.output <- function(flu_data_complex, start_year_ts,
                                      end_year_ts, xgb_pred){
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
  pred_timeseries
}

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
  my.at <- c(0,10,30,50,100,150,200,250,300,350,550)
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

#### historical average model ####
#' Historical avarage model to generate category predictions
# hist_dataform <- function(flu_data){

# year <- unique(substr(rownames(flu_data),0,4))
# year_min <- min(as.numeric(year))
# year_max <- max(as.numeric(year))

# if(year_max == 2018){
# index_2018 <- grep("2018", rownames(flu_data))
# flu_data <- flu_data[-index_2018,]
# }

# yr <- seq(year_min,year_max,1)
# week <- c(seq(27,52,1),seq(1,26,1))
# hist <- matrix(nrow = 52, ncol = length(yr))
# for (i in 1:length(yr)){
# year_week <- paste0(yr[i], "-", week)
# for (j in 1:length(year_week))
# if (year_week[j] %in% rownames(flu_data) == TRUE){
# row_index <- grep(year_week[j],rownames(flu_data))
# hist[j,i] <- flu_data$Y_week0[row_index]
# }else{
# hist[j,i] <- NA
# }
# }

# colnames(hist) <- year
# rownames(hist) <- paste0("week", week)

#na_index <- c()
# for (i in 1:nrow(hist)){
# if(length(which(is.na(hist[i,]))) >= 3){
# tmp <- i
# na_index <- append(na_index, tmp)
# }
# }
# hist <- hist[-na_index,]
# hist <- as.data.frame(hist)
# hist
# }

hist_average <- function(flu_data, country,num_category, numWeek_ahead){
  
  flu_data_complex <- adjust.data.size(flu_data,country,num_category,numWeek_ahead)
  
  flu_data_complex <- cbind(substr(rownames(flu_data_complex), 0,4), 
                            substr(rownames(flu_data_complex),6,7),
                            flu_data_complex[,1:3]) %>% as.data.frame()
  colnames(flu_data_complex) <- c("Year","Week","Y_week0","week_1","week_2")
  flu_data_complex$Week <- as.numeric(flu_data_complex$Week)
  
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
  
  score <- length(which(pred$Accurate == 1))/nrow(pred)
  
  result <- NULL
  result$pred <- pred
  result$score <- score
  
  return(result)
}

#### repeat model ####

#'
repeat_model <- function(flu_data,country, num_category, numWeek_ahead){
  require(dplyr)
  
  flu_data_complex <- adjust.data.size(flu_data,country,num_category,numWeek_ahead)                                
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

#' calculate the accuracy of repeat model
compare_accuracy_repeat <- function(flu_data,country,num_category,numWeek_ahead){
  pred <- NULL
  
  for (i in 1:length(country)){
    individual_pred <- repeat_model(flu_data,country[i],num_category, numWeek_ahead)
    individual_pred <- cbind(rep(country[i], nrow(individual_pred)),individual_pred)
    pred <- rbind(pred,individual_pred)
  }
  pred <- as.data.frame(pred)
  colnames(pred) <- c("Country",colnames(pred)[2:ncol(pred)])
  
  score <- length(which(pred$Accurate == 1))/nrow(pred)
  
  result <- NULL
  result$pred <- pred
  result$score <- score
  
  return(result)
  
}

#### functions for plotting ####

#' discrete data against raw data
#' Plot raw incidence data and categorical data in the same graph by dual axis is not clear for visualizing
#' So I will plot them into two seperate graphs but one is above another one
#' data: incidenc dataset; country: ISO3 code; countryName: full country name
rawData_TS_plot <- function(data, country,countryName){
  require(ggplot2)
  require(grid)
  
  count <- data[, which(colnames(data)==country)]
  week_time <- rownames(data)[-which(is.na(count))]
  count <- count[-(which(is.na(count)))]
  category <- cut_interval(count, 10)
  raw_data <- cbind(week_time, count, category)
  
  raw_data <- as.data.frame(raw_data)
  colnames(raw_data) <- c('week_time','count','category')
  
  # convert discrete week time to contious time series 
  raw_data$weekTS <- gsub("-", "-W", raw_data$week_time, fixed = TRUE)
  raw_data$weekTS <- week2date(raw_data$weekTS)
  raw_data$count <- as.numeric(as.character(raw_data$count))
  raw_data$category <- as.numeric(as.character(raw_data$category))
  
  # first plot the time series of numeric incidence
  grobTS <- grobTree(textGrob(countryName, x=0.85,  y=0.9, hjust=0,
                              gp=gpar(col="black", fontsize=14, fontface="bold")))
  p1 <- ggplot(raw_data, aes(x = weekTS))
  p1 <- p1 + geom_line(aes(y = count),colour = "salmon",size=1)
  p1 <- p1 + scale_y_continuous(expand = c(0, 0))
  # p1 <- p1 + scale_x_continuous(breaks = c(2010:2018))
  # p1 <- p1 + scale_colour_manual(values = "salmon")
  p1 <- p1 + labs(y = "Incidence",
                  x = "Time",
                  # colour = "Obsevation",
                  title = "Time series of numeric observations of incidence")
  p1 <- p1 + theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    annotation_custom(grobTS)
  
  # secondly plot the time series of categorical incidence
  p2 <- ggplot(raw_data, aes(x = weekTS))
  p2 <- p2 + geom_line(aes(y = category),colour = "darkslategray2",size=1)
  # p2 <- p2 + scale_x_continuous(breaks = c(2010:2018))
  p2 <- p2 + scale_y_continuous(breaks = c(1:10),expand = c(0, 0))
  # p2 <- p2 + scale_colour_manual(values = "darkslategray2")
  p2 <- p2 + labs(y = "Incidence",
                  x = "Time",
                  # colour = "Obsevation",
                  title = "Time series of categorical observations of incidence")
  p2 <- p2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  grid.newpage()
  grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
  
}

#' time series plot for prediction
#' pred: prediction results from xgboost.model.pred function
predTS_plot <- function(pred){
  # change week into format used for aweek package
  pred$weekTS <- gsub("-", "-W", pred$week_time, fixed = TRUE)
  
  # the aweek package is used to convert from epidemiological weeks to dates
  pred$weekTS <- week2date(pred$weekTS)
  
  # convert numbers from factors to numerics
  pred$Observation <- as.numeric(as.character(pred$Observation))
  pred$Prediction <- as.numeric(as.character(pred$Prediction))
  
  p <- ggplot(pred, aes(x = weekTS))
  
  # plot observational category 
  p <- p + geom_point(aes(y = Observation, color = "Observation"),size=3)
  p <- p + geom_line(aes(y = Observation, color = "Observation"),size=0.8)
  
  # plot predicted category
  p <- p + geom_point(aes(y = Prediction, color = "Prediction"),size=3)
  p <- p + geom_line(aes(y = Prediction, color = "Prediction"),size=0.8)
  
  # use segment to show the difference between prediction and obsercation more clearly
  # p <- p + geom_segment(aes(xend = weekTS, yend = Prediction), alpha = .2)
  
  # modifying colours and theme options
  p <- p + scale_colour_manual(values = c("darkslategray2", "salmon"))
  p <- p + scale_y_continuous(breaks = c(1:10))
  p <- p + labs(y = "Category",
                x = "Time",
                colour = "Category")
  p <- p + theme(legend.position = c(0.8,0.85))
  
  p
}


#' pick up European countries 
euro_countries <- function(countries){
  require(eurostat)
  data("eu_countries")
  euro <- c()
  for (i in 1:nrow(countries)){
    if (as.character(countries$Country[i]) %in% as.character(eu_countries$name)){
      tmp <- countries[i,]
      euro <- rbind(euro,tmp)
    }
  }
  euro
}

#' Check if countries names is same as names in "maps" package
map.country.name <- function(using_country_name){
  require(maps)
  
  world_map <- map_data("world")
  index <- which(using_country_name %in% unique(world_map$region)==FALSE)
  index
  
}

#' Tidy up accuracy scores of countries into a dataframe
accuracy_score <- function(forecast_result, country_code){
  dataframe <- matrix(NA, nrow = length(country_code), ncol = 2)
  
  for (i in 1:length(country_code)){
       dataframe[i,1] <- country_code[i]
        dataframe[i,2] <- forecast_result[[i*2]]
    }

  dataframe <- as.data.frame(dataframe)
  colnames(dataframe) <- c("Country","Accuracy_score")
  
  dataframe$Accuracy_score <- round(as.numeric(as.character(dataframe$Accuracy_score)),3)

  return(dataframe)
}

