
#' Make a note of when the report was generated.
Sys.time()

#' As always, remove all objects fromt the workspace before starting.
rm(list = ls(all = TRUE))


# library('idd') # WHO influenza-likeness-illness data; idd package is not available (for R version 3.6.1)
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
library("devtools") # to install package from github
# install_github("AppliedDataSciencePartners/xgboostExplainer")
# library("xgboostExplainer") # to interpret every single tree 
# library("pROC")


# load(file = "fluIliCountryData.rda")

source("C:/Users/hw3616/Desktop/Imperial/Project1_Forecasting/Project_Coding/iiag/forecasting_vivi/vivi_funcs.R")
source("C:/Users/hw3616/Desktop/Imperial/Project1_Forecasting/Project_Coding/iiag/forecasting_vivi/gbm_complex_funcs.R")


#' Load country list for WHO data.
countryISO <- read.csv("C:/Users/hw3616/Desktop/Imperial/Project1_Forecasting/Project_Coding/iiag/data_old/country_list_ISO.csv")

#' Load WHO FluID data set.
#' From 2010 week 1 to 2018 week 9 

fluWHO <- load.iiag.data.fluid(datadir = "C:/Users/hw3616/Desktop/Imperial/Project1_Forecasting/Project_Coding/iiag/data_old")

#' Extract weekly incidence for raw WHO FluID dataset
minprop <- 0.5
fluWHO.incidence <- extract.incidence.who(
  fluWHO,
  sel_iso3 = unique(fluWHO$ISO3),
  sel_ag = c("All"),
  sel_measure = c("ILI_CASES"),
  minYear=2010,
  maxYear = 2018)

sel_iso <- names(which(colSums(is.na(fluWHO.incidence))/dim(fluWHO.incidence)[1]<minprop))

#' countries whose data vailablity is over 50% are remained
fluWHO.incidence <- extract.incidence.who(
  fluWHO,
  sel_iso3 = sel_iso,
  sel_ag = c("All"),
  sel_measure = c("ILI_CASES"),
  minYear=2010,
  maxYear=2018
)

#' Check if the data of left countries is eligible for xgboost
#' check the data availablity in each year
#' exclude countries which do not have at least five consecutive weeks data in any year or do not have 
#' data records in 2010 or 2017 and 2018.
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

countryCode_no1718Or10 <- country_year$Country[which(country_year$end_year < 2017 | country_year$`2010`=="No"| 
                                                       country_year$`2018`=="No")]
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


#' Check if left countries have less than 5 weeks data in a year because I will do the 4-week ahead foreacast which
#' requires data of 5 consective weeks
#' Exclude countires whose datasets are uneligible to do the 4-week ahead forecast


# extract incidence of 41 eligible countries.
fluWHO.incidence <- extract.incidence.who(fluWHO,
                                          sel_iso_xgb,
                                          sel_ag = c("All"),
                                          sel_measure = c("ILI_CASES"),
                                          minYear = 2010,
                                          maxYear = 2018
                                          )


#### XGBooost ####

# example of data frame generated from gbm_complex() function
USA_complex1 <- gbm_complex(fluWHO.incidence, "USA", 10,1) # dataframe for 1-week ahead 
USA_complex2 <- gbm_complex(fluWHO.incidence, "USA", 10,2) # dataframe for 2-week ahead
USA_complex3 <- gbm_complex(fluWHO.incidence, "USA", 10,3) # dataframe for 3-week ahead
USA_complex4 <- gbm_complex(fluWHO.incidence, "USA", 10,4) # dataframe for 4-week ahead



#' check if data frames for 1,2,3,4-week ahead forecast are in the same size
length1 <- c()
for (i in 1:length(sel_iso_xgb)){
  a <- adjust.data.size(fluWHO.incidence, sel_iso_xgb[i], 10,1)
  tmp <- nrow(a)
  length1 <- append(length1,tmp)
  
}

length2 <- c()
for (i in 1:length(sel_iso_xgb)){
  a <- adjust.data.size(fluWHO.incidence, sel_iso_xgb[i], 10,2)
  tmp <- nrow(a)
  length2 <- append(length2,tmp)
}

length3 <- c()
for (i in 1:length(sel_iso_xgb)){
  a <- adjust.data.size(fluWHO.incidence, sel_iso_xgb[i], 10,3)
  tmp <- nrow(a)
  length3 <- append(length3,tmp)
}

length4 <- c()
for (i in 1:length(sel_iso_xgb)){
  a <- adjust.data.size(fluWHO.incidence, sel_iso_xgb[i], 10,4)
  tmp <- nrow(a)
  length4 <- append(length4,tmp)
}



# one-week ahead forecast
# 2010-2014 training, 2015 test
compare_pred15 <- compare_accuracy(sel_iso_xgb, 10, 0, 4 ,1)
length(which(compare_pred15$Accurate==1))/nrow(compare_pred15) # 0.7868296

# 2011-2015 training, 2016 test
compare_pred16 <- compare_accuracy(sel_iso_xgb,10, 1, 4, 1) 
length(which(compare_pred16$Accurate==1))/nrow(compare_pred16) # 0.7470383

# 2012-2016 training. 2017 test
compare_pred17 <- compare_accuracy(sel_iso_xgb,10, 2, 4, 1)
length(which(compare_pred17$Accurate==1))/nrow(compare_pred17)

# 2010-2016 traing, 2017 test
compare1016_pred17 <- compare_accuracy(sel_iso_xgb,10, 0, 6, 1)
length(which(compare1016_pred17$Accurate==1))/nrow(compare1016_pred17) # 0.6979038

# two-week ahead forecast
# 2010-2104 training, 2015 test
compare_pred15Two <- compare_accuracy(sel_iso_xgb, 10, 0, 4 ,2)

# 2011-2015 training, 2016 test
compare_pred16Two <- compare_accuracy(sel_iso_xgb, 10, 1, 4,2)

# 2012-2016 training. 2017&2018 test
compare_pred17Two <- compare_accuracy(sel_iso_xgb, 10, 2, 4, 2)

# 2010-2016 traing, 2017&2018 test
compare1016_pred17Two <- compare_accuracy(sel_iso_xgb, 10, 0, 6, 2)

# three week ahead
# 2010-2104 training, 2015 test
compare_pred15Three <- compare_accuracy(sel_iso_xgb, 10, 0, 4 ,3)

# 2011-2015 training, 2016 test
compare_pred16Three <- compare_accuracy(sel_iso_xgb, 10, 1, 4, 3)

# 2012-2016 training. 2017&2018 test
compare_pred17Three <- compare_accuracy(sel_iso_xgb, 10, 2, 4, 3)

# 2010-2016 traing, 2017&2018 test
compare1016_pred17Three <- compare_accuracy(sel_iso_xgb, 10, 0, 6, 3)

# four week ahead
# 2010-2104 training, 2015 test
compare_pred15Four <- compare_accuracy(sel_iso_xgb, 10, 0, 4 ,4)

# 2011-2015 training, 2016 test
compare_pred16Four <- compare_accuracy(sel_iso_xgb, 10, 1, 4, 4)

# 2012-2016 training. 2017&2018 test
compare_pred17Four <- compare_accuracy(sel_iso_xgb, 10, 2, 4, 4)

# 2010-2016 traing, 2017&2018 test
compare1016_pred17Four <- compare_accuracy(sel_iso_xgb, 10, 0, 6, 4)


# example of single tree
usa.tree.one <- xgboost.model.train(fluWHO.incidence,"USA",10,0,4,1)
usa.tr.one <- xgboost_dat(USA_complex1,2010,2014) 
usa.ts.one <- xgboost_dat(USA_complex1, 2015,2015)
usa.pred.one <- predict(usa.tree.one, usa.ts.one)
usa.xgboost.model.pred.output.one <- xgboost.model.pred.output(USA_complex1,2015,2015,usa.pred.one)

features <- attr(usa.ts.one, ".Dimnames")[[2]]
imp <- xgb.importance(features, usa.tree.one)
xgb.plot.importance(imp)


explainer <- buildExplainer(usa.tree.one,usa.tr.one, type="binary", base_score = 0.5, trees_idx = NULL)
pred.breakdown <- explainPredictions(usa.tree.one, explainer, usa.ts.one)
cat('Breakdown Complete','\n')

weights <- rowSums(pred.breakdown)
pred.xgb = 1/(1+exp(-weights))
cat(max(xgb.preds-pred.xgb),'\n')
idx_to_get = as.integer(802)
test[idx_to_get,-"left"]
showWaterfall(usa.tree.one, explainer, usa.ts.one, model.matrix(~.+0,USA_complex1[(195:242),]),1, type = "binary")



#### Heat plot for xgboost model ####

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

# example 
USA_hist_pred_one <- hist_average(fluWHO.incidence,"USA",10,1)
USA_hist_pred_two <- hist_average(fluWHO.incidence,"USA",10,2)
USA_hist_pred_three <- hist_average(fluWHO.incidence,"USA",10,3)
USA_hist_pred_four <- hist_average(fluWHO.incidence,"USA",10,4)


#' calculate the accuracy of historical average model
# 8050 / 11424 = 0.69
oneWeek_ahead_totalAccuracy_hist <- compare_accuracy_hist(fluWHO.incidence,sel_iso_xgb,10,1)
length(which(oneWeek_ahead_totalAccuracy_hist$Accurate==1))/nrow(oneWeek_ahead_totalAccuracy_hist)

# 0.68
twoWeek_ahead_totalAccuracy_hist <- compare_accuracy_hist(fluWHO.incidence,sel_iso_xgb,10,2)
length(which(twoWeek_ahead_totalAccuracy_hist$Accurate==1))/nrow(twoWeek_ahead_totalAccuracy_hist)

# 0.67
threeWeek_ahead_totalAccuracy_hist <- compare_accuracy_hist(fluWHO.incidence,sel_iso_xgb,10,3)
length(which(threeWeek_ahead_totalAccuracy_hist$Accurate==1))/nrow(threeWeek_ahead_totalAccuracy_hist)

# 0.66
fourWeek_ahead_totalAccuracy_hist <- compare_accuracy_hist(fluWHO.incidence,sel_iso_xgb,10,4)
length(which(fourWeek_ahead_totalAccuracy_hist$Accurate==1))/nrow(fourWeek_ahead_totalAccuracy_hist)

# plot the accuracy drop off as number of week-ahead increases
compareAccuracy_total_hist <- cbind(c(1:4),c(0.65,0.64,0.62,0.61)) %>% 
  as.data.frame()
colnames(compareAccuracy_total_hist) <- c("nWeek_ahead","percentage")

#### repeat model ####

# get the predictions by repeat model
USA_repeat_one <- repeat_model(fluWHO.incidence,'USA',10,1)
USA_repeat_two <- repeat_model(fluWHO.incidence,'USA',10,2)
USA_repeat_three <- repeat_model(fluWHO.incidence,'USA',10,3)
USA_repeat_four <- repeat_model(fluWHO.incidence,'USA',10,4)


# get the total accuract for 1,2,3,4-weeks ahead
oneWeek_ahead_totalAccuracy <- compare_accuracy_repeat(fluWHO.incidence,sel_iso_xgb,10,1)
length(which(oneWeek_ahead_totalAccuracy$Accurate==1))  # 9077
length(which(oneWeek_ahead_totalAccuracy$Accurate==1))/nrow(oneWeek_ahead_totalAccuracy) # 0.7945

twoWeek_ahead_totalAccuracy <- compare_accuracy_repeat(fluWHO.incidence,sel_iso_xgb,10,2)
length(which(twoWeek_ahead_totalAccuracy$Accurate==1)) # 8896
length(which(twoWeek_ahead_totalAccuracy$Accurate==1))/nrow(twoWeek_ahead_totalAccuracy) # 0.7440616

threeWeek_ahead_totalAccuracy <- compare_accuracy_repeat(fluWHO.incidence,sel_iso_xgb,10,3)
length(which(threeWeek_ahead_totalAccuracy$Accurate==1)) #8021
length(which(threeWeek_ahead_totalAccuracy$Accurate==1))/nrow(threeWeek_ahead_totalAccuracy) #0.7021183

fourWeek_ahead_totalAccuracy <- compare_accuracy_repeat(fluWHO.incidence,sel_iso_xgb,10,4)
length(which(fourWeek_ahead_totalAccuracy$Accurate==1)) # 7674
length(which(fourWeek_ahead_totalAccuracy$Accurate==1))/nrow(fourWeek_ahead_totalAccuracy) # 0.67

#' prepare the dataframe to plot the accuracy drop off as number of week-ahead increases
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