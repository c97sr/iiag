#' Make a note of when the report was generated.
Sys.time()

#' As always, remove all objects fromt the workspace before starting.
rm(list = ls(all = TRUE))

# CRAN
install.packages("cdcfluview")
# master branch
devtools::install_git("https://sr.ht/~hrbrmstr/cdcfluview")
devtools::install_git("https://gitlab.com/hrbrmstr/cdcfluview")
devtools::install_github("hrbrmstr/cdcfluview")

library("knitr")
library("xgboost")
library("stringr")

library(cdcfluview)
library(hrbrthemes)
library(tidyverse)

# current verison
packageVersion("cdcfluview")


fview_ILINet <- ilinet(region = "state")

# source("E:/Imperial/iiag/forecasting_vivi/fluView_funcs_vivi.R")
source("C:/Users/hw3616/Desktop/Imperial/Project1_Forecasting/Project_Coding/iiag/forecasting_vivi/fluView_funcs_vivi.R")

#' U.S. Outpatient Influenza-like Illness Surveillance Network (ILINet) consist of 
#' information on outpatient visits to health care providers for influenza-like illness
#' and the number of outpatinent confirmed as ILI.


#' Extract incidence data by states.
# region <- unique(fview_ILINet$REGION)

fview_incidence <- extract.incidence.fluView(fview_ILINet,
                                             sel_states = unique(fview_ILINet$region),
                                             minYear = 2010,
                                             maxYear = 2019)

minprop <- 0.5

#' None of satets contain entries between 2010-01 and 2010-39, so delete these rows before calculate the avaibility
#' of dataset
fview_incidence <- fview_incidence[-c(1:39),]

us_states <- names(which(colSums(is.na(fview_incidence))/dim(fview_incidence)[1]<minprop))

#' Check if the data of remained states are eligible for xgboost
#' check the data availablity in each year
#' exclude states which do not have at least five consecutive weeks data in any year or do not have 
#' entries in certain years.
state_year <- NULL
for (i in 1:length(us_states)){
  tmp <- duration(us_states[i])
  state_year <- rbind(state_year, tmp)
}
state_year <- as.data.frame(state_year)
colnames(state_year) <- c("State","2010","2011","2012","2013","2014","2015",
                            "2016","2017","2018","2019", "start_year","end_year")
state_year$end_year <- as.numeric(as.character(state_year$end_year))
state_year$start_year <- as.numeric(as.character(state_year$start_year))
rownames(state_year) <- c(1:nrow(state_year))

state_no19Or10 <- state_year$State[which(state_year$`2019`=="No" | state_year$`2010`=="No")]
                                            
state_no10 <- state_year$State[which(state_year$`2010`=="No")]

state_no19Or10 <- us_states[which(us_states %in% state_no19Or10)] 

#' List of US states will be included in the analysis.
us_xgb <- us_states
for (i in 1:length(state_no19Or10)){
  index <- which(us_xgb == state_no19Or10[i])
  us_xgb <- us_xgb[-(index)]
}


#' Check if left states have less than 5 weeks data in a year because I will do the 4-week ahead foreacast which
#' requires data of 5 consective weeks
#' Exclude states whose datasets are uneligible to do the 4-week ahead forecast
#' So far I (Vivi) did this manually, will write a function to check automatically.

# extract incidence of eligible states.
fview_incidence <- extract.incidence.who(fview_ILINet,
                                         us_xgb,
                                         minYear = 2010,
                                         maxYear = 2019
)

#' XGB model

