#' Make a note of when the report was generated.
Sys.time()

#' As always, remove all objects fromt the workspace before starting.
rm(list = ls(all = TRUE))


library("knitr")
library("xgboost")
library("stringr")

# source("E:/Imperial/iiag/forecasting_vivi/fluView_funcs_vivi.R")
source("C:/Users/hw3616/Desktop/Imperial/Project1_Forecasting/Project_Coding/iiag/forecasting_vivi/fluView_funcs_vivi.R")

#' U.S. Outpatient Influenza-like Illness Surveillance Network (ILINet) consist of 
#' information on outpatient visits to health care providers for influenza-like illness
#' and the number of outpatinent confirmed as ILI.

#' First of all, fix the headers of dataset to the form that won's cause error in the Rstudio
fview_ILINet <- load.iiag.data.fluView("C:/Users/hw3616/Desktop/Imperial/Project1_Forecasting/Project_Coding/iiag/forecasting_vivi/fluView_data/")

# fview_ILINet <- load.iiag.data.fluView("E:/Imperial/iiag/forecasting_vivi/fluView_data")

#' Extract incidence data by states.
region <- unique(fview_ILINet$REGION)
fview_ILINet2 <- fview_ILINet[order(fview_ILINet$REGION),]

