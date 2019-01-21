#' # Make a database of season incidences from the WHO data

#' ## Preamble to setup the session

#' Make a note of when the report was generated.
Sys.time()

#' Set the current working directory if needed when working interactively and
#' check the working directry for when spinning. To spin this, it needs the 
#' knitr package and then spin("thisfilename.r").
## setwd("~/Dropbox/git/iiag/notes")

#' As always, remove all objects fromt the workspace before starting.
rm(list=ls(all=TRUE))

#' Source the function files and libraries we will need. Install SR's idd package
#' from github if needed.
# install_github("c97sr/idd")
library(knitr)
library("idd")
library("devtools")
library("formatR")
source("../src/R/riley_funcs.r")
