#' ** Make a database of season incidences from the WHO data

#' Set the current working directory if needed
## setwd("~/Dropbox/git/iiag/")

#' As always, remove all objects fromt eh workspace before starting
rm(list=ls(all=TRUE))

#' Source the function files and libraries we will need
source("src/R/riley_funcs.r")

#' Follow Cecile Viboud's script to load up data only from countries with data reorted 
#' for at least 50% of weeks
minprop <- 0.5
dataflu=load.iiag.data()
x <- extract.incidence(
		dataflu$synd,
		minYear=2010,
		sel_iso3 = unique(dataflu$synd$ISO3),
		sel_ag = c("All"),
		sel_measure = c("ILI_CASES"))
sel_iso=names(which(colSums(is.na(x))/dim(x)[1]<minprop))
x <- extract.incidence(
		dataflu$synd,
		minYear=2010,
		sel_iso3 = sel_iso,
		sel_ag = c("All"),
		sel_measure = c("ILI_CASES")
)

#' Load up the country info file
countrydesc <- read.table("data/country_list_ISO.csv", sep=',', header=TRUE)

#' Then extract individual seasons taking care with the week 53 / 52 issue
#' A year name for a season is the year in which the mid week appears, so,
#' for example, the 2017 northern year contains week 1 20017.
tmp <- extract.seasons(x,countrydesc)
inc <- tmp$inc
lu <- tmp$lu


