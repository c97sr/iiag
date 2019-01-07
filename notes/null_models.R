#' # Make a database of season incidences from the WHO data

#' ## Preamble to setup the session

#' Set the current working directory if needed when working interactively and
#' check the working directry for when spinning.
## setwd("~/Dropbox/git/iiag/notes")
getwd()

#' As always, remove all objects fromt eh workspace before starting
rm(list=ls(all=TRUE))

#' Source the function files and libraries we will need
source("../src/R/riley_funcs.r")

#' ## Load up the datasets
#' Follow Cecile Viboud's script to load up data only from countries with data reorted 
#' for at least 50% of weeks, even before we look at the properties of individual seasons
#+ tidy = TRUE 
minprop <- 0.5
dataflu=load.iiag.data(datadir="../data")
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
countrydesc <- read.table("../data/country_list_ISO.csv", sep=',', header=TRUE)

#' ## Create a database of seasons
#' Then extract individual seasons taking care with the week 53 / 52 issue
#' A year name for a season is the year in which the mid week appears, so,
#' for example, the 2017 northern year contains week 1 20017.
tmp <- extract.seasons(x,countrydesc)
tmpinc <- tmp$inc
tmplu <- tmp$lu

#' Then look at the distribution of epidemic properties.
table(tmplu$seastype)
hist(tmplu$total, breaks = c(-0.5,0.5,9999999999),plot=FALSE)
plot(1+tmplu$total,1+tmplu$max,log="xy")
hist(tmplu$nozeros,breaks=seq(-0.5,26.5))
hist(tmplu$nonas,breaks=seq(-0.5,26.5))

#' And use these to decide on some cutoffs then create a filtered table
filtermask <- (tmplu$nonas < 11) & (tmplu$nozeros < 11) & (tmplu$total > 300) 
sum(filtermask)
tmp2lu <- tmplu[filtermask,]
tmp2inc <- tmpinc[filtermask,]
dim(tmp2inc)

#' Now check for the number of countries and filter for only those with enough
#' good seasons. This can e melded up into the other criteria above (maybe)
countrycounts <- table(tmp2lu$country)
ctryinc <- names(countrycounts[countrycounts > 5])
length(ctryinc)
fmask2 <- tmp2lu$country %in% ctryinc
lu <- tmplu[fmask2,]
obs <- tmpinc[fmask2,]
noepi <- dim(lu)[1]
noweeks <- dim(obs)[2]
dim(obs)

#' ## Define a null historical model
#' Using the general
#' approach that a forecast function only returns a point estimate at this
#' stage and the first two arguments must be country and week.
#' Quickly check that it works for one week!
fm.null.hist.vvcrude("ISL",2015,15,lu,obs)

#' Then create an empty matrix of the same shape as the observation matrix and
#' cycle through it making a forecast
cast <- matrix(nrow=nrow(obs),ncol=ncol(obs))
for (r in 1:noepi) {
  curctry <- lu$country[r]
  curyear <- lu$year[r]
  for (w in 1:noweeks) {
    cast[r,w] <- fm.null.hist.vvcrude(curctry,curyear,w,lu,obs) 
  }
}

#' ## Define a crude skill measure
#' Now need to assign accurate or not to every place where we made a forecast
#' Here the NAs have to come back
skill <- ac.skill.crude(obs,cast,lu)