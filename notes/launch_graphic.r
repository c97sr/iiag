#' # Summary heat chart of fluID data
#' This R file can be converted into a nicely formatted html page with a
#' figures subdirectory by using the command "spin("this_file_name")" from
#' the package knitr.

#' Taken almost directly from a script written by Cecile Viboud in early 2018.

#' Make a note of when the report was generated.
Sys.time()

#' Set the current working directory if needed when working interactively and
#' check the working directry for when spinning. To spin this, it needs the
#' knitr package and then spin("thisfilename.r").
## setwd("~/Dropbox/git/iiag/notes")

#' As always, remove all objects fromt the workspace before starting.
rm(list=ls(all=TRUE))

start_year <- 2010
end_year <- 2019
min_prop_reported <- 0.5

#' Source the function files and libraries we will need. Install SR's idd package
#' from github if needed.
# install_github("c97sr/idd")
library(knitr)
library("idd")
library("devtools")
library("formatR")
source("../src/R/riley_funcs.r")
source("../src/R/viboud_funcs.r")

#' Extract the data from the data directory
datadir <- "../../iiag_data/data"
dataflu <- load.iiag.data.new(datadir=datadir)
datafluold <- load.iiag.data(datadir=datadir)
## dataflu_new <- load.iiag.data.new(datadir="E:/Dropbox/data/who_flu")

#' Extract data for all countries
x <- extract.incidence(
    dataflu$synd,
    minYear = start_year,
    maxYear = end_year,
    sel_iso3 = unique(dataflu$synd$ISO3),
    sel_ag = c("All"),
    sel_measure = c("ILI_CASES"))

#' Identify countries with robust ILI by looking at sum of ILI cases by country.
#' Compare old and new data
colSums(is.na(x))

#' Calculate the proportion of missing weeks by country
hist(colSums(is.na(x))/dim(x)[1],breaks=seq(0,1,by=1/20))

#' List of countries with more than 50% of weeks of complete data
sel_iso=names(which(colSums(is.na(x))/dim(x)[1] < min_prop_reported))

#' Reextract data for set of more complete countries
x <- extract.incidence(
    dataflu$synd,
    minYear = start_year,
    maxYear = end_year,
    sel_iso3 = sel_iso,
    sel_ag = c("All"),
    sel_measure = c("ILI_CASES"))

#' Total number of ILI episodes reported here are
sum(x,na.rm=TRUE)

#' Setup a date timeline
datee=seq(from=start_year,by=1/52.17, length.out=dim(x)[1])

#' Load up the country data for latitude etc
countrydesc=read.table(
    paste(datadir,"/country_list_ISO.csv",sep=""), sep=',',header=TRUE
)
isos=colnames(x)
x2=x[,order(isos)] ## sorting countries in ILI matrix by ISO3  for later use
countries=countrydesc[ (countrydesc$ISO3 %in% isos), ]
idx.lat=order(countries$Latitude)

#' Setup the timeseries for the heat map
nwk=dim(x)[1]
x2.sc=apply(x2, 2, scale)
colnames(x2.sc)

#' Plot the mian chart
pdf_fig_proof(file="./figure/iiag_heatmap.pdf",pw=13.2,ph=18)
  nbreaks=11 ## nb of color breaks; seems to work for most flu examples
  breaksu=c(min(x2.sc,na.rm=T),-1.5,-1,-.5, 0, .5, 1,1.5, 2,4,6,10)
  image(x=datee, z=as.matrix(x2.sc[,idx.lat]),
        col=heat.colors(nbreaks)[nbreaks:1],
        breaks=breaksu,yaxt="n",xlab="Weeks")
  axis(2, at=seq(0,1,,dim(x2.sc)[2]),
       labels=colnames(x2.sc)[idx.lat],las=1,cex.axis=1.0)
  abline(v=seq(2010,2018),lty=1,col="grey",xpd=TRUE)
dev.off()

#' Plot the scale
pdf_fig_proof(file="./figure/iiag_heatmap_scale.pdf",pw=2,ph=7)
  image.scale(x2.sc, col=heat.colors(nbreaks)[nbreaks:1], breaks=breaksu,
    horiz=FALSE, yaxt="n", xaxt="n", xlab="", ylab="")
  axis(4)
  mtext("Scaled ILI cases", side=4, line=2)
  box()
dev.off()
