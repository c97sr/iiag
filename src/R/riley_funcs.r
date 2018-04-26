#' Load up the data as two large dataframes.
#' Possible refinements are to add in date ranges and to get it to construct
#' older snapshots of the data using the date string or prevvious version,
#' next-previous and so on.
load.iiag.data <- function(datadir="data") {
  
  ## Define the strings for all the files needed and read in the 6 files
  fid_this <- read.csv(paste(datadir,"/2017-2018_FluIDData.csv",sep=""),
                       header=FALSE)
  fid_old_0 <- read.csv(paste(datadir,"/2000-2009_FluIDData.csv",sep=""))
  fid_old_1 <- read.csv(paste(datadir,"/2010-2013_FluIDData.csv",sep=""))
  fid_old_2 <- read.csv(paste(datadir,"/2014-2016_FluIDData.csv",sep=""))
  fnet_this <- read.csv(paste(datadir,"/2017-2018_FluNetData.csv",sep=""))
  fnet_old_0 <- read.csv(paste(datadir,"/2000-2009_FluNetData.csv",sep=""))
  fnet_old_1 <- read.csv(paste(datadir,"/2010-2013_FluNetData.csv",sep=""))
  fnet_old_2 <- read.csv(paste(datadir,"/2014-2016_FluNetData.csv",sep=""))
  
  ## There is a slight issue with the name in the current versions of the data
  ## These lines make sure that the column title for the ISO3 column is
  ## consistent.
  ## Below also takes care of a missing header row for current snapshot
  names(fid_old_1)[1] <- names(fnet_this)[1]
  names(fid_old_2)[1] <- names(fnet_this)[1]
  names(fnet_old_1)[1] <- names(fnet_this)[1]
  names(fnet_old_2)[1] <- names(fnet_this)[1]
  names(fid_this) <- names(fid_old_1)
  
  ## Use rbind to make the large tables. Should throw an error if the column
  ## names change in the future.
  dfIdXX_dev <- rbind(fid_old_0,fid_old_1,fid_old_2,fid_this)
  dfNet_dev <- rbind(fnet_old_0,fnet_old_1,fnet_old_2,fnet_this)
  
  ## Not currently including the 00 to 09 data as it breaks the test plots
  dfIdXX <- rbind(fid_old_1,fid_old_2,fid_this)
  dfNet <- rbind(fnet_old_1,fnet_old_2,fnet_this)
  
  ## Return the two datasets as a list
  list(lab=dfNet,synd=dfIdXX,labDev=dfNet_dev,syndDev=dfIdXX_dev)
  
}

## Take the raw WHO country database and extract weekly incidence
extract.incidence <- function(
  dfId,
  sel_iso3,
  sel_ag,
  sel_measure,
  minYear = 2000,
  maxYear = 2018
) {
  
  ## Setup the week scale in a format consistent with the week format
  ## in the data and cope with 53-week years. Needs the list of 53 week years
  ## extending in both directions.
  ## Perhaps should have a few lines to get rid of data NAs and avoid a warning
  ## at the next line?
  dfId$yrweek <- paste(dfId$ISO_YEAR,sprintf("%02d",as.numeric(dfId$ISO_WEEK)),sep="-")
  min(dfId$ISO_YEAR)
  yrs53Weeks <- c(2009,2015,2020)
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

## Script used for development and testing
## Think about moving this to be in the notes directory and use spin
## Next need to get the historical data working and transfer these notes here to a
## separate stad alone script file
if (FALSE) {
  
  ## Clear objects form memory for debugging and source this file
  ## setwd("~/Dropbox/git/iiag/")
  rm(list=ls(all=TRUE))
  source("src/R/riley_funcs.r")
  
  ## Assumes running in src/R. Change datadir as needed. Select the
  ## syndromic data after loading
  tmp <- load.iiag.data(datadir="data")
  df <- tmp$synd
  dfDev <- tmp$syndDev
  
  ## Need a bit more of a play to see what is there from before 2010
  dim(df)
  dim(dfDev)
  
  ## Extract ILI cass for UK, USA and germany for all age groups
  ## Doesn't seem quite right at the moment because of the blanks and
  ## zeros for the UK and germany. But need to eyeball specific lines of
  ## the data to be sure its a real problem. Data for other countries seems
  ## to come out OK.
  
  ## This seems to be broken now for some reason?
  ## That must be current
  x <- extract.incidence(
    df,
    minYear=2010,
    sel_iso3 = c("GBR","USA","DEU"),
    sel_ag = c("All"),
    sel_measure = c("ILI_CASES")
  )
  
  ## Quick diagnostic plot
  plot(x[,"GBR"]+1,type="l",col="red",ylim=c(0,max(x,na.rm=TRUE)),ylog=TRUE)
  points(x[,"DEU"]+1,type="l",col="blue")
  points(x[,"USA"]+1,type="l",col="green")
  
}
