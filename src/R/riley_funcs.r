#' Load up the data as two large dataframes.
#' Possible refinements are to add in date ranges and to get it to construct
#' older snapshots of the data using the date string or prevvious version,
#' next-previous and so on.
load.iiag.data <- function(datadir="../iiag_data/data") {

    ## Define the strings for all the files needed and read in the 6 files
    fid_this <- read.csv(paste(datadir,"/2017-2018_FluIDData.csv",sep=""))
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

#' Setup to inport the latest data. Looks a little cumbersome, but is
#' probably better done like this to be explicit about every file. It will only need
#' changing ~ once a year.
load.iiag.data.new <- function(datadir="../iiag_data/data") {

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

    ## Define the strings for all the files needed and read in the 6 files
    fid_this <- read.csv(paste(datadir,"/2019-2020_FluIDData.csv",sep=""))
    fid_old_3 <- read.csv(paste(datadir,"/2017-2018_FluIDData_20190110.csv",sep=""))
    fid_old_2 <- read.csv(paste(datadir,"/2014-2016_FluIDData.csv",sep=""))
    fid_old_1 <- read.csv(paste(datadir,"/2010-2013_FluIDData.csv",sep=""))
    fid_old_0 <- read.csv(paste(datadir,"/2000-2009_FluIDData.csv",sep=""))
    fnet_this <- read.csv(paste(datadir,"/2017-2018_FluNetData_20190110.csv",sep=""))
    fnet_old_3 <- read.csv(paste(datadir,"/2017-2018_FluNetData.csv",sep=""))
    fnet_old_2 <- read.csv(paste(datadir,"/2014-2016_FluNetData.csv",sep=""))
    fnet_old_1 <- read.csv(paste(datadir,"/2010-2013_FluNetData.csv",sep=""))
    fnet_old_0 <- read.csv(paste(datadir,"/2000-2009_FluNetData.csv",sep=""))

    ## Fix names for flu id
    names(fid_this) <- fix_headers(fid_this)
    names(fid_old_1) <- fix_headers(fid_old_1)
    names(fid_old_2) <- fix_headers(fid_old_2)
    names(fid_old_3) <- fix_headers(fid_old_3)
    names(fnet_old_1) <- fix_headers(fnet_old_1)
    names(fnet_old_2) <- fix_headers(fnet_old_2)
    names(fnet_old_3) <- fix_headers(fnet_old_3)
    names(fnet_this) <- fix_headers(fnet_this)

    ## Use rbind to make the large tables. Should throw an error if the column
    ## names change in the future.
    ## dfId <- rbind(fid_old_0,fid_old_2,fid_old_2,fid_old_3,fid_this)
    ## dfNet <- rbind(fnet_old_0,fnet_old_2,fnet_old_2,fnet_old_3,fnet_this)
    dfId <- rbind(fid_old_1,fid_old_2,fid_old_3,fid_this)
    dfNet <- rbind(fnet_old_1,fnet_old_2,fnet_old_3,fnet_this)

    ## Sort both dataframe just incase
    dfId <- dfId[order(dfId$ISO2,dfId$ISO_YEAR,dfId$ISO_WEEK),]
    dfNet <- dfNet[order(dfNet$ISO2,dfNet$ISO_YEAR,dfNet$ISO_WEEK),]

    ## Return the two datasets as a list
    list(lab=dfNet,synd=dfId)

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

## Function to create a database of seasons from the incidence data
extract.seasons <- function(
		inc_arr,
		country_info,
		midwknorth=1,
		midwksouth=27,
		midleadweeks=8,
		midlagweeks=16) {

  ## Extract the number of seasons from the database


  tmp <- unlist(lapply(strsplit(rownames(inc_arr),split="-"),as.numeric))
  noelements <- length(tmp)
  noweeks <- noelements / 2
  yrwk <- matrix(tmp,nrow=noweeks,ncol=2,byrow=TRUE)
  colnames(yrwk) <- c("yr","wk")
  veccountries <- colnames(inc_arr)

  ## Setup the return data structures.
  ## We need to create a lookup table of seasons first, then create a 2D matrix of the
  ## correct form and then run through eac country putting in the
  ## correct data. It should have country ID, season type, with 1 for northern and 2
  ## for southern. and then
  weeksperseason <- 1 + midleadweeks + midlagweeks

  rtnarr <- matrix(nrow=0,ncol=weeksperseason)
  rtnlu <- data.frame(country=NULL, year=NULL, seastype=NULL, nonas=NULL, nozeros=NULL, mean=NULL, total=NULL, max=NULL)
  tmprow <- data.frame(country=NA, year=NA, seastype=NA, nonas=NA, nozeros=NA, mean=NA, total=NA, max=NA)

  for (lab in veccountries) {

	  ## Give a season type based on information about the country
	  ## for now just northern or southern
	  cind <- match(lab,country_info$ISO3,nomatch=-1)
	  if (cind < 0) stop("failed to match country code")
	  clat <- country_info$Latitude[cind]
	  if (clat > 0) {

      ## Northern hemisphere
      midweek <- midwknorth
      seasontype <- 1

	  } else if (clat < 0) {

      ## Southern hemisphere
      midweek <- midwksouth
      seasontype <- 2

	  } else {

      # But we could define custom weeks for each season
      stop("Nothing here right now.")

    }

    ## Main loop for the cutting out the seasons
	  for (w in 1:noweeks) {
		  if (yrwk[w,"wk"] == midweek) {

        if (((w - midleadweeks) > 0) && (w + midlagweeks <= noweeks)) {

          ## Extract the incidence for the correct weeks
          thisyrinc <- inc_arr[(w - midleadweeks):(w + midlagweeks), lab]

          ## Add the season details to the season lookup table
          ## Can easily add in more season stats here
          tmprow$country <- lab
          tmprow$year <- yrwk[w,"yr"]
          tmprow$seastype <- seasontype
          tmprow$nonas <- sum(is.na(thisyrinc))
          tmprow$nozeros <- sum(thisyrinc < 0.1, na.rm=TRUE)
          tmprow$mean <- mean(thisyrinc,na.rm=TRUE)
          tmprow$total <- sum(thisyrinc,na.rm=TRUE)
          if (tmprow$nonas < weeksperseason) {
            tmprow$max <- max(thisyrinc,na.rm=TRUE)
          } else {
            tmprow$max <- -1
          }

          rtnlu <- rbind(rtnlu,tmprow)

          ## Add the incidence itself to the incidence array
          ## Use rbind to make sure we do only one pass and to get match row sizes
          rtnarr <- rbind(rtnarr, thisyrinc)

        }
      }
	  }

  }

  list(inc = rtnarr, lu = rtnlu)

}

#' A function to take the two part database and return a forecast for that time
#' By crude, I mean, really crude. This needs to at least ave a precalc function
#' that creates a matrix that is then called. Ah well.
fm.null.hist.vvcrude <- function(ctry,yr,wk,lutab,inctab) {

  ## Extract a vector of incidence for this week
  tmpinc <- inctab[(lutab$country==ctry),wk]
  tmplu <- lutab[(lutab$country==ctry),]
  yrind <- match(yr,tmplu$year,nomatch=-1)
  if (yrind < 0) stop("didn't match the year properly")
  vecotheryrs <- tmpinc[-yrind]
  if (sum(is.na(vecotheryrs)==length(vecotheryrs))) {
    stop("Need to think about this a bit")
  }
	rtn <- mean(vecotheryrs,rm.na=TRUE)

  ## Return the simple estimate
  rtn

}

ac.skill.crude <- function(tabObs,tabCast,tablu,tol=0.2,thresh=0.05) {

  ## Up to debugging the lines below
  ## browser()

  # Set thresholds for the different countries
  inccountries <- names(table(tablu$country))

  0

}

## Code fragments for interactive sessions
if (FALSE) {

  ## Can put some dev code here, but need to move it to the notes sirectory as soon as
  ## possible. Nothing here right now.

}


