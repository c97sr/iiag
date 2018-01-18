load.iiag.data <- function(datadir="../../data") {

    ## Load up the key variables
    fid_this <- read.csv(paste(datadir,"/2017-2018_FluIDData.csv",sep=""))
    fid_old_1 <- read.csv(paste(datadir,"/2010-2013_FluIDData.csv",sep=""))
    fid_old_2 <- read.csv(paste(datadir,"/2014-2016_FluIDData.csv",sep=""))

    fnet_this <- read.csv(paste(datadir,"/2017-2018_FluNetData.csv",sep=""))
    fnet_old_1 <- read.csv(paste(datadir,"/2010-2013_FluNetData.csv",sep=""))
    fnet_old_2 <- read.csv(paste(datadir,"/2014-2016_FluNetData.csv",sep=""))

    ## There is a slight issue with the name, as read, of the first column in id.
    ## will use the name from the current data, to allow cbind. Woorth tidying up
    ## at some point in the data
    names(fid_old_1)[1] <- names(fid_this)[1]
    names(fid_old_2)[1] <- names(fid_this)[1]
    dfIdXX <- rbind(fid_old_1,fid_old_2,fid_this)
    names(fnet_old_1)[1] <- names(fnet_this)[1]
    names(fnet_old_2)[1] <- names(fnet_this)[1]
    dfNet <- rbind(fnet_old_1)

    list(net=dfNet,id=dfIdXX)
    
}

extract.incidence <- function(
                              dfId,
                              sel_iso3,
                              sel_ag,
                              sel_measure,
                              minYear = 2010,
                              maxYear = 2018
) {
    
    ## Below here can form a function to return a complete set of week labels
    dfId$yrweek <- paste(dfId$ISO_YEAR,sprintf("%02d",dfId$ISO_WEEK),sep="-")
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
    
    ## Now need to use the week labels, and a country code to extract sums
    ## sums of incidence
    sel_weeks <- vecWeekScale
    rtnmat <- matrix(data=NA,nrow=length(sel_weeks),ncol=length(sel_iso3))
    colnames(rtnmat) <- sel_iso3 
    
    for (cur_iso3 in sel_iso3) {
        
        crit1 <- (dfId$ISO3 == cur_iso3)
        crit2 <- (dfId$AGEGROUP_CODE %in% sel_ag) 
        crit3 <- (dfId$MEASURE_CODE %in% sel_measure)
        tmpdf <- dfId[crit1 & crit2 & crit3,]
        tmpdf <- tmpdf[order(tmpdf$yrweek),]
        max_ind_rtn <- dim(rtnmat)[1]    
        max_ind_df <- dim(tmpdf)[1]

        cur_ind_rtn <- 1
        cur_ind_df <- 1   
        while (cur_ind_df <= max_ind_df) {
            while (
                sel_weeks[cur_ind_rtn] != tmpdf$yrweek[cur_ind_df] &&
                cur_ind_rtn <= max_ind_rtn
            ) {
                cur_ind_rtn <- cur_ind_rtn + 1
            }
            if (cur_ind_rtn <= max_ind_rtn) {
                val_rtn <- rtnmat[cur_ind_rtn,cur_iso3]
                val_df <- tmpdf$ValueNumeric[cur_ind_df]
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
    }

    rtnmat

}

if (FALSE) {
    
    ## Clear objects form memory for debugging and source this file
    rm(list=ls(all=TRUE))
    source("riley_funcs.r")

    ## Assumes running in src/R. Change datadir arg if needed.
    tmp <- load.iiag.data()
    df <- tmp$id
    
    ## Extract ILI cass for UK, USA and germany for all age groups
    ## Doesn't seem quite right at the moment
    x <- extract.incidence(
        df,
        sel_iso3 = c("GBR","USA","DEU"),
        sel_ag = c("All"),
        sel_measure = c("ILI_CASES")
    )

    ## Quick diagnostic plot of the incidence
    plot(x[,"GBR"]+1,type="l",col="red",ylim=c(0,max(x,na.rm=TRUE)),ylog=TRUE)
    points(x[,"DEU"]+1,type="l",col="blue")
    points(x[,"USA"]+1,type="l",col="green")
    
    ## tmpdf1 <- dfIdXX[dfIdXX$ISO3=="DEU" & dfIdXX$MEASURE_CODE=="ILI_CASES",] 
    ## table(tmpdf1$AGEGROUP_CODE,tmpdf1$ISO_YEAR)
    
}
