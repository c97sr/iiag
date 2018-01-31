ave_seas_plots <- function(
                           df, # dataframe with output from other routines
                           zones, # dataframe defining transmission zones
                           plotCategory = "country", #"zone" or "country"
                           nameArea = "USA", #Select country
                           YearStart = 2016) {

    require(gdata)
    require(reshape2)
    require(ggplot2)
    require(plyr)
    require(zoo)
    require(epitools)
    require(Deducer)
    require(stringr)

    ## Path and File Definition #Change file directories according
    ## to your own settings
    baseFolder <- "~/Global Burden/Julia F"
    setwd(baseFolder)
    dataFile <- "dfNet.csv" #Option to import

    ##Fill these in as required:
    ##Need object linking country/ITZ/hemisphere
    ## SR These now taen as arguments
    # plotCategory <- "country" #"zone" or "country"
    # nameArea <- "USA" #Select country
    # YearStart <- 2016


    ## Zones needed below here
    if (plotCategory=="country"){
        ISO3ind <- nameArea #Select country
        ITZind <- zones$itz[zones$country==ISO3ind]
        HEMIS <- zones$hemis[zones$country==ISO3ind]
    } else {
        ITZind <- nameArea
        countries <- zones$country[zones$itz==ITZind]
        countries <- levels(countries)
        ISO3ind <- nameArea #Treat as single country in rest of code
    }

    if (HEMIS=="northern"){
        NOW1 <- toString(YearStart)
        NOW2 <- toString(YearStart+1)
        NOW <- paste(NOW1, "/", NOW2, sep="")
    } else {
        NOW <- toString(YearStart)
        NOW2 <- NOW
    }

    ################################
    ## Or just run from here, inputting the following manually:
    ##(country only)
    ## dataFile <- "dfNet.csv" #Option to import
    ## SR df and zones now assumed to be an argument
    ## df <- read.csv(dataFile,as.is=T)                 
    ## plotCategory <- "country" #"zone" or "country"
    ## ISO3ind <- "DEU" #Select country
    ## ITZind <- "ITZ" #Select ITZ
    ## HEMIS <- "northern" #"northern"/"southern"/"tropical"
    ## NOW <- "2016/2017"
    ## NOW2 <- "2017" #The later year of current season
    ## YearStart <- 2016 #First year in current season
    ################################

    ## measures <- c("SPEC_PROCESSED_NB", "DECTECTED_NB")
    ## res <- res[which(res$MEASURE_CODE %in% measures), ]
    df <- df[which(df$MEASURE_CODE=="SPEC_PROCESSED_NB" |
                   df$VIRUSTYPE_CODE=="ALL_INF"), ]

    if (plotCategory=="country"){
        df2 <- df[df$ISO3==ISO3ind, ]
        countryName <- df2[1,3]
    } else{ 
        df2 <- df[df$ISO3 %in% countries, ]
        countryName <- nameArea #Better versions of names?
    }
    
    df2 <- data.frame(df2$ISO_YEAR, df2$ISO_WEEK, df2$MEASURE_CODE, df2$ValueNumeric)
    colnames(df2) <- c("year", "week", "code", "value")
    
    if (plotCategory=="zone"){
        ## df2 < -aggregate(df2, by=list(year=df2$year, week=df2$week, code=df2$code), mean, na.rm=TRUE)
        dfx <- ddply(df2, c("year", "week", "code"), fun=mean)
    }
    
    df2 <- reshape(df2, idvar=c("year", "week"), timevar=c("code"), direction="wide")
                                        #Should land in chronological order ********
    df2$percPos <- df2[, 3]/df2[, 4]*100
    df2$percPosMov <- rollapply(df2$percPos, 3, mean, align="center", fill=T)
    df2 <- df2[, c(1,2,5,6)]
    ##Find peaks:
    if (HEMIS=="northern"){
        df2$syear <- df2$year
        these <- which(df2$week<40)
        df2$syear[these] <- df2$year[these]-1
        df2$sweek <- (df2$week-39) %% 53 #Year of 1 = season start year
    } else {
        df2$syear <- df2$year
        df2$sweek <- df2$week #Includes tropics for now ********
    }
    df2$sweek[df2$sweek==0] <- 53#sweek 14 not always there
                                        #
    df3 <- data.frame(df2$syear, df2$sweek, df2$percPos)
    colnames(df3) <- c("syear", "sweek", "mov")
                                        #df3 <- ddply(df3, "syear")
    df3 <- reshape(df3, idvar="syear", timevar="sweek", direction="wide")
                                        #
    df3[df3=="NaN"]=NA #Discovered in Fiji
                                        #
    rownames(df3) <- df3[, 1]
    df3 <- df3[, -1]
    nWeeks <- ncol(df3)
    if (HEMIS=="northern"){
        df3 <- df3[, c(seq(40,nWeeks), seq(1,39))] #Generalise, just in case ********
    }
    peakVal <- apply(df3, 1, max, na.rm=TRUE)
    peakInd <- max.col(df3[, seq(1,nWeeks)], "first")
    latest <- max(peakInd)
    ##Deal with week 53 (14):
    w14 <- which(is.na(df3[, 14]))
    df3[w14,seq(14,52)] <- df3[w14, seq(15,52)]
    df3[w14,53]=NA
    ##Allign to latest/mean over all:
    nyears <- nrow(df3)
    maxVal <- c(nyears,1,0)
    maxInd <- maxVal
    ##There must be a more elegant way than this loop
    for (i in 1:nyears){
        maxi <- max(df3[i, ], na.rm=TRUE)
        indi <- which(df3[i, ]==maxi)
        indi <- indi[1]
        maxVal[i] <- maxi
        maxInd[i] <- indi
    }
    maxShift <- max(maxInd)-min(maxInd)
    chop <- max(maxInd)-median(maxInd)
    peakVec <- maxInd#+median(maxInd)-min(maxInd) #as.data.frame(maxInd+min(maxInd))
    peakYears <- rownames(df3) #rownames(peakVec) <- rownames(df3)
    shiftWeeks <- max(maxInd)-maxInd
    ##df4 <- shift(df4, shiftWeeks, fill=NA, type="lag")#, give.names=FALSE)
    ##Clonky again - use "shift" as above? - Don't need "maxShift"
    df4 <- matrix(NA, nyears, (maxShift+53))
    for (i in 1:nyears){
        si <- shiftWeeks[i]
        df4[i, c(seq((si+1), (si+53)))] <- as.numeric(df3[i, ])
    }
    df4 <- as.data.frame(df4)
    rownames(df4) <- rownames(df3)
    yearList <- as.numeric(rownames(df4))
    currentYear <- which(yearList==YearStart)
    pastAv <- as.numeric(apply(df4[seq(2, (currentYear-1)), ], 2, mean, na.rm=TRUE))
    pastSd <- as.numeric(apply(df4[seq(2, (currentYear-1)), ], 2, sd, na.rm=TRUE))
    nVals <- apply(df4, 2, function(x) length(which(!is.na(x))))
    currentYear <- as.numeric(df4[currentYear, ])
    
                                        #Plot:
    endChop <- 5
    xmax <- length(currentYear)
    ymax <- min( (c( max(currentYear, na.rm=TRUE), max(pastAv, na.rm=TRUE))+20), 100, na.rm=TRUE)
    tvec <- seq((chop+1), (xmax-endChop))
    plot(pastAv[tvec], type="l",xlim=c(1, length(tvec)), ylim= c(0, ymax), ylab="Influenza positivity (%)",
         xlab="Weeks", main=paste(countryName, NOW), bty="l", cex.main=1, cex=2, cex.lab=1, cex.axis=1)#, xaxt="n", yaxs='i', xaxs='i')
    points(currentYear[tvec], type="l", lwd=2)
                                        #lvals=length(zoneOld)
                                        #Error bars:
    errorP <- pastAv+1.645*pastSd/nVals
    errorM <- pastAv-1.645*pastSd/nVals
    points(errorP[tvec], type="l", lty = "dashed")
    points(errorM[tvec], type="l", lty = "dashed")
    
    
                                        # From here - DHy1 (past average) - correct?
    boxplot(peakVec, horizontal=T, add=T, at=(ymax-10), axes=F)
    
    yVals <- jitter(rep((ymax-5),length(peakVec)), 2)
    points(peakVec, yVals, pch=18, cex=1)
    text(peakVec, yVals, labels=peakYears, pos=3)
                                        #axis(side=1, at=seq(1,xmax), labels=peakYears, cex.axis=1, par(xpd=T))
    
    legendtext <- c("Mean of the influenza positivity \nafter aligning at the median peak \nfor past years since 2010", 
                    "90% CI boundaries", "2017", 
                    "The boxplot visualizes spread \nof historical peaks")
    legLoc <- c(0,0)
    if (median(peakVec)>26){
        legLoc <- c(.6,0)
    }
    legend("topright", inset=legLoc, legendtext, bty="n", lty= c("solid", "dashed", "solid", NA), lwd=c(1,1,2,NA), cex=.8)
                                        #dev.off ()
    
}

if (FALSE) {

    ## setwd("~/Dropbox/git/iiag")
    rm(list=ls(all=TRUE))
    source("src/R/riley_funcs.r")
    source("src/R/haw_funcs.r")

    ## Assumes running in src/R. Change datadir as needed. Select the
    ## syndromic data after loading
    tmp <- load.iiag.data(datadir="data")
    df <- tmp$synd
    zones <- read.csv("data/country_list_ISO.csv")

    ## Plot country averages
    ## This line doesn't work, DH and SR working on the functions above at the moment
    ave_seas_plots(df,zones)
    
}
