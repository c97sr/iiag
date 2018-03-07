average_curve <- function(
                           df, # dataframe with output from other routines
                           zones, # dataframe defining transmission zones
                           plotCategory = "country", #"zone" or "country" - add whoregion?
                           nameArea = "FRA", #Select country (ISO3)
                           YearStart = 2016) { # First year of season, irresepctive of hemisphere

    require(gdata)
    require(reshape2)
    require(plyr)
    require(stringr)
    #require(ggplot2)
    #require(zoo)
    #require(epitools)
    #require(Deducer)
    

    ## Zones needed below here
    if (plotCategory=="country"){
        ISO3ind <- nameArea #Select country
        HEMIS <- zones$hemis[zones$ISO3==ISO3ind]
    } else {
        ITZind <- nameArea
        countries <- zones$ISO3[zones$itz==ITZind]
        #countries <- list(unique(countries))#levels(countries)
        ISO3ind <- nameArea #Treat as single country in rest of code
        HEMIS <- zones$hemis[zones$itz==nameArea]
        HEMIS <- HEMIS[1]
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

    ## measures <- c("SPEC_PROCESSED_NB", "DECTECTED_NB")
    ## res <- res[which(res$MEASURE_CODE %in% measures), ]
    df <- df[which(df$MEASURE_CODE=="SPEC_PROCESSED_NB" |
                   df$VIRUSTYPE_CODE=="ALL_INF"), ]

    if (plotCategory=="country"){
        df2 <- df[df$ISO3==ISO3ind, ]
        countryName <- zones$country[which(zones$ISO3==ISO3ind)]
        countryName <- countryName[1]
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
    
    ################################
    
    #Plot:
    col1 <- "darkred"#(.667,.224,.224)
    col2 <- "black"#rgb(0,0,0)#(.831,.416,.416)
    col3 <- rgb(1,1,1)#(1,.667,.667)
    endChop <- 5
    xmax <- length(currentYear)
    ymax <- min( (c( max(currentYear, na.rm=TRUE), max(pastAv, na.rm=TRUE))+20), 100, na.rm=TRUE)
    tvec <- seq((chop+1), (xmax-endChop))
    plot(pastAv[tvec], type="l",xlim=c(1, length(tvec)), ylim= c(0, ymax), ylab="Influenza positivity (%)",
         xlab="Weeks", main=paste(countryName, NOW), bty="l", cex.main=1, cex=2, cex.lab=1, cex.axis=1,#, xaxt="n", yaxs='i', xaxs='i')
         lwd=2)
    grid(NULL, NULL, col = "lightgray", lty = "dotted",
         lwd = par("lwd"), equilogs = TRUE)
    points(currentYear[tvec], type="l", lwd=3, col=col1)
    #Error bars:
    errorP <- pastAv+1.645*pastSd/nVals
    errorM <- pastAv-1.645*pastSd/nVals
    points(errorP[tvec], type="l", lty = "dashed", lwd=1.5)
    points(errorM[tvec], type="l", lty = "dashed", lwd=1.5)
    
    
    # From here - DHy1 (past average) - correct?
    boxplot(peakVec, horizontal=T, add=T, at=(ymax-10), axes=F, boxwex=6, pch=19, col=col3, flatten=1)
    
    numYears <- length(peakVec)
    pastYears <- seq(1, (numYears-1))
    yVals <- jitter(rep((ymax-5),length(peakVec)), 2)
    points(peakVec[pastYears], yVals[pastYears], pch=4, cex=1, lw=2, col=col2)#pch=18
    points(peakVec[numYears], max(yVals)+2, pch=4, cex=1, lw=2, col=col1)
    
    legendtext <- c(NOW, "Mean past positivity \n(alligned)", 
                    "90% CI", "Peak times \n(shifted)", #Past peaks
                    "Boxplot visualizes \nhistorical peaks \n")
    legLoc <- c(0,0)
    if (min(peakVec)>22){#(median(peakVec)>26){ #This condition can be changed - it's visual
      legLoc <- c(.7,0)
    }
    legend("topright", inset=legLoc, legendtext, lty= c("solid", "solid", "dashed", NA, NA), lwd=c(3,2,1.5,2,0), cex=.8,
           pch=c(NA,NA,NA,4,NA), col=c("darkred", "black", "black", "black", NA))
    #dev.off ()
    
}

if (FALSE) {

    ## You may need to change this for where you have downloaded the git repository 
    ## setwd("~/Dropbox/git/iiag")
    rm(list=ls(all=TRUE))
    source("src/R/riley_funcs.r")
    source("src/R/haw_funcs.r")

    ## Assumes running in src/R. Change datadir as needed. Select the
    ## syndromic data after loading
    tmp <- load.iiag.data(datadir="data")
    df <- tmp$lab#synd
    zones <- read.csv("data/CountryList.csv")
    zones <- data.frame(lapply(zones, as.character), stringsAsFactors=FALSE)
    colnames(zones) <- c("ISO3", "country", "itz", "whoreg", "hemis", "sov")
    zones$hemis[zones$hemis=="Northern hemishere"] <- "northern"
    zones$hemis[zones$hemis=="Southern hemisphere"] <- "southern"

    ## Plot country averages
    average_curve(df,zones)
    
}

