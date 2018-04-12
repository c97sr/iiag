average_curve <- function(df, zones, plotCategory, nameArea, YearStart){

    require(gdata)
    require(reshape2)
    require(plyr)
    require(stringr)
    require(zoo)
    require(ggplot2)
    #require(epitools)
    #require(Deducer)
    
    #### Extract desired country/zone from complete data. Retain hemosphere for season timing.
  
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
        HEMIS <- HEMIS[2]
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
    
    #### Create data frame with columns year/week/code/value, in chronological order.
  
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
    
    #### If using ITZ, aggregate over countries.
    
    if (plotCategory=="zone"){
        ## df2 < -aggregate(df2, by=list(year=df2$year, week=df2$week, code=df2$code), mean, na.rm=TRUE)
        dfx <- ddply(df2, c("year", "week", "code"), fun=mean)
    }
    
    #### Reshape data so there is one row per week.
    
    df2 <- reshape(df2, idvar=c("year", "week"), timevar=c("code"), direction="wide")
                                        #Should land in chronological order ********
    
    #### Use perventage positive (3 point moving average) as measure of incidence:
    
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
    
    #### Reshape data so there is one column per week.
    #### This is like the Excel pivot table in Julia's original code.
                                        #
    df3 <- data.frame(df2$syear, df2$sweek, df2$percPos)
    colnames(df3) <- c("syear", "sweek", "mov")
                                        #df3 <- ddply(df3, "syear")
    df3 <- reshape(df3, idvar="syear", timevar="sweek", direction="wide")
                                        #
    df3[df3=="NaN"]=NA #Some "NaN" where "NA" desired - discovered in Fiji
    
    #### Organise columns into chronological order and allign peaks.
                                        #
    rownames(df3) <- df3[, 1]#Row names as required
    df3 <- df3[, -1]
    nWeeks <- ncol(df3)
    if (HEMIS=="northern"){
        df3 <- df3[, c(seq(40,nWeeks), seq(1,39))] #Generalise, just in case ********
    }
    peakVal <- apply(df3, 1, max, na.rm=TRUE)
    peakInd <- max.col(df3[, seq(1,nWeeks)], "first")
    latest <- max(peakInd)
    
    #### Deal with week 53 (week 14 of season).
    
    w14 <- which(is.na(df3[, 14]))
    df3[w14,seq(14,52)] <- df3[w14, seq(15,52)]
    df3[w14,53]=NA
    
    #### Peaks are initially alligned to latest value, which can be stored more intuitively in a table.
    #### Peaks are shifted to median value in plot.
    
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
    
    #### Extract current year and compute past average.
    
    yearList <- as.numeric(rownames(df4))
    currentYear <- which(yearList==YearStart)
    pastAv <- as.numeric(apply(df4[seq(2, (currentYear-1)), ], 2, mean, na.rm=TRUE))
    pastSd <- as.numeric(apply(df4[seq(2, (currentYear-1)), ], 2, sd, na.rm=TRUE))
    nVals <- apply(df4, 2, function(x) length(which(!is.na(x))))
    currentYear <- as.numeric(df4[currentYear, ])
    
    ################################
    
    #Plot:
    
    #### Lots of this code just tweaks the plot asthetically.
    
    col1 <- "darkred"#(.667,.224,.224)
    col2 <- "black"#rgb(0,0,0)#(.831,.416,.416)
    col3 <- rgb(1,1,1)#(1,.667,.667)
    endChop <- 5
    xmax <- length(currentYear)
    ymax <- min( (c( max(currentYear, na.rm=TRUE), max(pastAv, na.rm=TRUE))+20), 100, na.rm=TRUE)
    tvec <- seq((chop+1), (xmax-endChop))
    if (HEMIS=="northern"){
     xvec <- c(seq(40,52),seq(1,39))
    } else {
      xvec <- seq(1,52)
    }
    lxvec <- length(xvec)
    plot(pastAv[tvec], type="l",xlim=c(1, length(tvec)), ylim= c(0, ymax), ylab="Influenza positivity (%)",
         xlab="Weeks", main=paste(countryName, NOW), bty="l", cex.main=1, cex=2, cex.lab=1, cex.axis=1,#, xaxt="n", yaxs='i', xaxs='i')
         lwd=2, xaxt="n")
    axis(1, at=seq(1, lxvec, 10), labels=xvec[seq(1, lxvec, 10)])
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
                    "90% CI", "Peak times", #Past peaks
                    "Boxplot visualizes \nhistorical peaks \n")
    legLoc <- c(0,0)
    if (min(peakVec)>22){#(median(peakVec)>26){ #This condition can be changed - it's visual
      legLoc <- c(.7,0)
    }
    legend("topright", inset=legLoc, legendtext, lty= c("solid", "solid", "dashed", NA, NA), lwd=c(3,2,1.5,2,0), cex=.8,
           pch=c(NA,NA,NA,4,NA), col=c("darkred", "black", "black", "black", NA))
    #dev.off ()
    
}

################################################################################################################################

#Actually subtype, not strain
prepStrain <- function(df, dfID, zones, proportion){ #df=dfNet
  
  library(reshape2)
  library(ggplot2)
  library(plyr)
  
  nh <- zones$ISO3[zones$hemis=="northern"]
  sh <- zones$ISO3[zones$hemis=="southern"]
  #tr <- zones$ISO3[zones$hemis=="tropics"]
  
  dfx <- data.frame(df$ISO3, df$ISO_YEAR, df$ISO_WEEK, df$VIRUS_SUB_CODE, df$ValueNumeric)
  colnames(dfx) <- c("code", "year", "week", "subtype", "value") #"Year" of season start
  
  ################################
  
  #### If incorporating ILI, writes a csv of time-series (without subtyping):
  if (proportion==FALSE){
    dfid <- dfID[which(dfID$MEASURE_CODE=="ILI_CASES" & dfID$AGEGROUP_CODE=="All"), ]
    dfid <- data.frame(dfid$ISO3, dfid$ISO_YEAR, dfid$ISO_WEEK, dfid$ValueNumeric)
    colnames(dfid) <- c("code", "year", "week", "ili")
    #Combines with dfx here:
    dfx <- merge(dfx, dfid, by=c("code", "year", "week"))
    dfx$ili[is.na(dfx$value)] <- 0
    dfx$ili[is.na(dfx$ili)] <- 0
    dfx$value <- dfx$value*dfx$ili
    dfx$ili <- NULL
    
    #ILI time series additional output:
    dfn <- dfid[dfid$code %in% nh, ]
    dfs <- dfid[dfid$code %in% sh, ]
    lateSeason <- which(as.numeric(dfn$week) <= 39)
    dfn$year[lateSeason] <- dfn$year[lateSeason]-1
    dfid <- rbind.data.frame(dfn, dfs)
    dfid$week <- NULL
    
    dfid <- ddply(dfid, c("code", "year"), function(dfid)sum(dfid$ili))
    colnames(dfid)[3] <- "ili"
    dfid <- reshape(dfid, idvar="code", timevar=c("year"), direction="wide")
    write.csv(dfid, "timeSerisILI.csv")
  }
  
  
  ################################
  
  dfn <- dfx[dfx$code %in% nh, ]
  dfs <- dfx[dfx$code %in% sh, ]
  #dft <- dfx[dfx$code %in% tr, ]
  
  lateSeason <- which(as.numeric(dfn$week) <= 39)
  dfn$year[lateSeason] <- dfn$year[lateSeason]-1
  
  subtypes <- levels(dfn$subtype)
  subtypes <- subtypes[!subtypes %in% c("", "TOTAL")]
  
  dfn <- ddply(dfn, c("code", "year", "subtype"), function(dfn)sum(dfn$value))
  dfs <- ddply(dfs, c("code", "year", "subtype"), function(dfs)sum(dfs$value))
  #dft <- ddply(dft, c("code", "year", "subtype"), function(dft)sum(dft$value))
  
  dfout <- rbind.data.frame(dfn, dfs)#, dft)
  dfout <- reshape(dfout, idvar=c("code", "year"), timevar=c("subtype"), direction="wide")
  dfout <- dfout[, c("code", "year", "V1.AH1", "V1.AH1N12009", "V1.AH3", "V1.AH5", 
                     "V1.ANOTSUBTYPED", "V1.AOTHER_SUBTYPE", 
                     "V1.BNOTDETERMINED", "V1.BVICTORIA", "V1.BYAMAGATA", 
                     "V1.AH7N9")]
  #To date, 2x H1N1, 0x non-subtypable/other subtype ********
  assign("dfout", dfout, .GlobalEnv) #Better as function output?
  write.csv(dfout, "strainYear.csv")
}

################################################################################################################################

extractCountry <- function(dfall, iso3, splitOther){
  
country <- dfall[dfall$code==iso3, ]
rownames(country) <- country$year
country <- country[, c(seq(3,12))]#(4,12),15)]
colnames(country) <- c("H1", "H1_2009", "H3", "H5", "A_NA", "A_other", "B_NA", "B_Vic", "B_Yam", "H7")#DH: matlab version excludes H7
country[is.na(country)] <- 0
#Add H5 & H7 to A other:
country$A_other <- country$A_other+country$H5+country$H7
#Add pandemic to H1:
country$H1 <- country$H1+country$H1_2009
country <- country[, c("H1", "H3", "A_other", "A_NA", "B_Vic", "B_Yam", "B_NA")]
#country$H1, country$H3, country$A_other, country$A_NA, country$B_Vic, country$B_Yam, country$B_other)]

if (splitOther==TRUE){
  country$H1 <- country$H1+.5*country$A_other
  country$H3 <-country$H3+.5*country$A_other
  country <- data.frame(country$H1, country$H3, country$A_NA, country$B_Vic, country$B_Yam, country$B_other) #Just delete column!
  palette <- c("#a6cee3",
               "#1f78b4",
               "#b2df8a",
               "#33a02c",
               "#fb9a99",
               "#e31a1c")
} else {
  palette <- c("#a6cee3",
               "#1f78b4",
               "#b2df8a",
               "#33a02c",
               "#fb9a99",
               "#e31a1c",
               "#fdbf6f")
}
nVirus <- ncol(country)

dat <- data.frame(t(country))
dat[is.na(dat)] <- 0
if (proportion==TRUE){
  dat <- sweep(dat, 2, colSums(dat), FUN="/")
}
dat$row <- seq_len(nrow(dat))
dat2 <- melt(dat, id.vars = "row")
dat2$variable <- substr(dat2$variable, 2, 5)
rnames <- rownames(dat)

dat2$row <- rnames[as.numeric(dat2$row)]

p <- ggplot(dat2, aes(x=variable, y=value, fill=as.factor(row))) + 
  geom_bar(stat="identity") +
  xlab("\nYear of season start") +
  ylab("Proportion of samples\n") +
  scale_fill_manual(values=palette)
p <- p + guides(fill=guide_legend(title="Subtype")) + ggtitle(paste(iso3, "historic profile", sep = " ", collapse = NULL))
p
}

################################################################################################################################

predictFlu <- function(dfNet, dfID, incidenceMeasure, zones, nameArea, plotACF, 
                       numWeeks, thisFar, recentData, 
                       arima1, arima2, arima3, 
                       seasonal1, seasonal2, seasonal3, seasonPeriod
){#, plotCategory, YearStart){
  
  options(warn=-1)
  #Country only for now:
  cName <- nameArea
  countryName <- zones$country[which(zones$ISO3==cName)]
  countryName <- countryName[1]
  
  if (incidenceMeasure=="sero"){
    dfx <- dfNet
    fluid1c <- dfx[dfx$ISO3==cName, ]
  }else{
    dfx <- dfID
    fluid1c <- dfx[dfx$ISO3==cName, ]
    fluid1c <- fluid1c[fluid1c$AGEGROUP_CODE=="All", ]
  }
  #
  if (incidenceMeasure=="sero"){
    df <- dfx[which(dfx$MEASURE_CODE=="SPEC_PROCESSED_NB" |
                      dfx$VIRUSTYPE_CODE=="ALL_INF"), ]
    df2 <- df[df$ISO3==cName, ]
    df2 <- data.frame(df2$ISO_YEAR, df2$ISO_WEEK, df2$MEASURE_CODE, df2$ValueNumeric)
    colnames(df2) <- c("year", "week", "code", "value")
    #### Reshape data so there is one row per week.
    df2 <- reshape(df2, idvar=c("year", "week"), timevar=c("code"), direction="wide")
    #Should land in chronological order ********
    #### Use perventage positive (3 point moving average) as measure of incidence:
    df2$percPos <- df2[, 3]/df2[, 4]*100
    #df2$percPosMov <- rollapply(df2$percPos, 3, mean, align="center", fill=T)
    df2 <- df2[, c(1,2,5)]
    names(df2)[1:3] <- paste(c("year", "week", "value"))
    fluid1is <- aggregate(df2$value, by=list(df2$year, df2$week), FUN=sum, na.rm=TRUE)
    names(fluid1is)[1:3] <- paste(c("year", "week", "value"))
    #
  }else{
    #Extract cases:
    toMatch <- c("ILI_CASES", "ILI_OUTPATIENTS", "SARIC_CASES", "SARI_DEATHS$", "SARI_INPATIENTS", "SARI_DEATHS_INPATIENTS")#All usable
    #toMatch <- c("ILI_CASES", "SARI _CASES")
    #toMatch <- c("DEATHS", "PATIENTS$")
    fluidc <- fluid1c[grep(paste(toMatch, collapse="|"), fluid1c$MEASURE_CODE),]
    #
    #fluidc <- data.frame(fluidc$YEAR_CODE,fluidc$WEEK_CODE,fluidc$MEASURE_CODE,fluidc$AGEGROUP_CODE,fluidc$DataValue)
    #Why doesn't work? - get "fluidc." in column names
    keeps <- c("ISO_YEAR","ISO_WEEK","MEASURE_CODE","AGEGROUP_CODE","ValueNumeric")
    fluidc <- fluidc[keeps]
    names(fluidc)[1:5] <- paste(c("year", "week", "case", "age", "value"))#case
    #
    fluidc$week <- as.numeric(fluidc$week)
    fluidc$year <- as.numeric(fluidc$year)
    fluidc$value <- as.numeric(as.character(fluidc$value))
    #fluisc <- fluidc[order(fluidc$week),]
    #fluisc <- order(fluidc$year),]#Retains sub-order from previous? #Col numbers retain old order
    #
    #Separate out ILI and SARI CASES only:
    cases <- c("ILI_CASES", "SARI_CASES")#Include sub-cases?
    fluidis <- fluidc[grep(paste(cases, collapse="|"), fluidc$case),]#"ILI/SARI"
    #ILI+SARI:
    fluid1is <- aggregate(fluidis$value, by=list(fluidis$year, fluidis$week), FUN=sum, na.rm=TRUE)
    names(fluid1is)[1:3] <- paste(c("year", "week", "value"))
  }
  #
  ##From here, the term "fluid" appears in variable names, even if data is fluNet. 
  ##This is because code was initially written with just fluID in mind.
  #
  #Order correctly:
  fluid1is <- fluid1is[order(fluid1is$week), ]
  fluid1is <- fluid1is[order(fluid1is$year), ]
  if (plotACF==TRUE){
    ##ACF and PACF in same figure:
    forACF <- forACF <- fluid1is$value
    par(mfrow=c(2,1))
    acf(forACF, main="Data: ACF")#, main="Autocorrelation function")
    #grid(nx = NULL, ny = nx, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
    pacf(forACF, main="Data: PACF")#, main="Partial autocorrelation function")
    #grid(nx = NULL, ny = nx, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
  }
  #
  #Disguard most recent:
  lf <- length(fluid1is$value)
  if(lf < recentData) stop("Not enough data: predict shorter time period")
  realInd <- seq(lf-recentData+1, lf)
  REAL <- fluid1is[realInd, ]
  fluid1is <- fluid1is[-realInd, ]
  #
  disguardStart <- lf-recentData-numWeeks
  if (disguardStart>0){
    fluid1is <- fluid1is[-c(1:disguardStart), ]
  }
  lf <- length(fluid1is$value)
  #
  #ARIMA:
  fluid.fit <- arima(fluid1is$value, order=c(arima1,arima2,arima3), seasonal=list(order=c(seasonal1,seasonal2,seasonal3), period=seasonPeriod), include.mean=FALSE)
  #print(fluid.fit)
  if (plotACF==TRUE){
    ##ACF and PACF in same figure:
    forACF <- fluid.fit$residuals
    par(mfrow=c(2,1))
    acf(forACF, main="Model: ACF")#, main="Autocorrelation function")
    #grid(nx = NULL, ny = nx, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
    pacf(forACF, main="Model: PACF")#, main="Partial autocorrelation function")
    #grid(nx = NULL, ny = nx, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
  }
  
  #Predictions:
  fluid.pred <- predict(fluid.fit, n.ahead=thisFar)
  xOld <- seq(1, lf)
  xReal <- seq(lf+1, lf+length(REAL$value))
  xNew <- seq(lf+1, lf+thisFar)
  y1 <- fluid.pred$pred-2*fluid.pred$se
  y2 <- fluid.pred$pred+2*fluid.pred$se
  xMax <- max(max(xOld), max(xReal))
  if (incidenceMeasure=="sero"){
    yMax <- 100
    yLabel <- "Percentage positive"
  }else{
    yMax <- max(max(fluid1is$value), max(y2))
    yLabel <- "Cases"
  }
  redTrans <- rgb(255, 0, 0255, max = 255, alpha = 50)
  
  par(mfrow=c(1,1))
  plot(xOld, fluid1is$value, type="l", xlim=c(0, xMax), ylim=c(0, yMax), xlab="Time (weeks)", ylab=yLabel, main=countryName, lwd=2)
  polygon(c(xNew,rev(xNew)),c(y2,rev(y1)),col=redTrans)
  lines(xNew, y2, col="red", lwd=2)
  lines(xNew, y1, col="red", lwd=2)
  lines(fluid.pred$pred, col="blue", lwd=2)
  lines(xReal, REAL$value, col="black", lty = "dotted", lwd=2)
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
}

################################################################################################################################

if (FALSE) {
  
  ## You may need to change this for where you have downloaded the git repository:
  ## Set working directory:
  #setwd("~/Dropbox/git/iiag")
  ## Clear workspace:
  rm(list=ls(all=TRUE))
  ## Source necessary scripts:
  source("src/R/riley_funcs.r")
  source("src/R/haw_funcs.r")
  
  ## You will need these packages - install if necessary:
  #install.packages("gdata")
  #install.packages("reshape2")
  #install.packages("plyr")
  #install.packages("stringr")
  #install.packages("zoo")
  #install.packages("ggplot2")
  
  ## Assumes running in src/R. Change datadir as needed, otherwise leave this bit:
  tmp <- load.iiag.data(datadir="data")
  dfNet <- tmp$lab
  dfID <- tmp$synd
  zones <- read.csv("data/CountryList.csv")
  zones <- data.frame(lapply(zones, as.character), stringsAsFactors=FALSE)
  colnames(zones) <- c("ISO3", "country", "itz", "whoreg", "hemis", "sov")
  zones$hemis[zones$hemis=="Northern hemishere"] <- "northern"
  zones$hemis[zones$hemis=="Southern hemisphere"] <- "southern"
  
  ################################
  # AVERAGE CURVES:
  # INPUTS - if above has run once, just run from here.
  
  ## For average curve:
  plotCategory <- "zone" #"zone" (ITZ) or "country" - add whoregion?
  nameArea <- "Southern Asia" #Select country (ISO3) or ITZ (ITZ name as in table)
  YearStart = 2017 #First year of required season
  ## Plot desired average curve:
  average_curve(dfNet, zones, plotCategory, nameArea, YearStart)
  
  
  #################################
  # SUBTYPES BY YEAR/SEASON:
  
  ## Option to plot distribution of subtypes by year:
  ## Comment out next 5 lines if no bar plot desired.
  proportion <- TRUE #False will cross-multiple with ILI incidence, true considers only "percentage positive"
  iso3 <- "USA" #Must be a country (ISO3) for now, not a zone
  splitOther <- FALSE #Divide non-H1/H3 between these 2 (for "tidier" plots)
  prepStrain(dfNet, dfID, zones, proportion)
  extractCountry(dfout, iso3, splitOther)
  
  ################################
  # SARIMA:
  # INPUTS - if above has run once, just run from here.
  
  ## For SARIMA:
  #plotCategory <- "country" #"zone" (ITZ) or "country" - add whoregion?
  nameArea <- "USA" #Select country (ISO3) or ITZ (ITZ name as in table)
  numYears <- 4
  numWeeks <- 52*numYears #Total number of weeks to use
  recentData <- 20 #Remove most recent weeks
  #Predict this many weeks:
  thisFar <- 28
  incidenceMeasure <- "case"#"sero" or "case"
  #
  #(S)ARIMA parameters:
  arima1 <- 1 #AR
  arima2 <- 1 #Diff
  arima3 <- 0 #MA
  seasonal1 <- 0 #Seasonal AR
  seasonal2 <- 1 #Seasonal diff
  seasonal3 <- 1 #Seasonal MA
  seasonPeriod <- 52 #Weeks - 52, 26 for tropics?
  #
  plotACF <- TRUE
  #Run the function:
  predictFlu(dfNet, dfID, incidenceMeasure, zones, nameArea, plotACF, 
             numWeeks, thisFar, recentData, 
             arima1, arima2, arima3, 
             seasonal1, seasonal2, seasonal3, seasonPeriod)
  ##Message from Sarah Kramer:
  #The 18 countries I considered as "tropical" (and that had decent-quality syndromic and virologic data) were:
  #Brazil, Bolivia, Colombia, Cuba, Ecuador, Honduras, Paraguay, Peru, Kenya, Madagascar, 
  #Oman, Pakistan, Bangladesh, Bhutan, Cambodia, Indonesia, Singapore, Thailand 
}