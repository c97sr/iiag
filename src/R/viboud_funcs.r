image.scale <- function(z, zlim, col = heat.colors(12),
breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
 if(!missing(breaks)){
  if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
 }
 if(missing(breaks) & !missing(zlim)){
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
 }
 if(missing(breaks) & missing(zlim)){
  zlim <- range(z, na.rm=TRUE)
  zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
  zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
 }
 poly <- vector(mode="list", length(col))
 for(i in seq(poly)){
  poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
 }
 xaxt <- ifelse(horiz, "s", "n")
 yaxt <- ifelse(horiz, "n", "s")
 if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
 if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
 if(missing(xlim)) xlim=XLIM
 if(missing(ylim)) ylim=YLIM
 plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
 for(i in seq(poly)){
  if(horiz){
   polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
  }
  if(!horiz){
   polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
  }
 }
}

if (FALSE) {

    rm(list=ls(all=TRUE))
    ## setwd() 
    source("src/R/riley_funcs.r")
    source("src/R/viboud_funcs.r")
    dataflu=load.iiag.data()

    ## Extract data for all countries
    x <- extract.incidence(
        dataflu$synd,
        minYear=2010,
        sel_iso3 = unique(dataflu$synd$ISO3),
        sel_ag = c("All"),
        sel_measure = c("ILI_CASES"))

    ## Identify countries with robust ILI
    ## Sum of ILI cases by country
    colSums(is.na(x))

    ## Proportion of missing weeks by country
    colSums(is.na(x))/dim(x)[1]

    ## List of countries with more than 50% of weeks of complete data
    sel_iso=names(which(colSums(is.na(x))/dim(x)[1]<.5))

    ## Reextract data for set of more complete countries
    x <- extract.incidence(
        dataflu$synd,
        minYear=2010,
        sel_iso3 = sel_iso,
        sel_ag = c("All"),
        sel_measure = c("ILI_CASES")
    )

    ## fluIliCountryData <- x
    ## save(fluIliCountryData,file="~/tmp/fluIliCountryData.rda")
    
    ## Quick diagnostic plot of the incidence
    ## numeric index for date
    datee=seq(from=2010,by=1/52.17, length.out=dim(x)[1])
    
    plot(datee, scale(x[,"USA"]),type="l",col="red",
         ylim=c(min(scale(x),na.rm=TRUE),max(scale(x),na.rm=TRUE)),
         ylab="ILI", xlab="Week")
    points(datee,scale(x[,"NLD"]),type="l",col="blue")
    points(datee,scale(x[,"IRL"]),type="l",col="green")
    points(datee,scale(x[,"NZL"]),type="l",col="purple")
    legend(x="topright",legend=c("USA","NLD","IRL","NZL"),
           col=c("red","blue","green","purple"),lwd=2)
    
    ## Getting basic geographic country info 
    ## All world countries
    ## setwd("H://FluNet/Data/")
    countrydesc=read.table("data/country_list_ISO.csv", sep=',',
                           header=TRUE)

    ## countries that have robust ILI data in x
    isos=colnames(x)
    x2=x[,order(isos)] ## sorting countries in ILI matrix by ISO3  for later use

    countries=countrydesc[ (countrydesc$ISO3 %in% isos), ]
    countries$COUNTRY=as.character(countries$Country)
    nco=dim(countries)[1]
    countries=countries[order(countries$ISO3),]

    idx.lat=order(countries$Latitude)

    ## Time series heatmap
    nwk=dim(x)[1]
    x2.sc=apply(x2, 2, scale)
    colnames(x2.sc)

    ## plot
    nbreaks=11 ##nb of color breaks; seems to work for most flu examples
    ## breaks based on quantiles
    breaksq=c(quantile(x2.sc,p=seq(0,1,length.out=nbreaks+1),na.rm=T)) 
    ## breaks defined by user --here putting more weight on high flu incidences
    breaksu=c(min(x2.sc,na.rm=T),-1.5,-1,-.5, 0, .5, 1,1.5, 2,4,6,10)  

    layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(4))
    layout.show(2)
    par(mar=c(4,6.5,1,1))
    image(x=datee, z=as.matrix(x2.sc[,idx.lat]),
          col=heat.colors(nbreaks)[nbreaks:1], 
          breaks=breaksu,yaxt="n",xlab="Weeks")

    axis(2, at=seq(0,1,,dim(x2.sc)[2]),
         labels=colnames(x2.sc)[idx.lat],las=1,cex.axis=0.5)
    abline(v=seq(2010,2018),lty=1,col="grey")

    par(mar=c(5,0,1,5))
    image.scale(x2.sc, col=heat.colors(nbreaks)[nbreaks:1], breaks=breaksu,
                horiz=FALSE, yaxt="n", xaxt="n", xlab="", ylab="")
    axis(4)
    mtext("Scaled ILI cases", side=4, line=2)
    box()

}
