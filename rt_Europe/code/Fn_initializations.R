## Script for initializing some parameters for the functions
## function to assign data (obs), find the first and last dates of a season
##  More information in 'Comparison of filtering methods for the modeling and retrospective forecasting of influenza epidemics' (PLoS Compute Biol)
##  by Wan Yang, Alicia Karspeck, and Jeffrey Shaman, 2014

Fn_dates=function(season){
  if (season=="2009-10"){
    weeks <- 1:75 # 2009 has 53 weeks
    start_date <- as.Date('2009/5/2') # end of week 18 # Wan used 2009/10/10
    end_date <- as.Date('2010/10/2') # end of week 39
    nsn <- 75
  } else if (season=="2010-11"){
    weeks <- 76:127
    start_date <- as.Date('2010/10/9') # end of week 40
    end_date <- as.Date('2011/10/1') # end of week 39
    nsn <- 52
  } else if (season=="2011-12"){
    weeks <- 128:179
    start_date <- as.Date('2011/10/8')
    end_date <- as.Date('2012/9/29')
    nsn <- 52
  } else if (season=='2012-13'){
    weeks <- 180:231
    start_date <- as.Date('2012/10/6')
    end_date <- as.Date('2013/9/28')
    nsn <- 52
  } else if (season=='2013-14') {
    weeks <- 232:283
    start_date <- as.Date('2013/10/5')
    end_date <- as.Date('2014/9/27')
    nsn <- 52
  } else if (season=='2014-15') {
    weeks <- 284:335
    start_date <- as.Date('2014/10/4')
    end_date <- as.Date('2015/9/26')
    nsn <- 52
  } else if (season=='2015-16') {
    weeks <- 336:388 # 2015 has 53 weeks
    start_date <- as.Date('2015/10/3')
    end_date <- as.Date('2016/10/1')
    nsn <- 53
  } else if (season=='2016-17') {
    weeks <- 389:440
    start_date <- as.Date('2016/10/8')
    end_date <- as.Date('2017/9/30')
    nsn <- 52
  } else if (season=='2017-18') {
    weeks <- 441:492
    start_date <- as.Date('2017/10/7')
    end_date <- as.Date('2018/9/29')
    nsn <- 52
  } else if (season=='2018-19') {
    weeks <- 493:544
    start_date <- as.Date('2018/10/6')
    end_date <- as.Date('2019/9/28')
    nsn <- 52
  } else {
    weeks=NA
  }
  rec <- list(weeks=weeks, start_date=start_date, end_date=end_date, nsn=nsn)
  #rec <- weeks
  rec;
}

# function to calculation the mode
MODE <- function(x) {
  md=count(x)
  mode=md$x[which.max(md$freq)];
  count=max(md$freq);
  c(mode,count)
}