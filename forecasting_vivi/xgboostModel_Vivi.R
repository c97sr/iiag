
#' Make a note of when the report was generated.
Sys.time()

#' As always, remove all objects fromt the workspace before starting.
rm(list = ls(all = TRUE))


# library('idd') # WHO influenza-likeness-illness data; idd package is not available (for R version 3.6.1)
library('ggplot2') # transform integer numbers to categorical variable
library('xgboost') # boosted regression tree

library('formattable') # formatting on data frames
library('gtools') # permutations
library("lattice") # heat plot
library("RColorBrewer") # colour platter
library("raster") # colour platter
library("rasterVis") # colour platter
# library("purrr") # Partial dependence plots
library("tidyverse")
library("dplyr") # mutate
library("lubridate")
library("mlr")
library("Matrix")
library("e1071")
library("DiagrammeR") # tree plot
library("rsvg")
library("DiagrammeRsvg") # multi tree plot
library("data.table") # convert data frame to data table
library("aweek") # convert discrete time variable into contious variable
library("knitr")
library("devtools") # to install package from github
# install_github("AppliedDataSciencePartners/xgboostExplainer")
# library("xgboostExplainer") # to interpret every single tree 
# library("pROC")
library("maps")
library("eurostat")


# load(file = "fluIliCountryData.rda")

source("C:/Users/hw3616/Desktop/Imperial/Project1_Forecasting/Project_Coding/iiag/forecasting_vivi/vivi_funcs.R")
source("C:/Users/hw3616/Desktop/Imperial/Project1_Forecasting/Project_Coding/iiag/forecasting_vivi/gbm_complex_funcs.R")


#' Load country list for WHO data.
countryISO <- read.csv("C:/Users/hw3616/Desktop/Imperial/Project1_Forecasting/Project_Coding/iiag/data_old/country_list_ISO.csv")

#' Load WHO FluID data set.
#' From 2010 week 1 to 2018 week 9 

fluWHO <- load.iiag.data.fluid(datadir = "C:/Users/hw3616/Desktop/Imperial/Project1_Forecasting/Project_Coding/iiag/data_old")

#' Extract weekly incidence for raw WHO FluID dataset
minprop <- 0.5
fluWHO.incidence <- extract.incidence.who(
  fluWHO,
  sel_iso3 = unique(fluWHO$ISO3),
  sel_ag = c("All"),
  sel_measure = c("ILI_CASES"),
  minYear=2010,
  maxYear = 2018)

sel_iso <- names(which(colSums(is.na(fluWHO.incidence))/dim(fluWHO.incidence)[1]<minprop))

#' countries whose data vailablity is over 50% are remained
fluWHO.incidence <- extract.incidence.who(
  fluWHO,
  sel_iso3 = sel_iso,
  sel_ag = c("All"),
  sel_measure = c("ILI_CASES"),
  minYear=2010,
  maxYear=2018
)

#' Check if the data of left countries is eligible for xgboost
#' check the data availablity in each year
#' exclude countries which do not have at least five consecutive weeks data in any year or do not have 
#' data records in 2010 or 2017 and 2018.
country_year <- NULL
for (i in 1:length(sel_iso)){
  tmp <- duration(sel_iso[i])
  country_year <- rbind(country_year, tmp)
}
country_year <- as.data.frame(country_year)
colnames(country_year) <- c("Country","2010","2011","2012","2013","2014","2015",
                            "2016","2017","2018","start_year","end_year")
country_year$end_year <- as.numeric(as.character(country_year$end_year))
country_year$start_year <- as.numeric(as.character(country_year$start_year))
rownames(country_year) <- c(1:nrow(country_year))

countryCode_no1718Or10 <- country_year$Country[which(country_year$end_year < 2017 | country_year$`2010`=="No"| 
                                                       country_year$`2018`=="No")]
countryCode_no10 <- country_year$Country[which(country_year$`2010`=="No")]
#' countries will not be used in xgboost model because of lack of 2017 and 2018 data or 2010 data
#' They are Barbados,Belarus,Bhutan, Honduras, New Zealand, Nigeria,Oman, Pakistan,Singapore,Tajikistan,
#' Thailand, Macedonia, the former Yugoslav Republic of, Northern Mariana Islands. 
country_no1718Or10 <- sel_iso[which(sel_iso %in%countryCode_no1718Or10)] # 20 in total

#' Exlcude the countries do not have data of 2010 or 2017 or 2018
#' 51 countries left
sel_iso_xgb <- sel_iso
for (i in 1:length(country_no1718Or10)){
  index <- which(sel_iso_xgb == country_no1718Or10[i])
  sel_iso_xgb <- sel_iso_xgb[-(index)]
}


#' Check if left countries have less than 5 weeks data in a year because I will do the 4-week ahead foreacast which
#' requires data of 5 consective weeks
#' Exclude countires whose datasets are uneligible to do the 4-week ahead forecast


# extract incidence of 41 eligible countries.
fluWHO.incidence <- extract.incidence.who(fluWHO,
                                          sel_iso_xgb,
                                          sel_ag = c("All"),
                                          sel_measure = c("ILI_CASES"),
                                          minYear = 2010,
                                          maxYear = 2018
                                          )


#### XGBooost ####

# example of data frame generated from gbm_complex() function
USA_complex1 <- gbm_complex(fluWHO.incidence, "USA", 10,1) # dataframe for 1-week ahead 
USA_complex2 <- gbm_complex(fluWHO.incidence, "USA", 10,2) # dataframe for 2-week ahead
USA_complex3 <- gbm_complex(fluWHO.incidence, "USA", 10,3) # dataframe for 3-week ahead
USA_complex4 <- gbm_complex(fluWHO.incidence, "USA", 10,4) # dataframe for 4-week ahead



#' check if data frames for 1,2,3,4-week ahead forecast are in the same size
length1 <- c()
for (i in 1:length(sel_iso_xgb)){
  a <- adjust.data.size(fluWHO.incidence, sel_iso_xgb[i], 10,1)
  tmp <- nrow(a)
  length1 <- append(length1,tmp)
  
}

length2 <- c()
for (i in 1:length(sel_iso_xgb)){
  a <- adjust.data.size(fluWHO.incidence, sel_iso_xgb[i], 10,2)
  tmp <- nrow(a)
  length2 <- append(length2,tmp)
}

length3 <- c()
for (i in 1:length(sel_iso_xgb)){
  a <- adjust.data.size(fluWHO.incidence, sel_iso_xgb[i], 10,3)
  tmp <- nrow(a)
  length3 <- append(length3,tmp)
}

length4 <- c()
for (i in 1:length(sel_iso_xgb)){
  a <- adjust.data.size(fluWHO.incidence, sel_iso_xgb[i], 10,4)
  tmp <- nrow(a)
  length4 <- append(length4,tmp)
}



# one-week ahead forecast
# 2010-2014 training, 2015 test
compare_pred15 <- compare_accuracy(sel_iso_xgb, 10, 0, 4 ,1)
length(which(compare_pred15$Accurate==1))/nrow(compare_pred15) # 0.7868296

# 2011-2015 training, 2016 test
compare_pred16 <- compare_accuracy(sel_iso_xgb,10, 1, 4, 1) 
length(which(compare_pred16$Accurate==1))/nrow(compare_pred16) # 0.7470383

# 2012-2016 training. 2017 test
compare_pred17 <- compare_accuracy(sel_iso_xgb,10, 2, 4, 1)
length(which(compare_pred17$Accurate==1))/nrow(compare_pred17)

# 2010-2016 traing, 2017 test
compare1016_pred17 <- compare_accuracy(sel_iso_xgb,10, 0, 6, 1)
length(which(compare1016_pred17$Accurate==1))/nrow(compare1016_pred17) # 0.6979038

# two-week ahead forecast
# 2010-2104 training, 2015 test
compare_pred15Two <- compare_accuracy(sel_iso_xgb, 10, 0, 4 ,2)

# 2011-2015 training, 2016 test
compare_pred16Two <- compare_accuracy(sel_iso_xgb, 10, 1, 4,2)

# 2012-2016 training. 2017&2018 test
compare_pred17Two <- compare_accuracy(sel_iso_xgb, 10, 2, 4, 2)

# 2010-2016 traing, 2017&2018 test
compare1016_pred17Two <- compare_accuracy(sel_iso_xgb, 10, 0, 6, 2)

# three week ahead
# 2010-2104 training, 2015 test
compare_pred15Three <- compare_accuracy(sel_iso_xgb, 10, 0, 4 ,3)

# 2011-2015 training, 2016 test
compare_pred16Three <- compare_accuracy(sel_iso_xgb, 10, 1, 4, 3)

# 2012-2016 training. 2017&2018 test
compare_pred17Three <- compare_accuracy(sel_iso_xgb, 10, 2, 4, 3)

# 2010-2016 traing, 2017&2018 test
compare1016_pred17Three <- compare_accuracy(sel_iso_xgb, 10, 0, 6, 3)

# four week ahead
# 2010-2104 training, 2015 test
compare_pred15Four <- compare_accuracy(sel_iso_xgb, 10, 0, 4 ,4)

# 2011-2015 training, 2016 test
compare_pred16Four <- compare_accuracy(sel_iso_xgb, 10, 1, 4, 4)

# 2012-2016 training. 2017&2018 test
compare_pred17Four <- compare_accuracy(sel_iso_xgb, 10, 2, 4, 4)

# 2010-2016 traing, 2017&2018 test
compare1016_pred17Four <- compare_accuracy(sel_iso_xgb, 10, 0, 6, 4)


# example of single tree
usa.tree.one <- xgboost.model.train(fluWHO.incidence,"USA",10,0,4,1)
usa.tr.one <- xgboost_dat(USA_complex1,2010,2014) 
usa.ts.one <- xgboost_dat(USA_complex1, 2015,2015)
usa.pred.one <- predict(usa.tree.one, usa.ts.one)
usa.xgboost.model.pred.output.one <- xgboost.model.pred.output(USA_complex1,2015,2015,usa.pred.one)

features <- attr(usa.ts.one, ".Dimnames")[[2]]
imp <- xgb.importance(features, usa.tree.one)
xgb.plot.importance(imp)


explainer <- buildExplainer(usa.tree.one,usa.tr.one, type="binary", base_score = 0.5, trees_idx = NULL)
pred.breakdown <- explainPredictions(usa.tree.one, explainer, usa.ts.one)
cat('Breakdown Complete','\n')

weights <- rowSums(pred.breakdown)
pred.xgb = 1/(1+exp(-weights))
cat(max(xgb.preds-pred.xgb),'\n')
idx_to_get = as.integer(802)
test[idx_to_get,-"left"]
showWaterfall(usa.tree.one, explainer, usa.ts.one, model.matrix(~.+0,USA_complex1[(195:242),]),1, type = "binary")



#### Heat plot for xgboost model ####

#' Print out heat chart 
# 2010-2014 training, 2015 test
for (i in 1:length(sel_iso_xgb)){
  #pdf(paste0(sel_iso_xgb[i],".pdf"))
  forecast_result <- xgboost.model.pred(fluWHO.incidence, sel_iso_xgb[i], 10, 0, 4, 1)
  frequency_table <- freq_table(forecast_result, 10)
  print(heat_plot(frequency_table, sel_iso_xgb[i]))
  #dev.off()
}

# 2011-2015 train, 2016 test
for (i in 1:length(sel_iso_xgb)){
  pdf(paste0(sel_iso_xgb[i],".pdf"))
  forecast_result <- xgboost.model.pred(fluWHO.incidence, sel_iso_xgb[i], 10, 1, 4, 1)
  frequency_table <- freq_table(forecast_result, 10)
  print(heat_plot(frequency_table, sel_iso_xgb[i]))
  dev.off()
}

# 2012-2016 train, 2017 test
for (i in 1:length(sel_iso_xgb)){
  pdf(paste0(sel_iso_xgb[i],".pdf"))
  forecast_result <-  xgboost.model.pred(fluWHO.incidence, sel_iso_xgb[i], 10, 2, 4, 1)
  frequency_table <- freq_table(forecast_result, 10)
  print(heat_plot(frequency_table, sel_iso_xgb[i]))
  dev.off()
}


#### historical average model ####

# example 
USA_hist_pred_one <- hist_average(fluWHO.incidence,"USA",10,1)
USA_hist_pred_two <- hist_average(fluWHO.incidence,"USA",10,2)
USA_hist_pred_three <- hist_average(fluWHO.incidence,"USA",10,3)
USA_hist_pred_four <- hist_average(fluWHO.incidence,"USA",10,4)


#' calculate the accuracy of historical average model
# 8050 / 11424 = 0.69
oneWeek_ahead_totalAccuracy_hist <- compare_accuracy_hist(fluWHO.incidence,sel_iso_xgb,10,1)
length(which(oneWeek_ahead_totalAccuracy_hist$Accurate==1))/nrow(oneWeek_ahead_totalAccuracy_hist)

# 0.68
twoWeek_ahead_totalAccuracy_hist <- compare_accuracy_hist(fluWHO.incidence,sel_iso_xgb,10,2)
length(which(twoWeek_ahead_totalAccuracy_hist$Accurate==1))/nrow(twoWeek_ahead_totalAccuracy_hist)

# 0.67
threeWeek_ahead_totalAccuracy_hist <- compare_accuracy_hist(fluWHO.incidence,sel_iso_xgb,10,3)
length(which(threeWeek_ahead_totalAccuracy_hist$Accurate==1))/nrow(threeWeek_ahead_totalAccuracy_hist)

# 0.66
fourWeek_ahead_totalAccuracy_hist <- compare_accuracy_hist(fluWHO.incidence,sel_iso_xgb,10,4)
length(which(fourWeek_ahead_totalAccuracy_hist$Accurate==1))/nrow(fourWeek_ahead_totalAccuracy_hist)

# plot the accuracy drop off as number of week-ahead increases
compareAccuracy_total_hist <- cbind(c(1:4),c(0.65,0.64,0.62,0.61)) %>% 
  as.data.frame()
colnames(compareAccuracy_total_hist) <- c("nWeek_ahead","percentage")

#### repeat model ####

# get the predictions by repeat model
USA_repeat_one <- repeat_model(fluWHO.incidence,'USA',10,1)
USA_repeat_two <- repeat_model(fluWHO.incidence,'USA',10,2)
USA_repeat_three <- repeat_model(fluWHO.incidence,'USA',10,3)
USA_repeat_four <- repeat_model(fluWHO.incidence,'USA',10,4)


# get the total accuract for 1,2,3,4-weeks ahead
oneWeek_ahead_totalAccuracy <- compare_accuracy_repeat(fluWHO.incidence,sel_iso_xgb,10,1)
length(which(oneWeek_ahead_totalAccuracy$Accurate==1))  # 9077
length(which(oneWeek_ahead_totalAccuracy$Accurate==1))/nrow(oneWeek_ahead_totalAccuracy) # 0.7945

twoWeek_ahead_totalAccuracy <- compare_accuracy_repeat(fluWHO.incidence,sel_iso_xgb,10,2)
length(which(twoWeek_ahead_totalAccuracy$Accurate==1)) # 8896
length(which(twoWeek_ahead_totalAccuracy$Accurate==1))/nrow(twoWeek_ahead_totalAccuracy) # 0.7440616

threeWeek_ahead_totalAccuracy <- compare_accuracy_repeat(fluWHO.incidence,sel_iso_xgb,10,3)
length(which(threeWeek_ahead_totalAccuracy$Accurate==1)) #8021
length(which(threeWeek_ahead_totalAccuracy$Accurate==1))/nrow(threeWeek_ahead_totalAccuracy) #0.7021183

fourWeek_ahead_totalAccuracy <- compare_accuracy_repeat(fluWHO.incidence,sel_iso_xgb,10,4)
length(which(fourWeek_ahead_totalAccuracy$Accurate==1)) # 7674
length(which(fourWeek_ahead_totalAccuracy$Accurate==1))/nrow(fourWeek_ahead_totalAccuracy) # 0.67

#' prepare the dataframe to plot the accuracy drop off as number of week-ahead increases
compareAccuracy_total <- cbind(c(1:4),c(10181, 9401, 8415,7996),c(0.75,0.69,0.62,0.59))
compareAccuracy_total <- as.data.frame(compareAccuracy_total)
colnames(compareAccuracy_total) <- c("nWeek_ahead","accuracy","percentage")
compareAccuracy_total$nWeek_ahead <- as.numeric(compareAccuracy_total$nWeek_ahead)

compareAccuracy_total$nWeek_ahead <- factor(compareAccuracy_total$nWeek_ahead,
                                            levels = c("one-week ahead","two-week ahead",
                                                       "three-week ahead","four-week ahead"))
compareAccuracy_total$accuracy <- as.numeric(as.character(compareAccuracy_total$accuracy))
compareAccuracy_total$percentage <- as.numeric(as.character(compareAccuracy_total$percentage))
compareAccuracy_total$percentage <- factor(compareAccuracy_total$percentage,levels = c(0,0.2,0.4,0.6,0.8,1))


#### plots ####
#' Plot world distribution of countries.
#' Extract geographic information of countries. 
#' To fill different colors in temperate and non-temperate, 
#' I will give every country a level to distinguish 
country_region <- c()
for (i in 1:length(sel_iso_xgb)){
  index <- grep(sel_iso_xgb[i],countryISO$ISO3)
  tmp <- countryISO[index,]
  country_region <- rbind(country_region, tmp)
}

for (i in 1:nrow(country_region)){
  if (country_region$Latitude[i] >= 23.25 | country_region$Latitude[i] <= -23.25){
    country_region$Region[i] <- "temperate"
  }
  if (country_region$Latitude[i] < 23.25 && country_region$Latitude[i] > -23.25){
    country_region$Region[i] <- "non-temperate"
  }
}

for (i in 1:nrow(country_region)){
  if (country_region$Latitude[i] > 0){
    country_region$Hemisphere[i] <- "Nothern hemisphere"
  }
  if (country_region$Latitude[i] < 0){
    country_region$Hemisphere[i] <- "Southern hemisphere"
  }
}
country_region$Region <- as.factor(country_region$Region)

length(which(country_region$Region == "temperate")) # 41
which(country_region$Region == "non-temperate") # 0
which(country_region$Hemisphere == "Southern hemisphere") # 0

#' change the country names in country_idd to the same as in world map
world_map <- map_data("world")
map.country.name(country_region$Country) # 26 35 40
# which(country_region$Country %in% unique(world_map$region)==FALSE) # 26 35 40

country_region$Country[c(26,35,40)] # Moldova, Republic of Russian Federation, United States 

mapCountry <- country_region$Country
country_region <- cbind(country_region$Country, mapCountry, country_region[,(2:6)])
colnames(country_region) <- c("Country","mapCountry",colnames(country_region)[3:7])
country_region$mapCountry <- as.character(country_region$mapCountry)

country_region$mapCountry[26] <- "Moldova"
country_region$mapCountry[35] <- "Russia" 
country_region$mapCountry[40] <- "USA"
country_region$mapCountry <- as.factor(country_region$mapCountry)

myCountries <- NULL
for (i in 1:nrow(country_region)){
  if((country_region$mapCountry[i] %in% unique(world_map$region)) == TRUE)
    tmp <- world_map[which(country_region$mapCountry[i] == world_map$region),]
  myCountries <- rbind(myCountries, tmp)
}

for (i in 1:nrow(myCountries)){
  if (myCountries$region[i] %in% country_region$mapCountry[which(country_region$Region == "temperate")] == TRUE){
    myCountries$Region[i] <- "temperate"
  }
  if(myCountries$region[i] %in% country_region$mapCountry[which(country_region$Region == "temperate")] == FALSE){
    myCountries$Region[i] <- "non-temperate"
  }
}


xgb_country_map <- borders("world", colour="gray50", fill="gray50")
ggplot() +  xgb_country_map +
  geom_polygon(data = myCountries, 
               aes(x = long, y = lat, group = group, fill = Region), color = "white")+
  scale_fill_manual(values = c("red","blue"))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank())

#' Visulizations of data
#' comparison between raw numeric incidence data and categorical data
# example
rawData_TS_plot(fluWHO.incidence,"ARM","Armenia")

#'save time series plot into work directory
for (i in 1:nrow(country_region)){
  pdf(paste0(country_region$Country[i],".pdf"))
  print(rawData_TS_plot(fluWHO.incidence,country_region$ISO3[i], country_region$Country[i]))
  dev.off()
}

#' visulization of predictions results against actual values
# example: USA
USA_pred15 <- xgboost.model.pred(fluWHO.incidence,"USA",10,0,4,1)
USA_pred16 <- xgboost.model.pred(fluIliCountryData,"USA",10,1,4,1)
USA_pred17 <- xgboost.model.pred(fluIliCountryData,"USA",10,2,4,1)
predTS_plot(USA_pred15)
predTS_plot(USA_pred16)
predTS_plot(USA_pred17)

#' time series plot for 6-year training and 1-year test
#' and save them to work directory
for (i in 1:length(sel_iso_xgb)){
  pdf(paste0(country_region$Country[i],".pdf"))
  forecast_result <- xgboost.model.pred(fluWHO.incidence, sel_iso_xgb[i],10,0,6,1)
  print(predTS_plot(forecast_result))
  dev.off()
}

#' plot the accuracy of historical model that the accuracy will keep roughlt stable as the number 
#' of week ahead forecast increases
# create a text
grob_hist <- grobTree(textGrob("B", x=0.8,  y=0.95, hjust=0,
                               gp=gpar(col="black", fontsize=24, fontface="bold")))

ggplot(data = compareAccuracy_total_hist, aes(x = nWeek_ahead, y = percentage, group = 1)) +
  geom_line(size=1,colour="darkgreen")+
  geom_point(size=3,colour="darkgreen")+
  scale_y_continuous("Percentage of accurate forecst", breaks = c(0.00,0.25,0.50,0.75,1), limits = c(0,1),
                     expand = c(0, 0))+
  xlab("n-week ahead")+ 
  theme_bw()+
  annotation_custom(grob_hist)

#' plot the accuracy of repeat model that the accuracy drops off as the numebr of week ahead forecasted 
#' increases

# Create a text
grob_repeat <- grobTree(textGrob("A", x=0.8,  y=0.95, hjust=0,
                                 gp=gpar(col="black", fontsize=24, fontface="bold")))

ggplot(data = compareAccuracy_total, aes(x = nWeek_ahead, y = percentage, group = 1)) +
  geom_line(size=1,colour="darkgreen")+
  geom_point(size=3,colour="darkgreen")+
  scale_y_continuous("Percentage of accurate forecst", breaks = c(0.00,0.25,0.50,0.75,1), limits = c(0,1),
                     expand = c(0, 0))+
  xlab("n-week ahead")+
  theme_bw()+
  annotation_custom(grob_repeat)

#' plot the distribution of categories
category <- NULL
for (i in 1:length(sel_iso_xgb)){
  tmp <- gbm_complex(fluWHO.incidence,sel_iso_xgb[i],10)
  tmp <- cbind(rep(sel_iso_xgb[i],nrow(tmp)),tmp)
  category <- rbind(category, tmp)
}
colnames(category) <- c("CountryCode",colnames(category)[2:6])

categorySummaryChart <- cbind(c(1:10),c(length(which(category$Y_week0==1)),length(which(category$Y_week0==2)),
                                        length(which(category$Y_week0==3)),length(which(category$Y_week0==4)),
                                        length(which(category$Y_week0==5)),length(which(category$Y_week0==6)),
                                        length(which(category$Y_week0==7)),length(which(category$Y_week0==8)),
                                        length(which(category$Y_week0==9)),length(which(category$Y_week0==10))))%>%
  as.data.frame()
colnames(categorySummaryChart) <- c("Category","Count")
categorySummaryChart$Density <- round(categorySummaryChart$Count/sum(categorySummaryChart$Count),3)

# Change the width of bars
ggplot(data=categorySummaryChart, aes(x=Category, y=Density)) +
  geom_bar(stat="identity", width=1, colour = "white", fill="steelblue")+
  scale_x_continuous("Category", breaks = c(1:10))+
  scale_y_continuous("Density",expand = c(0, 0),limits = c(0,1.00))+
  geom_text(aes(label=Density), vjust=-0.2, size=3.5)+
  # geom_text(aes(label=Density), vjust=1.6, color="white", size=3.5)+
  theme_bw()

ggplot(category, aes(x=Y_week0)) + 
  geom_histogram(aes(y=..density..), fill = "steelblue", binwidth = 1)+
  #scale_x_continuous("Category", breaks = c(1:10), limits = c(1,10))+
  # scale_y_continuous("Density",expand = c(0, 0))+
  xlab("Category")+
  ylab("Density")+
  theme_bw()

#' pick up European countries in the dataset
# geography information of Euro countries
eu_countries_region <- euro_countries(country_region)

#' ISO3 codes will be used through analysis
# ISO3 code of Euro countries
eu_xgb <- as.character(eu_countries_region$ISO3)


# incidence data of Euro countries
fluEuro.incidence <- extract.incidence.who(fluWHO,
                                           sel_iso3 = eu_xgb,
                                           sel_ag = c("All"),
                                           sel_measure = c("ILI_CASES"),
                                           minYear=2010,
                                           maxYear = 2018)
#' visualize Euro countries in a map
# Euro countries in package "maps"
# eu_xgb_counties will be only used in plotting Euro map
eu_xgb_countries <- eu_countries_region$mapCountry

eu_xgb_map <- NULL
for (i in 1:length(eu_xgb_countries)){
  if((eu_xgb_countries[i] %in% unique(world_map$region)) == TRUE)
    tmp <- world_map[which(eu_xgb_countries[i] == world_map$region),]
  eu_xgb_map <- rbind(eu_xgb_map, tmp)
  eu_xgb_map
}
eu_xgb_map$continent <- "Euro"

# Retrievethe map data
map.country.name(eu_countries$name) # 3, 28
eu_countries$name[c(3,28)] # Czechia, United Kingdom
eu_countries$mapName <- eu_countries$name
eu_countries$mapName[3] <- "UK"
eu_countries$mapName[28] <- "Czech Republic" 
eu.maps <- map_data("world", region = eu_countries$mapName)


# Compute the centroid as the mean longitude and lattitude
# Used as label coordinate for country's names
region.lab.data <- eu.maps %>%
  group_by(region) %>%
  summarise(long = mean(long), lat = mean(lat))

xgb_euro_map <- borders("world", regions = eu_countries$mapName, colour="gray50", fill="white")
ggplot() +  xgb_euro_map +
  geom_polygon(data = eu_xgb_map, 
               aes(x = long, y = lat, group = group,fill = continent), color = "white")+
  geom_text(aes(x = long, y = lat, label = region), data = region.lab.data,  size = 3, hjust = 0.5)+
  # scale_fill_manual(values = c("salmon"))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank())

#' results of Euro
#' Following is the average accuracy scores of Euro with different traing length
# 1-week ahead, 2010-2014 training, 2015 test
Euro_accuracy_15 <- compare_accuracy(eu_xgb,fluEuro.incidence,10,0,4,1)

# 1-week ahead, 2011-2015 training, 2016 test
Euro_accuracy_16 <- compare_accuracy(eu_xgb,fluEuro.incidence,10,1,4,1)

# 1-week ahead, 2012-2016 training, 2017 & 2018 test
Euro_accuracy_1718 <- compare_accuracy(eu_xgb,fluEuro.incidence,10,2,4,1)

# 1-week ahead, 2010-2016 training, 2016 test
Euro_accuracy_1016_1718 <- compare_accuracy(eu_xgb,fluEuro.incidence,10,0,6,1)


#' Following is accuracy score of each individual country in Euro

Euro_accuracyIndi_15 <- NULL
for (i in 1:length(eu_xgb)){
  tmp <- compare_accuracy_indi(eu_xgb[i],fluEuro.incidence, 10,0,4,1)
  Euro_accuracyIndi_15 <- append(Euro_accuracyIndi_15, tmp)
  
}
