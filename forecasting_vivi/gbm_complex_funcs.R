#' XGBoost uses a matrix of input data instead of a data frame,
#' so the output of gbm_complex should be a matrix of dataset.
gbm_complex <- function(data, country, num_category,nWeek_ahead){
  yr <- seq(2010, 2017, by = 1)
  initial_data <- c()
  for (i in 1:length(yr)){
    tmp <- extract_incidence.idd(data, country_code = country, yr[i])
    initial_data <- rbind(initial_data, tmp)
  }
  initial_data <- as.data.frame(initial_data)
  
  # Divide incidence into 10 categories
  initial_data2 <- cbind(initial_data, cut_interval(initial_data$incidence, n=num_category))
  colnames(initial_data2)[4] <- 'category'
  
  # match the intervals with numeric categories
  levels_index <- levels(initial_data2$category)
  for(i in 1:nrow(initial_data2)){
    x <- initial_data2$category[i]
    if (is.na(x)==FALSE){
      y <- which(levels(initial_data2$category)==x)
      initial_data2$category2[i] <- as.numeric(y)
    }else{
      initial_data2$category2[i] <- NA
    }
  }
  
  # convert into suitable stucture for gbm
  incidence_gbm <- matrix(NA, nrow = (dim(initial_data2)[1]-5), ncol=3)
  
  if(nWeek_ahead == 1){
    for (i in 1:nrow(incidence_gbm)){
      incidence_gbm[i,] <- initial_data2$category2[c((i+5),(i+4),(i+3))]
    }
    
    incidence_gbm <- as.data.frame(incidence_gbm)
    colnames(incidence_gbm) <- c('Y_week0','week_1','week_2')
    rownames(incidence_gbm) <- as.character(initial_data2$time_name[6:nrow(initial_data2)])
  }
  
  if(nWeek_ahead == 2){
    for (i in 1:nrow(incidence_gbm)){
      incidence_gbm[i,] <- initial_data2$category2[c((i+5),(i+3),(i+2))]
    }
    
    incidence_gbm <- as.data.frame(incidence_gbm)
    colnames(incidence_gbm) <- c('Y_week0','week_2','week_3')
    rownames(incidence_gbm) <- as.character(initial_data2$time_name[6:nrow(initial_data2)])
  }
  
  if(nWeek_ahead == 3){
    for (i in 1:nrow(incidence_gbm)){
      incidence_gbm[i,] <- initial_data2$category2[c((i+5),(i+2),(i+1))]
    }
    
    incidence_gbm <- as.data.frame(incidence_gbm)
    colnames(incidence_gbm) <- c('Y_week0','week_3','week_4')
    rownames(incidence_gbm) <- as.character(initial_data2$time_name[6:nrow(initial_data2)])
  }
  
  
  if(nWeek_ahead == 4){
    for (i in 1:nrow(incidence_gbm)){
      incidence_gbm[i,] <- initial_data2$category2[c((i+5),(i+1),i)]
    }
    
    incidence_gbm <- as.data.frame(incidence_gbm)
    colnames(incidence_gbm) <- c('Y_week0','week_4','week_5')
    rownames(incidence_gbm) <- as.character(initial_data2$time_name[6:nrow(initial_data2)])
  }
  
  
  # add potential covariates
  # month: Jan - Dec
  # incidence_gbm$week_time <- as.character(initial_data2$time_name[3:nrow(initial_data2)])
  week_time <- rownames(incidence_gbm)
  incidence_gbm$month <- NA
  
  jan_row <- c()
  jan_week <- paste("-0",1:5, sep = "")
  for (i in 1:length(jan_week)){
    tmp <- grep(jan_week[i], rownames(incidence_gbm))
    jan_row <- append(jan_row, tmp)
  }
  incidence_gbm$month[jan_row] <- "January"
  
  feb_row <- c()
  feb_week <- paste("-0",6:9, sep = "")
  for (i in 1:length(feb_week)){
    tmp <- grep(feb_week[i], rownames(incidence_gbm))
    feb_row <- append(feb_row, tmp)
  }
  incidence_gbm$month[feb_row] <- "Feburary"
  
  mar_row <- c()
  mar_week <- paste("-", 10:14, sep = "")
  for (i in 1:length(mar_week)){
    tmp <- grep(mar_week[i], rownames(incidence_gbm))
    mar_row <- append(mar_row, tmp)
  }
  incidence_gbm$month[mar_row] <- "March"
  
  apr_row <- c()
  apr_week <- paste("-", 15:18, sep = "")
  for (i in 1:length(apr_week)){
    tmp <- grep(apr_week[i], rownames(incidence_gbm))
    apr_row <- append(apr_row, tmp)
  }
  incidence_gbm$month[apr_row] <- "April"
  
  may_row <- c()
  may_week <- paste("-", 19:23, sep = "")
  for (i in 1:length(may_week)){
    tmp <- grep(may_week[i], rownames(incidence_gbm))
    may_row <- append(may_row, tmp)
  }
  incidence_gbm$month[may_row] <- "May"
  
  june_row <- c()
  june_week <- paste("-", 24:27, sep = "")
  for (i in 1:length(june_week)){
    tmp <- grep(june_week[i], rownames(incidence_gbm))
    june_row <- append(june_row, tmp)
  }
  incidence_gbm$month[june_row] <- "June"
  
  july_row <- c()
  july_week <- paste("-", 28:32, sep = "")
  for (i in 1:length(july_week)){
    tmp <- grep(july_week[i], rownames(incidence_gbm))
    july_row <- append(july_row, tmp)
  }
  incidence_gbm$month[july_row] <- "July"
  
  aug_row <- c()
  aug_week <- paste("-", 33:36, sep = "")
  for (i in 1:length(aug_week)){
    tmp <- grep(aug_week[i], rownames(incidence_gbm))
    aug_row <- append(aug_row, tmp)
  }
  incidence_gbm$month[aug_row] <- "August"
  
  sep_row <- c()
  sep_week <- paste("-", 37:40, sep = "")
  for (i in 1:length(sep_week)){
    tmp <- grep(sep_week[i], rownames(incidence_gbm))
    sep_row <- append(sep_row, tmp)
  }
  incidence_gbm$month[sep_row] <- "September"
  
  oct_row <- c()
  oct_week <- paste("-", 41:44, sep = "")
  for (i in 1:length(oct_week)){
    tmp <- grep(oct_week[i], rownames(incidence_gbm))
    oct_row <- append(oct_row, tmp)
  }
  incidence_gbm$month[oct_row] <- "October"
  
  nov_row <- c()
  nov_week <- paste("-", 45:48, sep = "")
  for (i in 1:length(nov_week)){
    tmp <- grep(nov_week[i], rownames(incidence_gbm))
    nov_row <- append(nov_row, tmp)
  }
  incidence_gbm$month[nov_row] <- "November"
  
  dec_row <- c()
  dec_week <- paste("-", 49:52, sep = "")
  for (i in 1:length(dec_week)){
    tmp <- grep(dec_week[i], rownames(incidence_gbm))
    dec_row <- append(dec_row, tmp)
  }
  incidence_gbm$month[dec_row] <- "December"
  
  # season 
  incidence_gbm$season <- NA
  for(i in 1:nrow(incidence_gbm)){
    if(incidence_gbm$month[i] %in% c("September","October","November")==TRUE){
      incidence_gbm$season[i] <-"autumn"
    }
    if(incidence_gbm$month[i] %in% c("June","July","August")==TRUE){
      incidence_gbm$season[i] <- "summer"
    }
    if(incidence_gbm$month[i] %in% c("December","January","February")==TRUE){
      incidence_gbm$season[i] <- "winter"
    } 
    if (incidence_gbm$month[i] %in% c("March","April","May")==TRUE){
      incidence_gbm$season[i] <- "spring"
    }
  }
  
  # deal with NA data
  b <- c()
  for(i in 1:nrow(incidence_gbm)){
    if(sum(is.na(incidence_gbm[i,])==TRUE) >= 1){
      tmp <- i
      b <- append(b,tmp)
    }
    b
  }
  
  incidence_gbm <- incidence_gbm[-b, ]
  incidence_gbm <- as.data.frame(incidence_gbm)
  incidence_gbm
}