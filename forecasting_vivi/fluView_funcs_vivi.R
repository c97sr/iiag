load.iiag.data.fluView <- function(datadir="../fluView_data") {
  
  ## Helper function to fix some header names to be used below
  fix_headers <- function(x){
    x <- as.character(x)
    if (length(unlist(strsplit(x," "))) > 1){
      if (length(unlist(strsplit(x, "-"))) > 1){
        x <- unlist(strsplit(x,"-"))
      }
      if (length(unlist(strsplit(x, ".", fixed = TRUE))) > 1){
        x <- unlist(strsplit(x, ".", fixed = TRUE))
      }
      
      if (length(unlist(strsplit(x, "%", fixed = TRUE))) > 1){
        x <- str_replace(x, "%", "PERCENTAGE OF")
        x <- unlist(strsplit(x, " "))
      }
      if (length(unlist(strsplit(x, " "))) > 1){
        x <- unlist(strsplit(x, " "))
      }
    }else{
      x <- x
    }
    x <- paste0(x, collapse = "_")
    x
  }
  ## Define the strings for all the files needed and read 
 
  fview_ILINet <- read.csv(paste0(datadir,"/ILINet.csv"))
  
  ## Fix names for fluView ILINet
  curnames <- unlist(fview_ILINet[1,], use.names=FALSE)
  newnames <- c()
  for (i in 1:length(curnames)){
    tmp <- fix_headers(curnames[i])
    newnames <- append(newnames,tmp)
  }
    
  names(fview_ILINet) <- newnames
  fview_ILINet <- fview_ILINet[-1,]
  
  fview_ILINet
}

fix_headers <- function(x){
  x <- as.character(x)
  if (length(unlist(strsplit(x," "))) > 1){
    if (length(unlist(strsplit(x, "-"))) > 1){
      x <- unlist(strsplit(x,"-"))
    }
    if (length(unlist(strsplit(x, ".", fixed = TRUE))) > 1){
      x <- unlist(strsplit(x, ".", fixed = TRUE))
    }
    
    if (length(unlist(strsplit(x, "%", fixed = TRUE))) > 1){
      x <- str_replace(x, "%", "PERCENTAGE OF")
      x <- unlist(strsplit(x, " "))
    }
    if (length(unlist(strsplit(x, " "))) > 1){
      x <- unlist(strsplit(x, " "))
    }
  }else{
    x <- x
  }
  x <- paste0(x, collapse = "_")
  x
}

header <- as.character(unlist(fluview[1,],use.names=FALSE))
headernew <- c()
for (i in 1:length(header)){
  tmp <- fix_headers(header[i])
  headernew <- append(headernew, tmp)
}

#' Extract incidence 
extract.incidence.fluView <- function(fluView_data,
                                      sel_states,
                                      # sel_ag,
                                      # sel_measure,
                                      minYear,
                                      maxYear) {
  ## reorder data by country alphabetically
  fluView_data <- fluView_data[order(fluView_data$REGION),]
  
  ## Setup the week scale in a format consistent with the week format
  ## in the data and cope with 53-week years. Needs the list of 53 week years
  ## extending in both directions.
  ## Perhaps should have a few lines to get rid of data NAs and avoid a warning
  ## at the next line?
  fluView_data$YRWEEK  <- paste(fluView_data$YEAR,sprintf("%02d",as.numeric(as.character(fluView_data$WEEK))),sep="-")

  min(as.numeric(as.character(fluView_data$YEAR)))
  yrs53Weeks <- c(2015,2020)
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
  rtnmat <- matrix(data=NA,nrow=length(sel_weeks),ncol=length(sel_states))
  colnames(rtnmat) <- sel_states
  rownames(rtnmat) <- sel_weeks
 
  ## Start outer loop over the country codes
  for (cur_iso3 in sel_states) {
    
    ## Define criteria and subset the data
    crit1 <- (fluView_data$REGION == cur_iso3)
    # if(!("AGEGROUP_CODE" %in% colnames(dfId))) {
      # crit2 <- TRUE
    # } else {
      # crit2 <- (dfId$AGEGROUP_CODE %in% sel_ag)
    # }
    
    # crit3 <- (dfId$MEASURE_CODE %in% sel_measure)
    tmpdf <- fluView_data[crit1,]
    tmpdf <- tmpdf[order(tmpdf$YRWEEK),]
    
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
        sel_weeks[cur_ind_rtn] != tmpdf$YRWEEK[cur_ind_df] &&
        cur_ind_rtn <= max_ind_rtn
      ) {
        cur_ind_rtn <- cur_ind_rtn + 1
      }
      
      if (cur_ind_rtn <= max_ind_rtn) {
        val_rtn <- rtnmat[cur_ind_rtn,cur_iso3]
        val_df <- tmpdf$ILITOTAL[cur_ind_df]
        if (as.character(val_df) == "X"){
          val_df <- as.character(val_df)
        }else{
          val_df <- as.numeric(as.character(val_df))
        }
        if (!is.na(val_df)) {
          if (is.na(val_rtn)) {
            rtnmat[cur_ind_rtn,cur_iso3] <- val_df
          } 
        }
      }
      cur_ind_df <- cur_ind_df + 1
    }
  
    ## Close the state-level loop
  }
  
  ## Return the populated incidence matrix as only result of function
  rtnmat
}

while (cur_ind_df <= max_ind_df) {
  while (
    sel_weeks[cur_ind_rtn] != tmpdf$YRWEEK[cur_ind_df] &&
    cur_ind_rtn <= max_ind_rtn
  ) {
    cur_ind_rtn <- cur_ind_rtn + 1
  }
  if (cur_ind_rtn <= max_ind_rtn) {
    val_rtn <- rtnmat[cur_ind_rtn,cur_iso3]
    val_df <- tmpdf$ILITOTAL[cur_ind_df]
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
