load.iiag.data.fluView <- function(datadir="../iiag_data/data_old/") {
  
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
  
  ## Define the strings for all the files needed and read 
  # fid_this <- read.csv(paste(datadir,"/2019-2020_FluIDData.csv",sep=""))
  fview_ILINet <- read.csv(paste0(datadir,"/ILINet.csv.csv",sep=""))
  
  ## Fix names for fluView ILINet
  names(fview_ILINe) <- fix_headers(fview_ILINe)
  
  ## Use rbind to make the large tables. Should throw an error if the column
  ## names change in the future.
  # dfId <- rbind(fid_old_0,fid_old_2,fid_old_2,fid_old_3,fid_this)
  # dfNet <- rbind(fnet_old_0,fnet_old_2,fnet_old_2,fnet_old_3,fnet_this)
  dfId <- rbind(fid_old_1,fid_old_2,fid_old_3)
  # dfNet <- rbind(fnet_old_1,fnet_old_2,fnet_old_3,fnet_this)
  
  ## Sort both dataframe just incase
  dfId <- dfId[order(dfId$ISO2,dfId$ISO_YEAR,dfId$ISO_WEEK),]
  # dfNet <- dfNet[order(dfNet$ISO2,dfNet$ISO_YEAR,dfNet$ISO_WEEK),]
  
  ## Return the two datasets as a list
  # list(lab=dfNet,synd=dfId)
  fview_ILINet
}