load.iiag.data.fluView <- function(datadir="../fluView_data") {
  
  ## Helper function to fix some header names to be used below
  fix_headers <- function(x){
    nchar <- sapply(strsplit(x, " "), length)
    if (nchar > 1){
      if (length(unlist(strsplit(x,"-"))) > 1){
        x <- unlist(strsplit(unlist(strsplit(x," ")),"-"))
      }
    }
  }
  
  ## Define the strings for all the files needed and read 
 
  fview_ILINet <- read.csv(paste0(datadir,"/ILINet.csv"))
  
  ## Fix names for fluView ILINet
  newnames <- unlist(fview_ILINet[1,], use.names=FALSE)
  names(fview_ILINet) <- newnames

  fview_ILINet <- fview_ILINet[-1,]
  
  fview_ILINet
}
