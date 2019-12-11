load.iiag.data.fluView <- function(datadir="../fluView_data") {
  
  ## Helper function to fix some header names to be used below
  fix_headers <- function(x){
    x <- as.character(x)
    nchar <- unlist(strsplit(x, " "))
    if (length(nchar) > 1){
      if (length(unlist(strsplit(nchar, "."))) > 1){
        x <- unlist(strsplit(nchar, ".", fixed = TRUE))
      }
      
      if (length(unlist(strsplit(nchar, "%", fixed = TRUE))) > 1){
        x <- str_replace(x, "%", "PERCENTAGE OF ")
        x <- unlist(strsplit(x, " "))
      }
      if (length(unlist(strsplit(nchar, "-"))) > 1){
        x <- unlist(strsplit(nchar,"-"))
      }
    }else{
      x <- nchar
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

#' Extract incidence 