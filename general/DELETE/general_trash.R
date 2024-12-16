
GetKdsErrorOrig <- function(mirna, experiment, start, stop, sitelist, min_method, ci,
                   scaled=TRUE) {
      sitesXcounts <- GetSitesXCountsCombined(mirna, experiment, start, stop,
                                              sitelist)
      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                         experiment, "/kds_PAPER/", start, "-", stop, "_", 
                         sitelist, "_error_PAPER.txt")
    # Remove the first column, which has the sequences of the site types.
    sitesXcounts <- sitesXcounts[,-1]
    # Remove any rows with no reads whoatsoever.
    sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
    # Remove any rows for which there are no input reads:
    sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]

    params <- read.table(params.file, header=TRUE, row.names=1,stringsAsFactors = FALSE)
    # print(params)
    kds <- rbind(Logistic(params[1:(nrow(sitesXcounts)-1),], 10),rep(1,length=4))
    rownames(kds)[nrow(kds)] <- "None"
    bgs <- 10^params["bg",]
    k.c.stockago <- 10^params["AGO",]
    pars_all <- rbind(kds, bgs, k.c.stockago)
    # print(pars_all)
    if (scaled == TRUE) {
      return(pars_all)
    } else {
      return(params)   
    }
}



GetKdsError <- function(mirna, experiment, start, stop, sitelist, min_method, ci,
                   scaled=TRUE,combined.input=TRUE,nosite=FALSE) {
      sitesXcounts <- GetSitesXCountsCombined(mirna, experiment, start, stop,
                                                 sitelist)
          if (nosite == TRUE){
        min_method <- paste0("nosite_", min_method)
      }

      if (combined.input==TRUE){
              params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                         experiment, "/kds_PAPER/", start, "-", stop, "_", 
                         sitelist, "_singlebg_combinedinput_", min_method, "_ci", ci,
                         "_PAPER.txt")

              } else {
                              params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                         experiment, "/kds_PAPER/", start, "-", stop, "_", 
                         sitelist, "_singlebg_", min_method, "_ci", ci,
                         "_PAPER.txt")

              }
    # Remove the first column, which has the sequences of the site types.
    sitesXcounts <- sitesXcounts[,-1]
    # Remove any rows with no reads whoatsoever.
    sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
    # Remove any rows for which there are no input reads:
    sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]

    params <- read.table(params.file, header=TRUE, row.names=1,stringsAsFactors = FALSE)
    print(params)
      if (nosite == FALSE){
      kds <- params[1:(nrow(sitesXcounts)-1),]   
    } else {
     kds <- params[1:nrow(sitesXcounts),]
     }
     print(kds)
    bgs <- params["bg",]
    k.c.stockago <- params["AGO",]
    pars_all <- rbind(kds, bgs, k.c.stockago)
    print(pars_all)
    if (scaled == TRUE) {
      return(pars_all)
    } else {
      return(params)   
    }
}

GetKdsErrorDist <- function(mirna, experiment, start, stop, sitelist, min_method,
                   scaled=TRUE) {
      sitesXcounts <- GetSitesXCountsCombined(mirna, experiment, start, stop,
                                              sitelist)
      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                         experiment, "/kds_PAPER/", start, "-", stop, "_", 
                         sitelist,
                         "_singlebg_logresiduals_bootstraps_PAPER.txt")
    # Remove the first column, which has the sequences of the site types.
    sitesXcounts <- sitesXcounts[,-1]
    # Remove any rows with no reads whoatsoever.
    sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
    # Remove any rows for which there are no input reads:
    sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]

    params <- read.table(params.file, header=TRUE, row.names=1,stringsAsFactors = FALSE)
    # print(params)
    kds <- rbind(Logistic(params[1:(nrow(sitesXcounts)-1),], 10),rep(1,length=4))
    rownames(kds)[nrow(kds)] <- "None"
    bgs <- 10^params["bg",]
    k.c.stockago <- 10^params["AGO",]
    pars_all <- rbind(kds, bgs, k.c.stockago)
    # print(pars_all)
    if (scaled == TRUE) {
      return(pars_all)
    } else {
      return(params)   
    }
}
