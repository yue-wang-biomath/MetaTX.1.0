#' Calculate multiple groups of features.
#' @export remapCoordMulti
remapCoordMulti <-
function(features_list,
           txdb               = NA,
           num_bin            = 10,
           includeNeighborDNA = TRUE){
    txdb1 <- txdb
    includeNeighborDNA1 <- includeNeighborDNA
    num_bin1 <- num_bin
    Iter <- length(features_list)
    
    remapresults_list <-list() 
    for(i in 1: Iter){
      remapresults_list[[i]] <- 
        remapCoord(features = features_list[[i]],
                   txdb = txdb1,
                   num_bin = num_bin1,
                   includeNeighborDNA = includeNeighborDNA1)
    }
    return(remapresults_list)
}
