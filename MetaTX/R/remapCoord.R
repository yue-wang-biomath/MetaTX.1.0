remapCoord <-
  function(features,
           txdb               = NA,
           num_bin            = 10,
           includeNeighborDNA = FALSE,
           cds_by_tx0         = NA, 
           fiveUTR_tx0        = NA,
           threeUTR_tx0       = NA){
    
    if(includeNeighborDNA){

      num_bin_fiveUTR        <- num_bin
      num_bin_cds            <- num_bin
      num_bin_threeUTR       <- num_bin
      num_bin_promoter       <- num_bin
      num_bin_tail           <- num_bin  
  
      }else{

      num_bin_fiveUTR        <- num_bin
      num_bin_cds            <- num_bin
      num_bin_threeUTR       <- num_bin
            
    }
    
      # remap features to genome
      methyl                 <- features
      trans                  <- mapToTranscripts(methyl, txdb)
      trans_info             <- get_trans_function(trans, methyl)
    
    # get 5'UTR, CDS, 3'UTR alignment
    num_bin_sum            <- num_bin_cds + num_bin_fiveUTR + num_bin_threeUTR
    cds_align              <- get_align_function(num_bin_cds, trans_info, cds_by_tx0)
    fiveUTR_align          <- get_align_function(num_bin_fiveUTR, trans_info, fiveUTR_tx0)
    threeUTR_align         <- get_align_function(num_bin_threeUTR, trans_info, threeUTR_tx0)
    # get alignment matrix
    align_mtr_5cds3        <- data.frame(fiveUTR_align[, paste0('coordinate.', 1:num_bin_fiveUTR)]
                                         , cds_align[, paste0('coordinate.', 1:num_bin_cds)]
                                         , threeUTR_align[, paste0('coordinate.', 1:num_bin_threeUTR)])
    names(align_mtr_5cds3) <- paste0('coordinate.', 1:num_bin_sum)
    
    ## get width matrix
    width_mtr_5cds3        <- data.frame(coordinate = cbind(replicate(num_bin_fiveUTR, fiveUTR_align[, 'width'] / num_bin_fiveUTR)
                                                            , replicate(num_bin_cds, cds_align[, 'width'] / num_bin_cds)
                                                            , replicate(num_bin_threeUTR, threeUTR_align[, 'width'] / num_bin_threeUTR)))
    
    if(includeNeighborDNA){
      ## get promoter and tail
      width_promoter         <- 1000
      width_tail             <- 1000
      num_bin_sum            <- num_bin_sum + num_bin_promoter + num_bin_tail 
      promoter_align         <- get_promoter_function(width_promoter, num_bin_promoter, fiveUTR_align, fiveUTR_tx0, methyl)
      tail_align             <- get_tail_function(width_tail, num_bin_tail, threeUTR_align, threeUTR_tx0, methyl)
      
      # get alignment matrix (with promoter and tail)
      align_mtr              <- data.frame(coordinate = cbind(promoter_align[, paste0('coordinate.', 1:num_bin_promoter)]
                                                              , fiveUTR_align[, paste0('coordinate.', 1:num_bin_fiveUTR)]
                                                              , cds_align[, paste0('coordinate.', 1:num_bin_cds)]
                                                              , threeUTR_align[, paste0('coordinate.', 1:num_bin_threeUTR)]
                                                              , tail_align[, paste0('coordinate.', 1:num_bin_tail)]))
      names(align_mtr)       <- paste0('coordinate.', 1:num_bin_sum)
      # get width matrix (with promoter and tail)
      width_mtr              <- data.frame(coordinate = cbind(replicate(num_bin_promoter, promoter_align[, 'width'] / num_bin_promoter)
                                                              , replicate(num_bin_fiveUTR, fiveUTR_align[, 'width'] / num_bin_fiveUTR)
                                                              , replicate(num_bin_cds, cds_align[, 'width'] / num_bin_cds)
                                                              , replicate(num_bin_threeUTR, threeUTR_align[, 'width'] / num_bin_threeUTR)
                                                              , replicate(num_bin_tail, tail_align[, 'width'] / num_bin_tail)))
    }else{
      align_mtr              <- align_mtr_5cds3
      width_mtr              <- width_mtr_5cds3
    }
    
    return(list(align_mtr, width_mtr, trans_info))
  }
