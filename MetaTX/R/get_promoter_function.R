get_promoter_function <-
  function(width_promoter, num_bin_promoter, fiveUTR_align, fiveUTR_tx0, methyl){
    
    fiveUTR_align_non_exist   <- fiveUTR_align[ which(fiveUTR_align[, 'width'] == 0), ]
    fiveUTR_align_exist       <- fiveUTR_align[ which(fiveUTR_align[, 'width'] != 0), ]
    fiveUTR_align_pstv        <- fiveUTR_align_exist[fiveUTR_align_exist[, 'strand'] == '+', ]
    fiveUTR_align_ngtv        <- fiveUTR_align_exist[fiveUTR_align_exist[, 'strand'] == '-', ]
    fiveUTR_pstv              <- fiveUTR_tx0[as.character(fiveUTR_align_pstv[, 'trans_ID'])]
    fiveUTR_ngtv              <- fiveUTR_tx0[as.character(fiveUTR_align_ngtv[, 'trans_ID'])]
    
    promoter_pstv_start <- GRanges()
    promoter_ngtv_start <- GRanges()
    if(nrow(fiveUTR_align_pstv) != 0){
      promoter_pstv_start       <- GRanges(seqnames = fiveUTR_align_pstv[, 'seqnames'], 
                                           ranges = IRanges(min(start(ranges(fiveUTR_pstv))), width = 1), 
                                           strand = '+')
    }
    
    if(nrow(fiveUTR_align_ngtv) != 0){
      promoter_ngtv_start       <- GRanges(seqnames = fiveUTR_align_ngtv[, 'seqnames'], 
                                           ranges = IRanges(max(end(ranges(fiveUTR_ngtv))), width = 1), 
                                           strand = '-')
    }
    
    promoter_pstv_rgs         <- flank(promoter_pstv_start, width_promoter, start = TRUE)
    promoter_ngtv_rgs         <- flank(promoter_ngtv_start, width_promoter, start = TRUE)
    
    methyl_promoter_map_pstv  <- data.frame(findOverlaps(methyl, promoter_pstv_rgs))
    methyl_promoter_map_ngtv  <- data.frame(findOverlaps(methyl, promoter_ngtv_rgs))
    
    #pstv
    dist_from_start           <- start(methyl[methyl_promoter_map_pstv[, 'queryHits']]) - start(
      promoter_pstv_rgs)[methyl_promoter_map_pstv[, 'subjectHits']]
    align_mtr_index           <- ceiling(dist_from_start / width_promoter * num_bin_promoter)
    align_mtr_pstv            <- matrix(0, nrow(fiveUTR_align_pstv), num_bin_promoter)
    
    for (i in 1:nrow(align_mtr_pstv)){
      align_mtr_pstv[methyl_promoter_map_pstv[, 'subjectHits'][i], align_mtr_index[i]] <- 1
    }
    
    promoter_align_pstv       <- data.frame(fiveUTR_align_pstv[, c('index_trans', 'index_methyl', 'seqnames', 'methyl_pos', 'strand', 'trans_ID')]
                                            , width = replicate(nrow(align_mtr_pstv), width_promoter)
                                            , coordinate = align_mtr_pstv)
    
    #ngtv
    dist_from_start           <- end(promoter_ngtv_rgs)[methyl_promoter_map_ngtv[, 'subjectHits']] - start(
      methyl[methyl_promoter_map_ngtv[, 'queryHits']]) 
    align_mtr_index           <- ceiling(dist_from_start / width_promoter * num_bin_promoter)
    align_mtr_ngtv            <- matrix(0, nrow(fiveUTR_align_ngtv), num_bin_promoter)
    
    for (i in 1:nrow(align_mtr_ngtv)){
      align_mtr_ngtv[methyl_promoter_map_ngtv[, 'subjectHits'][i], align_mtr_index[i]] <- 1
    }
    
    promoter_align_ngtv       <- data.frame(fiveUTR_align_ngtv[, c('index_trans', 'index_methyl', 'seqnames', 'methyl_pos', 'strand', 'trans_ID')]
                                            , width = replicate(nrow(align_mtr_ngtv), width_promoter)
                                            , coordinate = align_mtr_ngtv)
    
    #promoter_non_exsit
    promoter_align_non_exist  <- data.frame(fiveUTR_align_non_exist[, c('index_trans', 'index_methyl', 'seqnames', 'methyl_pos', 'strand', 'trans_ID', 'width')]
                                            , coordinate = matrix(0, nrow(fiveUTR_align_non_exist), num_bin_promoter))
    
    # Sort
    promoter_align            <- rbind(promoter_align_pstv, promoter_align_ngtv, promoter_align_non_exist)
    promoter_align            <- promoter_align[order(promoter_align[, 'index_trans']),]
    
    return(promoter_align)
  }
