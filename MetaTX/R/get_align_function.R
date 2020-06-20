get_align_function <-
  function(num_bin, trans_info, tx0){
    cds_name_all         <- names(tx0)
    cds_start_all        <- start(tx0)
    cds_end_all          <- end(tx0)
    cds_width_all        <- width(tx0)
    
    cds_index            <- trans_info[, 'index_trans']
    cds_exist_index      <- match(intersect(trans_info[, 'trans_ID'], as.double(cds_name_all)),
                                  trans_info[, 'trans_ID'])
    cds_non_exist_index  <- setdiff(cds_index, cds_exist_index)
    
    trans_cds_exist      <- trans_info[cds_exist_index, ]
    trans_cds_non_exist  <- trans_info[cds_non_exist_index, ]
    
    trans_cds_ID         <- trans_cds_exist[, 'trans_ID']
    trans_cds_ID         <- as.character(trans_cds_ID)
    cds_start            <- cds_start_all[trans_cds_ID]
    cds_end              <- cds_end_all[trans_cds_ID]
    cds_width            <- cds_width_all[trans_cds_ID]
    trans_cds_exist      <- data.frame(trans_cds_exist,
                                       width = sum(cds_width))
    methyl_pos           <- as.vector(trans_cds_exist[, 'methyl_pos'])
    
    pstv_index           <- which(trans_cds_exist[, 'strand'] == '+')
    ngtv_index           <- which(trans_cds_exist[, 'strand'] == '-')
    methyl_pos_pstv      <- methyl_pos[pstv_index]
    methyl_pos_ngtv      <- methyl_pos[ngtv_index]
    trans_cds_exist_pstv <- trans_cds_exist[pstv_index, ]
    trans_cds_exist_ngtv <- trans_cds_exist[ngtv_index, ]
    cds_start_pstv       <- cds_start[pstv_index, ]
    cds_end_pstv         <- cds_end[pstv_index, ]
    cds_start_ngtv       <- cds_start[ngtv_index, ]
    cds_end_ngtv         <- cds_end[ngtv_index, ]
    cds_width_pstv       <- cds_width[pstv_index, ]
    cds_width_ngtv       <- cds_width[ngtv_index, ]
    cds_width_pstv_sum   <- sum(cds_width_pstv)
    cds_width_ngtv_sum   <- sum(cds_width_ngtv)
    
    # pstv
    
    align_index_temp1    <- lapply(cds_start_pstv <= methyl_pos_pstv, '+')
    align_index_temp2    <- lapply(methyl_pos_pstv <= cds_end_pstv, '+')
    align_index_temp     <- mapply('+', align_index_temp1, align_index_temp2)
    align_index_1        <- lapply(align_index_temp, function(x){return(which(x==2))})
    align_index_1        <- IntegerList(align_index_1)
    dist_from_start      <- trans_cds_exist_pstv[, 'methyl_pos'] - cds_start_pstv[align_index_1] + 1
    
    align_index_2        <- lapply(cds_start_pstv <= methyl_pos_pstv, function(x){return(which(x==1))})
    align_index_before   <- IntegerList(align_index_2) - 1
    
    width_from_start     <-  IntegerList(mapply(function(x,y){return(y[x])}, align_index_before, cds_width_pstv))
    
    dist_from_start      <- as.numeric(as.vector(sum(width_from_start))+ dist_from_start)
    align_mtr_pstv       <- matrix(0, nrow(trans_cds_exist_pstv), num_bin)
    align_mtr_index      <- ceiling(dist_from_start / trans_cds_exist_pstv[, 'width'] * num_bin)
    
    
    if(nrow(align_mtr_pstv) != 0){
      for (i in 1:nrow(align_mtr_pstv)){
        align_mtr_pstv[i, align_mtr_index[i]] <- 1
      }
      align_mtr_pstv       <- data.frame(trans_cds_exist_pstv, coordinate = align_mtr_pstv)
    }
    
    # ngtv
    
    align_index_temp1    <- lapply(cds_start_ngtv <= methyl_pos_ngtv, '+')
    align_index_temp2    <- lapply(cds_end_ngtv >= methyl_pos_ngtv, '+')
    align_index_temp     <- mapply('+', align_index_temp1, align_index_temp2)
    align_index_1        <- lapply(align_index_temp, function(x){return(which(x==2))})
    align_index_1        <- IntegerList(align_index_1)
    dist_from_start      <- cds_end_ngtv[align_index_1] - trans_cds_exist_ngtv[, 'methyl_pos'] +1
    
    align_index_temp3    <- lapply(cds_start_ngtv >= methyl_pos_ngtv, '+')
    align_index_2        <- lapply(align_index_temp3, function(x){return(which(x==1))})
    align_index_before   <- IntegerList(align_index_2) - 1
    
    width_from_start     <- IntegerList(mapply(function(x,y){return(y[x])}, align_index_before, cds_width_ngtv))
    
    dist_from_start      <- as.numeric(as.vector(sum(width_from_start))+ dist_from_start)
    align_mtr_ngtv       <- matrix(0, nrow(trans_cds_exist_ngtv), num_bin)
    align_mtr_index      <- ceiling(dist_from_start / trans_cds_exist_ngtv[, 'width'] * num_bin)
    
    if(nrow(align_mtr_ngtv) != 0){
      for (i in 1:nrow(align_mtr_ngtv)){
        align_mtr_ngtv[i, align_mtr_index[i]] <- 1
      }
      align_mtr_ngtv       <- data.frame(trans_cds_exist_ngtv, coordinate = align_mtr_ngtv)
      
    }
    
    # trans_cds_non_exist
    align_mtr_non_cds    <- matrix(0, nrow(trans_cds_non_exist), num_bin)
    align_mtr_non_cds    <- data.frame(trans_cds_non_exist,
                                       width = matrix(0, nrow(trans_cds_non_exist), 1),
                                       coordinate = align_mtr_non_cds)
    # Sort
    align_cds_mtr        <- rbind(align_mtr_pstv, align_mtr_ngtv, align_mtr_non_cds)
    align_cds_mtr        <- align_cds_mtr[order(align_cds_mtr[, 'index_trans']),]
    
    return(align_cds_mtr)
  }

