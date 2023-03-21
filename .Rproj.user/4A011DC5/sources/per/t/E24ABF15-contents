################################################################################
# UTILITY FUNCTIONS
################################################################################
# These functions remain interior to the package (not exported)

#' @importFrom GenomicRanges start end width ranges GRanges mcols
#' @importFrom IRanges IntegerList flank findOverlaps IRanges
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
    methyl_pos           <- as.vector(trans_cds_exist[, 'feature_pos'])

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
    dist_from_start      <- trans_cds_exist_pstv[, 'feature_pos'] - cds_start_pstv[align_index_1] + 1

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
    dist_from_start      <- cds_end_ngtv[align_index_1] - trans_cds_exist_ngtv[, 'feature_pos'] +1

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

    promoter_align_pstv       <- data.frame(fiveUTR_align_pstv[, c('index_trans', 'index_feature', 'seqnames', 'feature_pos', 'strand', 'trans_ID')]
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

    promoter_align_ngtv       <- data.frame(fiveUTR_align_ngtv[, c('index_trans', 'index_feature', 'seqnames', 'feature_pos', 'strand', 'trans_ID')]
                                            , width = replicate(nrow(align_mtr_ngtv), width_promoter)
                                            , coordinate = align_mtr_ngtv)

    #promoter_non_exsit
    promoter_align_non_exist  <- data.frame(fiveUTR_align_non_exist[, c('index_trans', 'index_feature', 'seqnames', 'feature_pos', 'strand', 'trans_ID', 'width')]
                                            , coordinate = matrix(0, nrow(fiveUTR_align_non_exist), num_bin_promoter))

    # Sort
    promoter_align            <- rbind(promoter_align_pstv, promoter_align_ngtv, promoter_align_non_exist)
    promoter_align            <- promoter_align[order(promoter_align[, 'index_trans']),]

    return(promoter_align)
  }

get_tail_function <-
  function(width_tail, num_bin_tail, threeUTR_align, threeUTR_tx0, methyl){

    threeUTR_align_non_exist   <- threeUTR_align[which(threeUTR_align[, 'width'] == 0), ]
    threeUTR_align_exist       <- threeUTR_align[which(threeUTR_align[, 'width'] != 0), ]
    threeUTR_align_pstv        <- threeUTR_align_exist[threeUTR_align_exist[, 'strand'] == '+', ]
    threeUTR_align_ngtv        <- threeUTR_align_exist[threeUTR_align_exist[, 'strand'] == '-', ]
    threeUTR_pstv              <- threeUTR_tx0[as.character(threeUTR_align_pstv[, 'trans_ID'])]
    threeUTR_ngtv              <- threeUTR_tx0[as.character(threeUTR_align_ngtv[, 'trans_ID'])]


    tail_pstv_start <- GRanges()
    tail_ngtv_start <- GRanges()
    if(nrow(threeUTR_align_pstv) != 0){
      tail_pstv_start            <- GRanges(seqnames = threeUTR_align_pstv[, 'seqnames'],
                                            ranges = IRanges(max(end(ranges(threeUTR_pstv))), width = 1),
                                            strand = '+')
    }
    if(nrow(threeUTR_align_ngtv) != 0){
      tail_ngtv_start            <- GRanges(seqnames = threeUTR_align_ngtv[, 'seqnames'],
                                            ranges = IRanges(min(start(ranges(threeUTR_ngtv))), width = 1),
                                            strand = '-')
    }

    tail_pstv_rgs              <- flank(tail_pstv_start, width_tail, start = FALSE)
    tail_ngtv_rgs              <- flank(tail_ngtv_start, width_tail, start = FALSE)

    methyl_tail_map_pstv       <- data.frame(findOverlaps(methyl, tail_pstv_rgs))
    methyl_tail_map_ngtv       <- data.frame(findOverlaps(methyl, tail_ngtv_rgs))

    #pstv
    dist_from_start            <- start(methyl[methyl_tail_map_pstv[, 'queryHits']]) - start(
      tail_pstv_rgs)[methyl_tail_map_pstv[, 'subjectHits']]
    align_mtr_index            <- ceiling(dist_from_start / width_tail * num_bin_tail)
    align_mtr_pstv             <- matrix(0, nrow(threeUTR_align_pstv), num_bin_tail)

    for (i in 1:nrow(align_mtr_pstv)){
      align_mtr_pstv[methyl_tail_map_pstv[, 'subjectHits'][i], align_mtr_index[i]] <- 1
    }

    tail_align_pstv            <- data.frame(threeUTR_align_pstv[, c('index_trans', 'index_feature', 'seqnames', 'feature_pos', 'strand', 'trans_ID')]
                                             , width = replicate(nrow(align_mtr_pstv), width_tail)
                                             , coordinate = align_mtr_pstv)

    #ngtv
    dist_from_start           <- end(tail_ngtv_rgs)[methyl_tail_map_ngtv[, 'subjectHits']] - start(
      methyl[methyl_tail_map_ngtv[, 'queryHits']])
    align_mtr_index           <- ceiling(dist_from_start / width_tail * num_bin_tail)
    align_mtr_ngtv            <- matrix(0, nrow(threeUTR_align_ngtv), num_bin_tail)

    for (i in 1:nrow(align_mtr_ngtv)){
      align_mtr_ngtv[methyl_tail_map_ngtv[, 'subjectHits'][i], align_mtr_index[i]] <- 1
    }

    tail_align_ngtv           <- data.frame(threeUTR_align_ngtv[, c('index_trans', 'index_feature', 'seqnames', 'feature_pos', 'strand', 'trans_ID')]
                                            , width = replicate(nrow(align_mtr_ngtv), width_tail)
                                            , coordinate = align_mtr_ngtv)

    #promoter_non_exsit
    tail_align_non_exist      <- data.frame(threeUTR_align_non_exist[, c('index_trans', 'index_feature', 'seqnames', 'feature_pos', 'strand', 'trans_ID', 'width')]
                                            , coordinate = matrix(0, nrow(threeUTR_align_non_exist), num_bin_tail))

    # Sort
    tail_align        <- rbind(tail_align_pstv, tail_align_ngtv, tail_align_non_exist)
    tail_align        <- tail_align[order(tail_align[, 'index_trans']),]

    return(tail_align)
  }

get_trans_function <-
  function(trans, methyl){
    sample_methyl_pos      <- start(ranges(methyl))
    sample_methyl_seqnames <- data.frame(seqnames(methyl))[[1]]
    trans_metadata         <- mcols(trans)

    index_trans            <- 1:length(trans)
    index_methyl           <- data.frame(trans)[, 'xHits']
    seqnames               <- sample_methyl_seqnames[data.frame(trans)[, 'xHits']]
    methyl_pos             <- sample_methyl_pos[data.frame(trans)[, 'xHits']]
    strand                 <- data.frame(trans)[, 'strand']
    trans_ID               <- data.frame(trans)[, 'transcriptsHits']

    trans_info             <- data.frame(index_trans   = index_trans,
                                         index_feature = index_methyl,
                                         seqnames      = seqnames,
                                         feature_pos   = methyl_pos,
                                         strand        = strand,
                                         trans_ID      = trans_ID)
    return(trans_info)
  }

get_correct_prob_function <-
  function(num_bin_sum, align_mtr, weight_mtr, trans_info){

    alpha                     <- matrix(1/num_bin_sum, 1, num_bin_sum)
    index_methyl              <- trans_info[, 'index_feature']

    for (j in 1:20){
      numerator_prob          <- align_mtr * weight_mtr

      for (k in 1:num_bin_sum){
        numerator_prob[, k]   <- numerator_prob[, k] * alpha[k]
      }
      row_sum_prob            <- as.vector(rowSums(numerator_prob))
      row_sum_prob            <- data.frame(group_name = as.character(index_methyl), value = row_sum_prob)
      sum_prob                <- aggregate(row_sum_prob[,'value'],
                                           by = list(group_name = factor(trans_info[,'index_feature'], levels = unique(trans_info[,'index_feature']))),
                                           FUN = sum)
      names(sum_prob)         <- c('index_feature', 'value')
      denominator_prob        <- sum_prob[, 'value']
      names(denominator_prob) <- sum_prob[, 'index_feature']
      denominator_prob        <- denominator_prob[as.character(trans_info[, 'index_feature'])]
      denominator_prob        <- replicate(num_bin_sum, denominator_prob)
      denominator_prob[denominator_prob == 0] <- 1

      prob_mtr                <- numerator_prob / denominator_prob
      alpha_numerator         <- colSums(prob_mtr * weight_mtr)
      alpha                   <- alpha_numerator/sum(alpha_numerator)
    }


    return(list(alpha, prob_mtr))
  }

get_prob_function <-
  function(align_mtr, width_mtr){
    alpha <- colSums(align_mtr)
    alpha <- alpha/sum(alpha)
    return(alpha)
  }


