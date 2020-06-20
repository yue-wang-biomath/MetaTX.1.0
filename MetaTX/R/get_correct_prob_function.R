get_correct_prob_function <-
  function(num_bin_sum, align_mtr, weight_mtr, trans_info){
    
    alpha                     <- matrix(1/num_bin_sum, 1, num_bin_sum)
    index_methyl              <- trans_info[, 'index_methyl']
    
    for (j in 1:20){
      numerator_prob          <- align_mtr * weight_mtr
      
      for (k in 1:num_bin_sum){
        numerator_prob[, k]   <- numerator_prob[, k] * alpha[k]
      }
      row_sum_prob            <- as.vector(rowSums(numerator_prob))
      row_sum_prob            <- data.frame(group_name = as.character(index_methyl), value = row_sum_prob)
      sum_prob                <- aggregate(row_sum_prob[,'value'], 
                                           by = list(group_name = factor(trans_info[,'index_methyl'], levels = unique(trans_info[,'index_methyl']))), 
                                           FUN = sum)
      names(sum_prob)         <- c('index_methyl', 'value')
      denominator_prob        <- sum_prob[, 'value']
      names(denominator_prob) <- sum_prob[, 'index_methyl']
      denominator_prob        <- denominator_prob[as.character(trans_info[, 'index_methyl'])]
      denominator_prob        <- replicate(num_bin_sum, denominator_prob)
      denominator_prob[denominator_prob == 0] <- 1
      
      prob_mtr                <- numerator_prob / denominator_prob 
      alpha_numerator         <- colSums(prob_mtr * weight_mtr)
      alpha                   <- alpha_numerator/sum(alpha_numerator)
    }
    
    
    return(list(alpha, prob_mtr))
  }
