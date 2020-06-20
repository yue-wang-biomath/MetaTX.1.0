isoformProb <-
function(remap_results,
         num_bin = 10,
        includeNeighborDNA = FALSE,
        lambda = 2     
    ){
        align_mtr          <- remap_results[[1]]
        width_mtr          <- remap_results[[2]]
        trans_info         <- remap_results[[3]]
        num_bin_sum        <- ncol(align_mtr)
        
        if(includeNeighborDNA){
          weight_start     <- num_bin + 1
          weight_end       <- num_bin * 4
          num_bin_sum      <- num_bin * 5
          weight_mtr       <- replicate(num_bin_sum, rowSums(width_mtr[, weight_start:weight_end]))    
          prob_mtr         <- get_correct_prob_function(num_bin_sum, align_mtr, weight_mtr^lambda, trans_info)[[2]]
         }else{
          num_bin_sum      <- num_bin * 3
          weight_mtr       <- replicate(num_bin_sum, rowSums(width_mtr))
          prob_mtr         <- get_correct_prob_function(num_bin_sum, align_mtr, weight_mtr^lambda, trans_info)[[2]]
          } 

        prob           <- as.numeric(rowSums(prob_mtr))
        prob_isoforms  <-  cbind(trans_info,
                                 isoform_prob = prob)
        return(prob_isoforms)
    }
