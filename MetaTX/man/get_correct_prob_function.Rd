\name{get_correct_prob_function}
\alias{get_correct_prob_function}

\title{
get_correct_prob_function
}
\description{
Calculate probability density distribution for features.
}
\usage{
get_correct_prob_function(num_bin_sum, align_mtr, weight_mtr, trans_info)
}

\arguments{
  \item{num_bin_sum}{
The number of bins into which a RNA model is divided.
}
  \item{methyl}{
Genomic features, which should be a GRanges object.
}
 \item{weight_mtr}{
Weight matrix.
}
 \item{trans_info}{
A data frame object containing methylation sites and associated transcripts.
}
}

\value{
a list object.
}

\examples{

## load remap_results_m6A
data(remap_results_m6A_1)

align_mtr       <- remap_results_m6A_1[[1]]
width_mtr       <- remap_results_m6A_1[[2]]
num_bin_sum     <- ncol(align_mtr)
weight_start    <- num_bin_sum / 5 + 1
weight_end      <- num_bin_sum / 5 * 4
weight_mtr      <- replicate(num_bin_sum, rowSums(width_mtr[, weight_start : weight_end]))    
alpha           <- get_correct_prob_function(num_bin_sum, align_mtr, weight_mtr*weight_mtr, trans_info)[[1]]

}
