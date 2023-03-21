#' Map genomic features/regions to an mRNA model.
#' @export remapCoord
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @importFrom GenomicFeatures mapToTranscripts
#' @importFrom GenomicFeatures cdsBy fiveUTRsByTranscript threeUTRsByTranscript
#' @description Map genomic features/regions to an mRNA model. An mRNA model contains three main components: 5UTR, CDS, 3UTR (probably promoter and tail). It returns a remap.results (list) object contains an alignment matrix, a width matrix and an annotation matrix.
#' @usage remapCoord(features,txdb = NA, num_bin = 10, includeNeighborDNA = TRUE)
#' @param features Genomic features of interest, which should be a \code{GRanges} object. Each range should be single-based.
#' @param txdb A TxDb object.
#' @param num_bin The number of bins each mRNA component (5UTR, CDS or 3UTR) is divided into.
#' @param includeNeighborDNA Whether include neighborhood regions (promoter and tail).
#' @return A remap.results (list) object.
#' \itemize{
#' \item \bold{\code{align_mtr}} Alignment matrix. Each row represents the mapping result of a particular feature on a particular transcript. The information of this feature and this transcript can be found on the same row in the \code{trans_info}.
#' \item \bold{\code{width_mtr}} Width matrix. Each row represents the width of different bins on a particular transcript.
#' \item \bold{\code{trans_info}} Annotation data.frame.
#' }
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' file <- system.file(package="MetaTX", "extdata/m6A_methyl.rds")
#' m6A_methyl<- readRDS(file)
#' remap_results_m6A_1 <- remapCoord(features = m6A_methyl[1:100], txdb = txdb, num_bin = 10)
remapCoord <-
  function(features,
           txdb               = NA,
           num_bin            = 10,
           includeNeighborDNA = TRUE){

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
    trans                  <- mapToTranscripts(features, txdb)
    trans_info             <- get_trans_function(trans, features)

    # get annotation
    cds_by_tx0   <- cdsBy(txdb, "tx")
    fiveUTR_tx0  <- fiveUTRsByTranscript(txdb,use.names=FALSE)
    threeUTR_tx0 <- threeUTRsByTranscript(txdb,use.names=FALSE)

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

    # get width matrix
    width_mtr_5cds3        <- data.frame(coordinate = cbind(replicate(num_bin_fiveUTR, fiveUTR_align[, 'width'] / num_bin_fiveUTR)
                                                            , replicate(num_bin_cds, cds_align[, 'width'] / num_bin_cds)
                                                            , replicate(num_bin_threeUTR, threeUTR_align[, 'width'] / num_bin_threeUTR)))

    if(includeNeighborDNA){
      # get promoter and tail
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
    remap_results <- list(align_mtr, width_mtr, trans_info)
    class(remap_results) <- "remap.results"
    names(remap_results) <- c("align_mtr", "width_mtr", "trans_info")
    return(remap_results)
  }
