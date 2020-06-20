\name{get_promoter_function}
\alias{get_promoter_function}

\title{
get_promoter_function
}
\description{
Remap genomic features to the promoter component of RNA.
}
\usage{
get_promoter_function(width_promoter, num_bin_promoter, fiveUTR_align, fiveUTR_tx0, methyl)
}

\arguments{
  \item{width_promoter}{
Promoter length.
}
  \item{num_bin_promoter}{
The number of bins into which a promoter component is divided.
}
  \item{fiveUTR_align}{
Mapped results of features over the 5'UTR component.
}
  \item{fiveUTR_tx0}{
5'UTR-related annotation.
}
  \item{methyl}{
Genomic features, which should be a GRanges object.
}
}

\value{
A data frame object.
}

\examples{
## download the txdb object from UCSC
txdb           <- TxDb.Hsapiens.UCSC.hg19.knownGene
fiveUTR_tx0_1  <- fiveUTRsByTranscript(txdb,use.names=FALSE)

## read genomic features:
data("m6A_methyl_1")

trans                  <- mapToTranscripts(m6A_methyl_1, txdb)
trans_info_1           <- get_trans_function(trans, m6A_methyl_1)
fiveUTR_align_1        <- get_align_function(10, 
                                             trans_info = trans_info_1, 
                                             tx0 = fiveUTR_tx0_1)
promoter_align         <- get_promoter_function(width_promoter = 1000, 
                          num_bin_promoter = 10, 
                          fiveUTR_align = fiveUTR_align_1, 
                          fiveUTR_tx0 = fiveUTR_tx0_1, 
                          methyl = m6A_methyl_1)
   
}
