\name{get_tail_function}
\alias{get_tail_function}

\title{
get_tail_function
}
\description{
Remap genomic features to the tail component of RNA.
}
\usage{
get_tail_function(width_tail, num_bin_tail, threeUTR_align, threeUTR_tx0, methyl)
}

\arguments{
  \item{width_tail}{
Tail length.
}
  \item{num_bin_tail}{
The number of bins into which a tail component is divided.
}
  \item{threeUTR_align}{
Mapped results of features over the 3'UTR component.
}
  \item{threeUTR_tx0}{
3'UTR-related annotation.
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
threeUTR_tx0_1  <- threeUTRsByTranscript(txdb,use.names=FALSE)

## read genomic features:
data("m6A_methyl_1")

trans                  <- mapToTranscripts(m6A_methyl_1, txdb)
trans_info_1           <- get_trans_function(trans, m6A_methyl_1)
threeUTR_align_1        <- get_align_function(10, 
                                             trans_info = trans_info_1, 
                                             tx0 = threeUTR_tx0_1)

tail_align             <- get_tail_function(width_tail = 1000, 
                          num_bin_tail = 10, 
                          threeUTR_align = threeUTR_align_1, 
                          threeUTR_tx0 = threeUTR_tx0_1, 
                          methyl = m6A_methyl_1)
    
}
