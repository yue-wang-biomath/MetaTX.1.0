\name{get_align_function}
\alias{get_align_function}

\title{
get_align_function
}
\description{
Remap genomic features to 5'UTR/CDS/3'UTR RNA component.
}
\usage{
get_align_function(num_bin, trans_info, cds_by_tx0)
}
\arguments{
  \item{num_bin}{
The number of bins into which an mRNA component is divided.
}
  \item{trans_info}{
A data frame object containing methylation sites and associated transcripts.
}
  \item{tx0}{
5'UTR/CDS/3'UTR-related annotation.
}
}

\value{
A data frame object.
}

\examples{
## download the txdb object from UCSC
txdb           <- TxDb.Hsapiens.UCSC.hg19.knownGene
cds_by_tx0_1   <- cdsBy(txdb, "tx")

## read genomic features:
data("m6A_methyl_1")

trans                  <- mapToTranscripts(m6A_methyl_1, txdb)
trans_info_1           <- get_trans_function(trans, m6A_methyl_1)
cds_align              <- get_align_function(num_bin = 10, 
                                             trans_info = trans_info_1, 
                                             tx0 = cds_by_tx0_1)

}
