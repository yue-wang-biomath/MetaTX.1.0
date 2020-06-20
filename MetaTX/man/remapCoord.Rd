\name{remapCoord}
\alias{remapCoord}

\title{
remapCoord
}
\description{

Remap genomic features represented by genome-based coordinate to an mRNA model. An alignment matrix, a width matrix and other associated annotations can be returned.

}
\usage{
remapCoord(site_info, txdb = NA, num_bin_fiveUTR = 10, num_bin_cds = 10, num_bin_threeUTR = 10, num_bin_promoter = 10, num_bin_tail = 10, includeNeighborDNA = FALSE, cds_by_tx0 = NA, fiveUTR_tx0 = NA, threeUTR_tx0 = NA)
}

\arguments{
  \item{features}{
Genomic features, which should be a GRanges object.
}
  \item{txdb}{
A TxDb object, which should be a container for storing transcript annotations.
}
  \item{num_bin}{
The number of bins into which an mRNA component is divided.
}
  \item{includeNeighborDNA}{
Whether neighborhood DNA regions should be included.
}
  \item{cds_by_tx0}{
CDS-related annotation.
}
  \item{fiveUTR_tx0}{
5'UTR-related annotation.
}
  \item{threeUTR_tx0}{
3'UTR-related annotation.
}
}

\details{
This function returns a list contains three components, which are an alignment matrix, a width matrix, and a matrix of other associated annotations.

}
\value{
A list object. 
}

\examples{
## download the txdb object from UCSC
txdb           <- TxDb.Hsapiens.UCSC.hg19.knownGene
cds_by_tx0_1   <- cdsBy(txdb, "tx")
fiveUTR_tx0_1  <- fiveUTRsByTranscript(txdb,use.names=FALSE)
threeUTR_tx0_1 <- threeUTRsByTranscript(txdb,use.names=FALSE)

## read genomic features:
data("m6A_methyl_1")

## remap genomic features to an mRNA model:
remap_results_m6A_1 <- remapCoord(features = m6A_methyl_1, txdb = txdb, num_bin = 10,                                        includeNeighborDNA = TRUE,
                                  cds_by_tx0         = cds_by_tx0_1, 
                                  fiveUTR_tx0        = fiveUTR_tx0_1,
                                  threeUTR_tx0       = threeUTR_tx0_1) 
}
