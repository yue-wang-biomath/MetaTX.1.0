\name{MetaTX-package}
\alias{MetaTX-package}
\alias{MetaTX}
\docType{package}
\title{
\packageTitle{MetaTX}
}
\description{
\packageDescription{MetaTX}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{MetaTX}
\packageIndices{MetaTX}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
\packageAuthor{MetaTX}

Maintainer: \packageMaintainer{MetaTX}
}
\

\examples{
## load remap_results_m6A
data(remap_results_m6A_1)

## download the txdb object from UCSC
txdb  <- TxDb.Hsapiens.UCSC.hg19.knownGene

## plot
p1 <- metaTXplot(remap_results_m6A_1,
                 num_bin              = 10,
                 includeNeighborDNA   = TRUE,
                 relativeProportion   = c(1, 1, 1, 1),
                 lambda = 2,
                 adjust = 0.15,
                 title  = '',
                 legend = '',
                 returnPlot = TRUE,
                 type = 'relative')
p1

}
