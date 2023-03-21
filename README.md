# MetaTX.1.0
## Introduction

MetaTX is for visualizing the distribution of RNA-related genomic features at the mRNA-level. The density distribution is an EM solution that corrects isoform ambiguity among mRNAs.

- The issue related to installation has been solved. (2023-03-21)

## 1. Install

To install MetaTX from Github, please use the following codes.

```{r introduction}
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("yue-wang-biomath/MetaTX.1.0")
library('MetaTX')
```

## 2. A quick start
### - remapCoord
We use `m6A_methyl.rds`, derived from the miCLIP-seq dataset (Linder, et al., 2015) to illustrate how to use MetaTX to sketch feature distribution. We take the first 1000 m6A sites as an example. Users can also apply it to visualize other RNA-related genomic feature datasets.

First, use `remapCoord` function in MetaTX to map the features to a specified transcriptome, which requires a feature set in `GRanges` format and a txdb object. 

```R
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
file <- system.file(package="MetaTX", "extdata/m6A_methyl.rds")
m6A_methyl<- readRDS(file)
remap_results_m6A <- remapCoord(features = m6A_methyl[1:1000], txdb = txdb)
``` 

The  `remapCoord` function returns a `remap.results` (list) object, which contains an Alignment matrix, a Width matrix and an Annotation data.frame.

### - metaTXplot

Next, `metaTXplot` can calculate and visualize a density distribution of genomic features on an mRNA model. The density distribution it returns is an EM solution that corrects isoform ambiguity among mRNAs.

```R
p1 <- metaTXplot(remap_results_m6A)
```
<img src = 'https://github.com/yue-wang-biomath/MetaTX.1.0/blob/master/inst/extdata/Figures/Figure1.jpg' width = '500px'>

### - isoformProb

The package also provides an `isoformProb` function that can return the probability of a particular feature being located on different isoforms. 

```R
isoform_probs <- isoformProb(remap_results_m6A)
```


## 3. More about visualization

`metaTXplot` supports two different ways of visualization. The ‘absolute’ method (a, b) provides absolute density (with the unit: number of features per bp exon transcript), which will not be affected by the relative length of different RNA components defined by the user. The ‘relative’ method (c, d) provides probability density function (with the area under the curve equals to 1), which can be affected by the relative length of different RNA components specified by user. Furthermore, the ratio of each mRNA component can be adjusted. The following example shows m6A pattern is visualized with the promoter/5’UTR/CDS/3’UTR/tail ratio of 1:1:1:1:1 (a, c) and 3:1:3:2:3 (b, d).
```
p1 <-  metaTXplot(remap_results_m6A,
                 relativeProportion   = c(1, 1, 1, 1),
                 title  = '(a)',
                 legend = 'absolute',
                 type = 'absolute'
    )

p2 <-  metaTXplot(remap_results_m6A,
                 relativeProportion   = c(1, 3, 2, 3),
                 title  = '(b)',
                 legend = 'absolute',
                 type = 'absolute'
    )

p3 <-  metaTXplot(remap_results_m6A,
                 relativeProportion   = c(1, 1, 1, 1),
                 title  = '(c)',
                 legend = 'relative',
                 type = 'relative'
    )

p4 <-  metaTXplot(remap_results_m6A,
                 relativeProportion   = c(1, 3, 2, 3),
                 title  = '(d)',
                 legend = 'relative',
                 type = 'relative'
    )

ggdraw() +
    draw_plot(p1, 0, .5, .5, .5) +
    draw_plot(p2, .5, .5, .5, .5) +
    draw_plot(p3, 0, 0, .5, .5) +     
    draw_plot(p4, .5, 0, .5, .5) 
``` 

<img src = 'https://github.com/yue-wang-biomath/MetaTX.1.0/blob/master/inst/extdata/Figures/Figure.png' width = '1000px'>


## 4. Why should we correct isoform ambiguity?
Lets illustrate why we should correct isoform ambiguity by an EM solution.
To generate simulated data, first we need to install a RgnTX package.
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("RgnTX")
library(RgnTX)
```
We use the following codes to randomly pick 2000 transcripts and pick a site on the CDS region of each transcript.
```R
set.seed(01231)
cds.tx0 <- cdsBy(txdb)
trans.ids <- names(cds.tx0)
regions <- cds.tx0[sample(trans.ids,2000)]
sites_cds <- randomizeTransByOrder(regions, random_length = 1)
sites_cds <- GRangesList2GRanges(sites_cds)
```
If we directly plot the mapping results (counts/width of each bin) of these CDS sites, we can see stronger bias can be observed in the results of this direct estimation.
```R
remap_results_1 <- remapCoord(features = sites_cds, txdb = txdb)
p2 <-  directPlot(remap_results_1)
```
<img src = 'https://github.com/yue-wang-biomath/MetaTX.1.0/blob/master/inst/extdata/Figures/Figure2.jpg' width = '500px'>

We can compare with the real distribution of these CDS sites using the following codes. (Users can skip the details.) 
```R
align_mtr          <- remap_results_1[[1]]
width_mtr          <- remap_results_1[[2]]
trans_info         <- remap_results_1[[3]]
num_bin_sum        <- ncol(align_mtr)
trans_id_real <- sites_cds$transcriptsHits
index_real <- unlist(lapply(1:max(trans_info[, 'index_feature']), function(x){
  index_trans_x <- which(trans_info[, 'index_feature'] == x)
  which_x <- which(trans_info[index_trans_x, 'trans_ID'] == trans_id_real[x])
  return(index_trans_x[which_x])
}))
align_mtr          <- align_mtr[index_real, ]
width_mtr          <- width_mtr[index_real, ]
trans_info         <- trans_info[index_real, ]
remap_results_real <- list(align_mtr, width_mtr, trans_info)
p3 <-  directPlot(remap_results_real)
```
<img src = 'https://github.com/yue-wang-biomath/MetaTX.1.0/blob/master/inst/extdata/Figures/Figure3.jpg' width = '500px'>

Now we correct the bias by the MetaTX model (maximizing the likelihood of each site being located on different isoforms, solved by EM algorithm), we would see the obtained result is much more closer to the real distribution.
```R
p4 <-  metaTXplot(remap_results_1)
```
<img src = 'https://github.com/yue-wang-biomath/MetaTX.1.0/blob/master/inst/extdata/Figures/Figure4.jpg' width = '500px'>

# References
Wang, Y. et al. (2021) MetaTX: deciphering the distribution of mRNA-related features in the presence of isoform ambiguity, with applications in epitranscriptome analysis. Bioinformatics. 37(9), 1285–1291. (https://doi.org/10.1093/bioinformatics/btaa938)

Linder,B. et al. (2015) Single-nucleotide-resolution mapping of m6A and m6Am throughout the transcriptome. Nat. Methods, 12, 767–772. (https://doi.org/10.1038/nmeth.3453)


