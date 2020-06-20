# MetaTX.1.0
## Introduction
The MetaTX is aimed for plotting the transcriptomic distribution of RNA-related genomic features.


## 1. Quick Start with MetaTX

To install MetaTX from Github, please use the following codes.

```{r introduction}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GenomicAlignments", "GenomicRanges", "GenomicFeatures", "ggplot2",
                       "TxDb.Hsapiens.UCSC.hg19.knownGene", "cowplot"))
               
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("yue-wang-biomath/MetaTX/MetaTX.1.0")
library('MetaTX')
```
Or you can download directly from ```MetaTX_1.0.tar.gz```. 

## 2. Data preprocessing 

First, the TxDb object and other information about mRNA components need to be downloaded.

```
txdb           <- TxDb.Hsapiens.UCSC.hg19.knownGene
cds_by_tx0_1   <- cdsBy(txdb, "tx")
fiveUTR_tx0_1  <- fiveUTRsByTranscript(txdb,use.names=FALSE)
threeUTR_tx0_1 <- threeUTRsByTranscript(txdb,use.names=FALSE)
```
It requires basic information of target feature set, involving the genomic locations, seqnames and strand types of each feature. The input feature set is required to be provided as a GRanges object.

In the MetaTX package, we provide three example feature sets stored in the file ```m6A_methyl_1.rda```, ```m6A_methyl_2.rda``` and ```m6A_methyl_3.rda```. They are m6A datasets derived from different high-throughput sequencing approaches, including an miCLIP-seq dataset (Linder, et al., 2015; Olarerin-George and Jaffrey, 2017), a PA-m6A-seq dataset (Chen, et al., 2015) and an m6A-seq dataset (Schwartz, et al., 2014). 

We will use the example ```m6A_methyl_1.rda``` to illustrate how to use MetaTX sketching feature distribution. Users can also use other RNA-related genomic feature datasets.

Load the example dataset provided in the MetaTX package. 
```
data("m6A_methyl_1")
```

Please see the following example, which will read m6A methylation sites from the file ```m6A_methyl_1.rda``` into R and map these features to an mRNA model. 


```
remap_results_m6A_1 <- remapCoord(features = m6A_methyl_1, txdb = txdb, num_bin = 10, includeNeighborDNA = TRUE,
                                  cds_by_tx0         = cds_by_tx0_1, 
                                  fiveUTR_tx0        = fiveUTR_tx0_1,
                                  threeUTR_tx0       = threeUTR_tx0_1) 
``` 

We provide this result and store it in the file ```remap_results_m6A_1.rda```.

## 3. Visualization of the transcriptomic distribution 

Use previously generated results or provided file ```remap_results_m6A_1.rda```.

```
data("remap_results_m6A_1")
```
The ```metaTXplot``` function enables the visualization of RNA-related genomic features. 

As shown in the following codes, m6A pattern is visualized with the promoter/5’UTR/CDS/3’UTR/tail ratio of 1:1:1:1:1 (a, c) and 3:1:3:2:3 (b, d). MetaTX R package supports two different ways of visualization. The ‘absolute’ method (a, b) provides absolute density (with the unit: number of features per bp exon transcript), which will not be affected by the relative length of different RNA components defined by the user. The ‘relative’ method (c, d) provides probability density function (with the area under the curve equals to 1), which can be affected by the relative length of different RNA components specified by user.
```
p1 <-  metaTXplot(remap_results_m6A_1,
                 num_bin              = 10,
                 includeNeighborDNA   = TRUE,
                 relativeProportion   = c(1, 1, 1, 1),
                 title  = '(a)',
                 legend = 'absolute',
                 type = 'absolute'
    )

p2 <-  metaTXplot(remap_results_m6A_1,
                 num_bin              = 10,
                 includeNeighborDNA   = TRUE,
                 relativeProportion   = c(1, 3, 2, 3),
                 title  = '(b)',
                 legend = 'absolute',
                 type = 'absolute'
    )

p3 <-  metaTXplot(remap_results_m6A_1,
                 num_bin              = 10,
                 includeNeighborDNA   = TRUE,
                 relativeProportion   = c(1, 1, 1, 1),
                 title  = '(c)',
                 legend = 'relative',
                 type = 'relative'
    )

p4 <-  metaTXplot(remap_results_m6A_1,
                 num_bin              = 10,
                 includeNeighborDNA   = TRUE,
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

![image](https://github.com/yue-wang-biomath/MetaTX/blob/master/Fig.png)

## 4. Resolving ambiguity problem

The package also provides an ```isoformProb``` function that can return the probabilities of a particular feature being located on different isoforms. 

```
# load remap_results_m6A_1
data(remap_results_m6A_1)

isoform_probs <- isoformProb(remap_results_m6A_1, num_bin = 10, includeNeighborDNA = TRUE, lambda = 2)

```

The probabilities of a particular feature being located on different isoforms (the last column) can be returned.

```
      index_trans index_methyl seqnames methyl_pos strand trans_ID isoform_prob
1             1            1    chr19     581474      +    65776 2.022115e-01
2             2            1    chr19     581474      +    65777 2.022115e-01
3             3            1    chr19     581474      +    65778 2.030733e-01
4             4            1    chr19     581474      +    65779 1.894304e-01
5             5            1    chr19     581474      +    65780 2.030733e-01
6             6            2     chr6  122744806      +    25251 5.000000e-01
7             7            2     chr6  122744806      +    25252 5.000000e-01
8             8            3     chr3  195250580      -    17278 1.000000e+00
```


# References
Chen, K. et al. (2014) High-resolution N6-methyladenosine (m6A) map using photo-crosslinking-assisted m6A sequencing. Angew. Chem. Int. Ed., 54, 1587-1590. (https://doi.org/10.1002/anie.201410647)

Linder,B. et al. (2015) Single-nucleotide-resolution mapping of m6A and m6Am throughout the transcriptome. Nat. Methods, 12, 767–772. (https://doi.org/10.1038/nmeth.3453)

Olarerin-George, A.O. and Jaffrey, S.R. MetaPlotR: a Perl/R pipeline for plotting metagenes of nucleotide modifications and other transcriptomic sites. Bioinformatics 2017;33(10):1563-1564 (https://doi.org/10.1093/bioinformatics/btx002)

Schwartz, S. et al. (2014) Perturbation of m6A writers reveals two distinct classes of mRNA methylation at internal and 5'sites. Cell Reports, 8(1), 284-296. (https://doi.org/10.1016/j.celrep.2014.05.048)
