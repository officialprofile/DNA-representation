# GNR-DNA
Graphical and Numerical Representation of DNA Sequences

R library for characterizing biological sequences graphically and numerically. If you would like to get the gist of implemented approach please the following article 
> A new graphical representation and analysis of DNA sequence structure: 1. Methodology and application to globin genes, Curr Sci 66:309-314

## Install
In order to install and load the package run the following code in the R console

```r
library(devtools)
install_github('officialprofile/GNR-DNA')
library(drep)
```

## How to use it
In most cases a single line of code can yield satisfactory results. For example

```r
plot2DGraph(c('KX369547', 'HQ234498'), genbank = TRUE, main = 'Comparison of two sequences')
```
returns a ready-to-use graph of the ZIKV genomes

<img src="img/example1.png" width="50%" />

```r
plot2DGraph(c('KX369547', 'HQ234498', 'MH063265'), genbank = TRUE, main = 'Comparison of three sequences')
```

<img src="img/example2.png" width="50%" />

Notice that plotting such graphs can be achieved solely by putting the GenBank accession number.

In similar fashion one can obtain numerical characteristics by employing dRep function, e.g.
```r
dRep(c('KX369547', 'HQ234498'), genbank = TRUE)
```
returns the following dataframe
|         |len  |mi_x  |mi_y   |sqrt   | I_xx    | I_yy    |  I_xy    |...
|---------|-----|------|-------|-------|---------|---------|----------|---
|KX369547 |10769|84.660|-16.061| 86.170| 11371982| 60032538|-17116036 |...
|HQ234498 |10269|75.171|-17.691| 77.224|  5137647| 31072416|  -1928373|...

Naturally, instead of using data from GenBank, one can apply implemented method to one's own sequence or vector of sequences, e.g.
```r
seq <- 'ACCCTCGCGCCGCGATTCTACGGACCCTGAAAATG'
dRep(seq)
```
