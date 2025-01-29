# ZINQ-L: A Zero-Inflated Quantile Approach for Differential Abundance Analysis of Longitudinal Microbiome Data
The manuscript is publicly available at https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2024.1494401/full.

## Instructions for use

From an `R` session, install `ZINQ-L` by:
```
devtools::install_github('https://github.com/AlbertSL98/ZINQ-L')
```
And find the vignettes at: https://albertsl98.github.io/ZINQ-L-Tutorial/ .

To include the vignettes during installing the package:
```
devtools::install_github('https://github.com/AlbertSL98/ZINQ-L', build_vignettes = TRUE)
```
To view the vignettes, from the `R` terminal, type: 
```
browseVignettes("ZINQL")
```
or
```
vignette("ZINQL_Tutorial", package = "ZINQL")
```

From an `R` session, library the package by:
```
library(ZINQL)
```


Details can be found in the manual: https://github.com/AlbertSL98/ZINQ-L/blob/main/ZINQL.pdf.
