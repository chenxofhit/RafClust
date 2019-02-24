# RafClust

## Introduction

A novel random forest  similarity learning based clustering method for single cell RNA sequencing data

Author: Xiang Chen [au]; Min Li [cre]

Maintainer: Xiang Chen <chenxofhit@gmail.com>

## Installation

```R
library(devtools)
install_github("chenxofhit/RafClust")
```

If [devtools](https://github.com/hadley/devtools) package is not avaliable, please invoke R and then type the following command first: 

```R
install.packages("devtools")
```

## Usage

```R
res <- RafClust(texpr,verbose = TRUE) #cluster number is determined by the program
res <- RafClust(texpr,NumC = 10)
```
