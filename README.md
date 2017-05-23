# Quack

## Overview ##
Quack is a machine learning algorithm to predict biological pathway relationships based on functional genomics networks:

* `makeGraph()` creates compatible graph network file using user uploaded network or InWeb3 network file
* `trainModel()` trains pathway models using Quack
* `quackPrediction()` predicts potential candidate genes in a network based on user input genes.

## Installation ##
```r
devtools::install_github("lagelab/quack")
```

The Quack machine learning method can also be accessed through the Broad Institute Web Platfor for Genome Networks ([GeNets](https://apps.broadinstitute.org/genets)) 
