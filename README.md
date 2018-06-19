# Quack (improved performance with parallel computing)

## Overview ##
Quack is a machine learning algorithm to predict biological pathway relationships based on functional genomics networks:

* `makeGraph()` creates compatible graph network file using user uploaded network or InWeb3 network file
* `trainModel()` trains pathway models using Quack
* `quackPrediction()` predicts potential candidate genes in a network based on user input genes.

Quack requires R and RStudio on either Mac or Windows. For best results, use R version 3.3.2.<br>
R is available at https://www.r-project.org/<br>
RStudio is available at https://www.rstudio.com/products/rstudio/download/<br>
No additional hardware is required.


## Installation ##
```r
install.packages("devtools")
devtools::install_github("lagelab/quack", ref = "ranger")
```
Depending on your internet speed and packages already installed, the installation may take up to 10 minutes. 

The Quack machine learning method can also be accessed through the Broad Institute Web Platfor for Genome Networks ([GeNets](https://apps.broadinstitute.org/genets)) 

Quack has been tested on and is stable using the following version of each packages:
```r
R 3.3.2
boot 1.3-20
changepoint 2.2.2
data.table 1.10.4-3
doMC 1.3.5
doSNOW 1.0.16
doParallel 1.0.11
dplyr 0.7.4
entropy 1.2.1
hash 2.2.6
igraph 1.1.2
lattice 0.20-35
MASS 7.3-48
mclust 5.4
poweRlaw 0.70.1
pROC 1.10.0
randomForest 4.6-12
ranger 0.9
stringr 1.2.0
```

## Demo ##
If the user would like to predict Quack candidate genes using the default InWeb3 network and 853 pathways, user may skip ahead to `quackPrediction()` function. Please download the Quack gene model file found <a href="http://www.lagelab.org/wp-content/uploads/2018/02/QuackGeneModel_Quackv1.3-InWeb3-General853StringentPathways.RData_.zip">here</a> and unzip.

### makeGraphs() ### 
`makeGraphs()` function creates Quack compatible graph network file. `makeGraphs` function takes three arguments: `usr_network_file`, `usr_network_name`, and `usr_save_loc`. `usr_network_file` is the path to user's network file, which user may upload a tab delimited network file of three columns representing source, target and numeric score associated with each edge (do not include headers) or use the InWeb3 network file in `data` folder. `user_network_name` is characters referring to the user defined name of the network. `usr_save_loc` is the directory where the network output will be saved.

`makeGraphs()` example:
```r
makeGraphs("~/data/InWeb3_network.txt", "InWeb3", "~/data/output")
```

The runtime varies depending on the size of your network. This will generate an RData file in the name of `usr_network_name` argument.

### trainModel() ###
`trainModel()` function trains pathway model using the user defined network and known pathways. This function takes three arguments: `usr_modeling_file`, `usr_pathway_file`, and `usr_output_loc`. `usr_modeling_file` is the path to network file generated from `makeGraphs()` function. `usr_pathway_file` is a user defined pathway file, which is in the form of biological pathway name in first column, and gene names identified in a given pathway in the subsequent columns with no headers. User may also use the provided stringent853Pathways.RData file in the `data` folder. `usr_output_loc` is the directory where the trained pathway model will be saved.

`trainModel()` example:
```r
trainModel("~/data/output/InWeb3.RData", "~/data/stringent853Pathways.RData", "~/data/output")
```

The runtime varies depending on the size of the network and the number of pathways. Parallel computing is highly recommended, and the runtime may take more than 12 hours on a normal computer. This will generate an RData file in the name of QuackGeneModel_Quackv1.3 - `usr_network_name` - General - `usr_pathway_file`.

### quackPrediction() ###
`quackPrediction()` function uses the information from network and trained pathway model to predict potential candidate genes in a network based on user input genes. This function takes six arguments: `usr_network_name`, `usr_seed_file`, `trained_network_file`, `quack_gene_model`, `usr_output_loc`, `usr_pathway_file`. `usr_network_name` is the network name used in `makeGraphs()` function. `usr_seed_file` is the path to user's file containing list of genes with no headers. `trained_network_file` is the output from `makeGraphs()` function. `quack_gene_model` is the output from `trainModel()` function. `usr_output_loc` is the directory where the predicted genes and Quack p-value will be saved, and lastly `usr_pathway_file` is the pathway file used in `trainModel()` function.

`quackPrediction()` example:
```r
quackPrediction("InWeb3", "~/data/SandersGenes.txt", "~/data/output/InWeb3.RData", "QuackGeneModel_Quackv1.3-InWeb3-General853StringentPathways.RData", "~/data/output", "~/data/stringent853Pathways.RData")
```

The runtime varies depending on the size of seed genes and may exceed 10 minutes. This will generate a text file in the name of `usr_seed_file` - `usr_network_name` - QuackPredictions.
