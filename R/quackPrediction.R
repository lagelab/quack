#' Quack prediction function
#' 
#' Using the trained model from user uploaded network, using seed genes, predict genes that should be in a network.
#' @param usr_network_name Characters of user's network name, used in make graphs and train model functions.
#' @param usr_seed_file Characters of path and name of user's seed file containing list of genes.
#' @param trained_network_file Characters of path and name to network file from make graphs function.
#' @param quack_gene_model Characters of path and name to quack gene model created from train model function.
#' @param usr_output_loc Characters of location to save output file.
#' @param usr_pathway_file Characters of name of pathway files, which can be found in data/ or user uploaded.
#' @return text file
#' @import hash igraph randomForest
#' @export
#' @examples 
#' quackPrediction("InWeb3", "data/SandersGenes.txt", "data/InWeb3.RData", "QuackGeneModel_Quackv1.3-InWeb3-General853StringentPathways.RData", "data/", "data/stringent853Pathways.RData")
quackPrediction <- function(usr_network_name, usr_seed_file, trained_network_file, quack_gene_model, usr_output_loc, usr_pathway_file){
  
  ### Load user's trained network file (from makeGraphs.R)
  load(file=(trained_network_file))
  ### Stored in quack::data() OR upload user's pathway file
  load(usr_pathway_file)
  
  pathwayName <- strsplit(basename(usr_pathway_file), "\\.")[[1]][1]
  pathwayClassKey <- "General"
  
  quackVersion <- "Quackv1.3"
  predictorList <- c("degreeRatio","wdegreeRatio","ecRatio", "btwnRatio", "closRatio", "lccRatio","giEigenv","giDegree","giWDegree","giCloseness", 
                     "giBetweenness","giLCC","giDegreeInG","giWDegreeInG","giECInG","giClosenessInG","giBetweenessInG","giLCCInG")
  
  ### Change  network name with user input
  networkKey <- usr_network_name
  modelSpecs <- paste(quackVersion, "-",networkKey,"-",pathwayClassKey,pathwayName,sep="")
  
  ### Input user's seed gene lists
  seed <- scan(usr_seed_file, what = character())
  ### Set seed.name with user input
  seed.name <- strsplit(basename(usr_seed_file), "\\.")[[1]][1]
  
  network.name <- networkKey
  
  ### Load quack gene model trained from trainModel.R
  load(quack_gene_model)
  
  g <<- usrNetworkGraphs[[keys(usrNetworkGraphs)[[1]]]]
  remove(usrNetworkGraphs)
  
  predictorList <- c("degreeRatio","wdegreeRatio","ecRatio", "btwnRatio", "closRatio", 
                     "lccRatio","giEigenv","giDegree","giWDegree","giCloseness", 
                     "giBetweenness","giLCC","giDegreeInG","giWDegreeInG","giECInG",
                     "giClosenessInG","giBetweenessInG","giLCCInG")
  
  CalculateModelingFile <- function(pathwayGeneSymbols, graph, QuackModel){
    otherFile <- NULL
    #downSampleSize <- 200
    modelingFile <- data.frame()
    pathwayConnectivityInfo <- getConnectivityInfo(pathwayGeneSymbols)
    thisPathwayData <- computeStatsForAllConnectedNodes(pathwayGeneSymbols,
                                                        downSampleFlag = FALSE,
                                                        usePathwaySize = FALSE,
                                                        multFactor = 0,
                                                        sampleSize = 0)
    
    ### Create obs for the unconnected and unconvered nodes
    nd <- pathwayConnectivityInfo[[4]]$disconnected
    unconnectedDF <- data.frame(gi=rep(0,nd), isInPathway = as.factor(rep(1,nd)))
    rowsUnc <- nrow(unconnectedDF)
    if(is.na(rowsUnc) || is.null(rowsUnc)) rowsUnc <- 0
    
    un <- pathwayConnectivityInfo[[4]]$uncovered
    uncoveredDF <- data.frame(gi=rep(-1,un), isInPathway = as.factor(rep(1,un)))
    rowsUnv <- nrow(uncoveredDF)
    if(is.na(rowsUnv) || is.null(rowsUnv)) rowsUnv <- 0
    
    otherObs <- data.frame()
    if( rowsUnc > 0 && rowsUnv > 0 ){
      #unconnectedDF$pathway <- pathway
      #uncoveredDF$pathway <- pathway
      otherObs <- rbind(unconnectedDF,uncoveredDF)
    }else if ( rowsUnc > 0 && ! (rowsUnv > 0) ){
      #unconnectedDF$pathway <- pathway
      otherObs <- unconnectedDF
    }else if (! (rowsUnc > 0) && rowsUnv > 0 ){
      #uncoveredDF$pathway <- pathway
      otherObs <- uncoveredDF
    }
    
    returnList <- c(list(thisPathwayData),list(otherObs),list(pathwayConnectivityInfo[[4]]))
    
    modelingFile <- data.frame(thisPathwayData)
    rm(thisPathwayData)
    
    if(nrow(otherObs)>0) otherFile <- data.frame(otherObs)
    
    connnectivityInfoFile <- data.frame(pathwayConnectivityInfo[[4]]) 
    
    modelingFile$degreeRatio <- ifelse(modelingFile$giDegreeInG==0,0,modelingFile$giDegree / modelingFile$giDegreeInG)
    modelingFile$wdegreeRatio <- ifelse(modelingFile$giWDegreeInG==0,0,modelingFile$giWDegree / modelingFile$giWDegreeInG)
    
    modelingFile$ecRatio <- ifelse(modelingFile$giECInG==0, 0,modelingFile$giEigenv / modelingFile$giECInG)
    modelingFile$btwnRatio <- ifelse(modelingFile$giBetweenessInG==0, 0,modelingFile$giBetweenness / modelingFile$giBetweenessInG)
    modelingFile$closRatio <- ifelse(modelingFile$giClosenessInG==0, 0,modelingFile$giCloseness / modelingFile$giClosenessInG)
    modelingFile$lccRatio <- ifelse(modelingFile$giLCCInG==0, 0,modelingFile$giLCC / modelingFile$giLCCInG)
    
    modelingFile$ecRatio <- ifelse(modelingFile$ecRatio==0,0,log(modelingFile$ecRatio))
    modelingFile$btwnRatio <- ifelse(modelingFile$btwnRatio==0,0,log(modelingFile$btwnRatio))
    modelingFile$closRatio <- ifelse(modelingFile$closRatio==0,0,log(modelingFile$closRatio))
    modelingFile$lccRatio <- ifelse(modelingFile$lccRatio==0,0,log(modelingFile$lccRatio))
    
    ### Omit missing values and change the type of the response for modeling
    modelingFile$isInPathway <- suppressWarnings(as.factor(modelingFile$isInPathway))
    modelingFile <- na.omit(modelingFile)
    modelingFile$isInPathway <- modelingFile$isInPathway[,drop=TRUE]
    modelingFile$isInPathway <- suppressWarnings(as.factor(modelingFile$isInPathway))
    
    rowsInOther <- nrow(otherFile)
    if( is.null(rowsInOther) || is.na(rowsInOther) ) rowsInOther <- 0
    
    if(rowsInOther > 0){
      
      ### Add predictor columns
      addTheseColumns <- colnames(modelingFile)
      for(col in addTheseColumns){
        if(col!="isInPathway" && col!="gi" && col != "pathway"){
          otherFile[col] <- rep(0,nrow(otherFile))
        }
      }
      modelingFile <- rbind(modelingFile,otherFile) # combine connected and rest
    }
    modelingFile <- modelingFile[sample(nrow(modelingFile)),] #shuffle
    return(modelingFile)
  }
  
  
  computeStatsForAllConnectedNodes <- function(inSymbols,downSampleFlag,usePathwaySize,multFactor,sampleSize){
    
    allIndices <- getGeneIndices(inSymbols)
    connectedIndices <- getConnectedGeneIndices(inSymbols)
    disconnectedIndices <- getDisconnectedGeneIndices(inSymbols)
    
    neighbors <- neighborhood(g, 1, nodes=allIndices) #cant ignore mode mode=c("all", "out", "in") if edges are directed
    neighbors <- unique(unlist(neighbors))
    
    context <- neighbors[ ! neighbors %in% allIndices ]
    contextSub <- c()
    contextSize <- length(context)
    if(!downSampleFlag){
      contextSub <- context
    } else {
      if(usePathwaySize){
        contextSub <- sample(context, length(allIndices)*multFactor)
      }else{
        if ( contextSize > sampleSize ){
          contextSub <- sample(context, sampleSize)
        }else{
          contextSub <- context
        }
      }
    }
  
    ### Connected indices and those that are in context (but not unconnected pathway nodes)
    finalList <- c(contextSub, connectedIndices)
    
    allScoredNodes <- data.frame()
    firstGI <- TRUE
    
    for(gi in finalList){
      
      isInPathway <- 0
      if(gi %in% allIndices){
        isInPathway <- 1
      }
      
      tempPathwayGeneIndices <- unique(c(allIndices,gi))
      
      temppg <- induced.subgraph(g, which(V(g) %in% tempPathwayGeneIndices))
      indexOfgi <- match(gi, V(temppg)$indexInG)
      
      tempdegree <- degree(temppg,loops = FALSE)
      giDegree <- tempdegree[indexOfgi]
      
      tempec <- evcent(temppg)$vector
      giEigenv <- tempec[indexOfgi]
      
      if(giDegree==0) giEigenv <- 0
      
      tempwdegree <- graph.strength(temppg,loops = FALSE)
      giWDegree <- tempwdegree[indexOfgi]
      
      closeness <- closeness(temppg, normalized = FALSE)
      giCloseness <- closeness[indexOfgi]
      
      betweenness <- betweenness(temppg,normalized = FALSE)
      giBetweenness <- betweenness[indexOfgi]
      
      lcc <- transitivity(temppg, type="local") #local clustering coefficient
      lcc[is.na(lcc)] <- 0  #correct NA's with 0
      giLCC <- lcc[indexOfgi]
      
      giDegreeInG     <- get.vertex.attribute(g,"degree",index=gi)
      giWDegreeInG    <- get.vertex.attribute(g,"wDegree",index=gi)
      giClosenessInG  <- get.vertex.attribute(g,"closeness",index=gi)
      giBetweenessInG <- get.vertex.attribute(g,"betweenness",index=gi)
      giLCCInG        <- get.vertex.attribute(g,"lcc",index=gi)
      giECInG         <- get.vertex.attribute(g,"ec",index=gi)
      
      pathwaySize <- vcount(temppg)
      
      df <- data.frame(gi,
                       giEigenv,
                       giDegree,
                       giWDegree,
                       giCloseness,
                       giBetweenness,
                       giLCC,
                       
                       giDegreeInG,
                       giWDegreeInG,
                       giECInG,
                       giClosenessInG,
                       giBetweenessInG,
                       giLCCInG, 
                       
                       pathwaySize,
                       isInPathway)
      
      if(firstGI){
        allScoredNodes <- df
        firstGI <- FALSE
      }else{
        allScoredNodes <- rbind(allScoredNodes, df)   
      }
    }
    return (allScoredNodes)
  }
  
  #####################
  ## Important
  ## Overwrite original
  #####################
  
  splitHoldoutForAUC <- function(holdoutFile, model, response, predList){
    
    tempP <- predict(model,newdata = holdoutFile[,c(predList,response)],type = "prob")
    probIsPathway <- tempP[,2]
    holdoutFile$QuackP <- probIsPathway
    
    ### correct
    holdoutFile$QuackP <- ifelse(holdoutFile$gi < 1,0,holdoutFile$QuackP) #adjust unconnected
    
    con <- subset(holdoutFile, gi>0)
    conAndUnc <- subset(holdoutFile, gi>=0)
    conAndUncAndUncv <- subset(holdoutFile, gi>=-1)
    
    ### connected only
    con <- con[,c("isInPathway", "QuackP", "gi")]
    con <- con[c('gi', 'QuackP','isInPathway')]
    con <- con[order(con$QuackP,decreasing=TRUE),] #sort
    
    conAndUnc <- conAndUnc[,c("isInPathway", "QuackP", "gi")]
    conAndUnc <- conAndUnc[c('gi', 'QuackP','isInPathway')]
    conAndUnc <- conAndUnc[order(conAndUnc$QuackP,decreasing=TRUE),] #sort
    
    conAndUncAndUncv <- conAndUncAndUncv[,c("isInPathway", "QuackP", "gi")]
    conAndUncAndUncv <- conAndUncAndUncv[c('gi', 'QuackP','isInPathway')]
    conAndUncAndUncv <- conAndUncAndUncv[order(conAndUncAndUncv$QuackP,decreasing=TRUE),] #sort
    
    return (c(list(con), list(conAndUnc), list(conAndUncAndUncv)))
  }
  
  #splitFile <- splitTrainingAndHoldout(modelingFile, 0.7, byPathway = T)
  #modelingFile.holdout <- splitFile[[1]]

  #pathways <- c()
  #aucs <- c()
  
  allIndices <- getGeneIndices(seed)
  # allNeighbors <- c()
  # for(index in allIndices){
  #   allNeighbors <- c(allNeighbors, as.numeric(neighbors(g, index)))
  # }
  # allNeighbors <- unique(allNeighbors)
  # allGenes <- c(V(g)$label[allNeighbors], seed)
  this.modelingFile <- CalculateModelingFile(seed)
  
  holdoutsForAUCs <- splitHoldoutForAUC(this.modelingFile, QuackGeneModel, "isInPathway", predictorList) 
  aucFile <- holdoutsForAUCs[[1]]
  
  score.file <- subset(aucFile, !gi == -1)
  score.file$Gene <- V(g)$label[score.file$gi]
  score.file$Seed <- ifelse(score.file$Gene%in%seed, 1, 0)
  score.file <- score.file[order(-score.file$QuackP),]
  score.file <- subset(score.file, Seed == 0 & !QuackP == 0)
  score.file <- score.file[,c('Gene', 'QuackP')]

  ### Save quack predicted genen list as text file
  write.table(score.file, file = paste0(usr_output_loc, seed.name, networkKey, 'QuackPredictions.txt'),
              sep = '\t', row.names = F, col.names = T, quote = F)
  
}
