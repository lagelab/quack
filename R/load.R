#' @export
getCustomSampleForTrainingPathways <- function(phash,inSampleList,restrictionType,sampleSize, doDownSample, downSampleSize, pathwayBounds, sortByConnected){
  #phash <- useThisPathwayHash
  #inSampleList <- pl
  # restrictionType = "restrictTo"
  # sampleSize = maxFileSize
  #doDownSample = TRUE
  #downSampleSize = 1500
  #pathwayBounds = c(10,100)
  #sortByConnected =  FALSE

  infoFile <- getConnectivityInfoForPathwayList(phash)

  #subset if a list is provided
  nop <- length(inSampleList)
  if(nop > 0){
    if(restrictionType == "restrictTo"){
      infoFile <- infoFile[infoFile$pathway %in% inSampleList,]
    }else{
      infoFile <- infoFile[ !infoFile$pathway %in% inSampleList,]
    }
  }

  #infoFile <- subset(infoFile, connectedPerc > connPercThreshold) # has to have something!
  infoFile <- subset(infoFile, connected > pathwayBounds[1] & connected <= pathwayBounds[2] ) # subset to max pathway size criterion

  if(sortByConnected){
    infoFile <- infoFile[order(infoFile$connectedPerc,decreasing=TRUE),]
  }else{
    infoFile <- infoFile[sample(nrow(infoFile)),] #shuffle
  }

  totalObs <- sum(infoFile$context) + sum(infoFile$total)

  infoFile$gsPlusContext <- infoFile$context + infoFile$total

  if(totalObs <= sampleSize){

    infoFile$takeThisContext <- infoFile$context
    return (infoFile)

  }else{

    if(doDownSample){

      # split to genesets with less than max amount versus more context than max amount
      infoFile.below <- subset(infoFile,context <= downSampleSize)
      infoFile.above <- subset(infoFile,context  > downSampleSize)

      np <- nrow(infoFile.below)
      totalSample <- 0
      finalPathwayList <- NULL
      foundEnough <- FALSE
      if(np > 0){

        i <- 1
        while(!foundEnough && i <= np){

          thisRow <- infoFile.below[i,]
          addThis <- thisRow$context + thisRow$connected

          totalSample <- totalSample + addThis

          thisRow$takeThisContext <- thisRow$context

          if( totalSample >= sampleSize ){
            foundEnough <- TRUE
          }

          finalPathwayList <- rbind(finalPathwayList, thisRow)
          i <- i + 1
        }
      }

      if(foundEnough){
        return (finalPathwayList)
      }else{

        ####################################################################
        ## If the sample size isn't met, get getThisManyMore from contexts #
        ####################################################################
        getThisManyMore <- sampleSize - totalSample

        # sample infoFile.above
        infoFile.above$contextCap <- downSampleSize

        infoFile.above$gsPlusContextCap <- infoFile.above$connected + infoFile.above$contextCap
        totalAdditional <- sum(infoFile.above$gsPlusContextCap)

        if(totalAdditional <= getThisManyMore){
          infoFile.above$takeThisContext <- infoFile.above$contextCap
          infoFile.above <- subset(infoFile.above, select=-c(gsPlusContextCap,contextCap))
          finalPathwayList <- rbind(finalPathwayList, infoFile.above)
        }else{

          increase <- 100
          while(totalAdditional <= getThisManyMore){

            increaseDownSampleSize <- downSampleSize + increase
            infoFile.above$contextCap <- ifelse(infoFile.above$context > increaseDownSampleSize , increaseDownSampleSize, infoFile.above$context)

            infoFile.above$gsPlusContextCap <- infoFile.above$connected + infoFile.above$contextCap
            totalAdditional <- sum(infoFile.above$gsPlusContextCap)

            increase <- increase + 100

          }

          # compute cumulative sum and subset of above's until I have enough
          infoFile.above <- within(infoFile.above, gsPlusContextCapCUM <- cumsum(gsPlusContextCap))
          infoFile.above.sub <- subset(infoFile.above, gsPlusContextCapCUM <= getThisManyMore)

          infoFile.above.sub$takeThisContext <- infoFile.above.sub$contextCap
          infoFile.above.sub <- subset(infoFile.above.sub, select=-c(gsPlusContextCap, gsPlusContextCapCUM,contextCap))

          finalPathwayList <- rbind(finalPathwayList, infoFile.above.sub)

        }

        return (finalPathwayList)
      }

    }else{  # if ! doDownSample

      totalSample <- 0
      finalPathwayList <- NULL
      for(i in 1:nrow(infoFile)){

        thisRow <- infoFile[i,]

        thisRow$takeThisContext <- thisRow$context
        addThis <- thisRow$takeThisContext + thisRow$connected

        totalSample <- totalSample + addThis

        if( totalSample >= sampleSize )
          break

        finalPathwayList <- rbind(finalPathwayList, thisRow)
      }

      return (finalPathwayList)

    }

  } # end else

}

#' @export
getConnectivityInfoForPathwayList <- function(phash){

  PathwayList <- keys(phash)
  upperBound <- length(PathwayList)

#  firstPathway <- TRUE

#  infoFile <- data.frame()

  pb <- txtProgressBar(max = upperBound, style = 3)
  progress <- function(n) { setTxtProgressBar(pb, n) }
  opts <- list(progress = progress)

  r <- foreach(i         = 1:upperBound,
               .combine  = 'rbind',
               .export   = ls(envir=globalenv()),
               .packages = c("hash","igraph"),
               .options.snow = opts) %dopar% {

                   pathway <- PathwayList[i]
                   pathwayGeneSymbols <- phash[[pathway]]

                   pathwayConnectivityInfo <- getConnectivityInfo(pathwayGeneSymbols)

                   dfInfo <- pathwayConnectivityInfo[[4]]
                   dfInfo$pathway <- pathway
                   returnList <-  dfInfo
                 }

  #now append results to data frame
#  for (i in 1:length(r)){
#    if(firstPathway == TRUE){
#      infoFile <-  r[[i]]
#      firstPathway <- FALSE
#    }else{
#      infoFile <- rbind(infoFile, r[[i]])
#    }
#  }

  #infoFile <- infoFile[order(infoFile$connectedPerc,decreasing=TRUE),]
  infoFile <- r[order(r$connectedPerc,decreasing=TRUE),]

  return (infoFile)
}

#' @export
getConnectivityInfo <- function(pathwaySymbols){

  initIndices <- which(V(g)$label %in% pathwaySymbols)

  neighbors <- neighborhood(g, 1, nodes=initIndices) #cant ignore mode mode=c("all", "out", "in") if edges are directed
  neighbors <- unique(unlist(neighbors))

  context <- neighbors[ ! neighbors %in% initIndices ]

  xn <- length(context)

  subgraph <- induced.subgraph(g, which(V(g) %in% initIndices))

  connectedSymbols <- V(subgraph)$label[degree(subgraph)>0]
  disconnectedSymbols <- V(subgraph)$label[degree(subgraph)==0]
  uncoveredSymbols <- pathwaySymbols[ ! pathwaySymbols %in% V(g)$label ]

  n <- length(pathwaySymbols)
  cn <- length(connectedSymbols)
  dn <- length(disconnectedSymbols)
  un <- length(uncoveredSymbols)

  df <- data.frame(total = n,
                   connected = cn,
                   disconnected = dn,
                   uncovered = un,
                   connectedPerc = round(cn/n,2),
                   disconnectedPerc = round(dn/n,2),
                   uncoveredPerc = round(un/n,2),
                   context = xn)

  returnThis <-    c(list(connectedSymbols),
                     list(disconnectedSymbols),
                     list(uncoveredSymbols),
                     list(df))
  return (returnThis)
}

#' @export
GetQuackModelingFile <- function(genesetHash,inSamplePathways, logOn=TRUE,stepSize=50,normMetrics){

  # genesetHash = useThisPathwayHash;inSamplePathways = samplePathways;logOn = TRUE;stepSize = 50; normMetrics = TRUE

  # inSamplePathways: pathway, takeThisContext
  PathwayList <- unique(inSamplePathways$pathway)
  upperBound <- length(PathwayList)

  modelingFile          <- data.frame()
  otherFile             <- NULL
  connnectivityInfoFile <- data.frame()
  firstPathway          <- TRUE

  pb <- txtProgressBar(max = upperBound, style = 3)
  progress <- function(n) { setTxtProgressBar(pb, n) }
  opts <- list(progress = progress)

  r <- foreach(i = 1:upperBound,
#               .combine  = 'rbind',
               .export=ls(envir=globalenv()),
               .packages=c("hash","igraph","pROC"),
               .options.snow = opts) %dopar% {

               thisPathwayRow <- inSamplePathways[i,]
               pathway <- as.character(thisPathwayRow$pathway)
               contextSampleSize <- as.numeric(thisPathwayRow$takeThisContext)

               pathwayGeneSymbols <- genesetHash[[pathway]]

               pathwayConnectivityInfo <- getConnectivityInfo(pathwayGeneSymbols)
               pathwayConnectivityInfo[[4]]$pathway <- pathway

               #get the scored pathway and context genes
               thisPathwayData <- computeStatsForAllConnectedNodes(pathwayGeneSymbols,
                                                                   usePathwaySize = FALSE,
                                                                   multFactor = 0,
                                                                   sampleSize=contextSampleSize,
                                                                   normalizeMetrics = normMetrics)

               if(nrow(thisPathwayData)>0){
                 thisPathwayData$pathway <- pathway
               }else{
                 thisPathwayData <- data.frame()
               }

               #create obs for the unconnected and unconvered nodes
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
                 unconnectedDF$pathway <- pathway
                 uncoveredDF$pathway <- pathway
                 otherObs <- rbind(unconnectedDF,uncoveredDF)
               }else if ( rowsUnc > 0 && ! (rowsUnv > 0) ){
                 unconnectedDF$pathway <- pathway
                 otherObs <- unconnectedDF
               }else if (! (rowsUnc > 0) && rowsUnv > 0 ){
                 uncoveredDF$pathway <- pathway
                 otherObs <- uncoveredDF
               }
               returnList <- c(list(thisPathwayData),list(otherObs),list(pathwayConnectivityInfo[[4]]))
  }

  #now append results to data frame
  for (i in 1:length(r)){
    if(firstPathway == TRUE){
      modelingFile <- data.frame(r[[i]][[1]])
      if(nrow(data.frame(r[[i]][[2]]))>0) otherFile <- rbind(otherFile, data.frame(r[[i]][[2]]))
      connnectivityInfoFile <- data.frame(r[[i]][[3]])
      firstPathway <- FALSE
    }else{
      modelingFile <- rbind(modelingFile, r[[i]][[1]])
      if(nrow(data.frame(r[[i]][[2]]))>0) otherFile <- rbind(otherFile, data.frame(r[[i]][[2]]))
      connnectivityInfoFile <- rbind(connnectivityInfoFile, r[[i]][[3]])
    }
  }

#  modelingFile          <- data.frame(r[[1]])
#  otherFile             <- data.frame(r[[2]])
#  connnectivityInfoFile <- data.frame(r[[3]])

  modelingFile$degreeRatio  <- ifelse(modelingFile$giDegreeInG==0,0,modelingFile$giDegree / modelingFile$giDegreeInG)
  modelingFile$wdegreeRatio <- ifelse(modelingFile$giWDegreeInG==0,0,modelingFile$giWDegree / modelingFile$giWDegreeInG)

  modelingFile$ecRatio   <- ifelse(modelingFile$giECInG==0, 0,modelingFile$giEigenv / modelingFile$giECInG)
  modelingFile$btwnRatio <- ifelse(modelingFile$giBetweenessInG==0, 0,modelingFile$giBetweenness / modelingFile$giBetweenessInG)
  modelingFile$closRatio <- ifelse(modelingFile$giClosenessInG==0, 0,modelingFile$giCloseness / modelingFile$giClosenessInG)
  modelingFile$lccRatio  <- ifelse(modelingFile$giLCCInG==0, 0,modelingFile$giLCC / modelingFile$giLCCInG)

  modelingFile$ecRatio <- ifelse(modelingFile$ecRatio==0,0,log(modelingFile$ecRatio))
  modelingFile$btwnRatio <- ifelse(modelingFile$btwnRatio==0,0,log(modelingFile$btwnRatio))
  modelingFile$closRatio <- ifelse(modelingFile$closRatio==0,0,log(modelingFile$closRatio))
  modelingFile$lccRatio <- ifelse(modelingFile$lccRatio==0,0,log(modelingFile$lccRatio))

  #omit missing values and change the type of the response for modeling
  modelingFile$isInPathway <- suppressWarnings(as.factor(modelingFile$isInPathway))
  modelingFile <- na.omit(modelingFile)
  modelingFile$isInPathway <- modelingFile$isInPathway[,drop=TRUE]
  modelingFile$isInPathway <- suppressWarnings(as.factor(modelingFile$isInPathway))

  #add predictor columns
  addTheseColumns <- colnames(modelingFile)

  rowsInOther <- nrow(otherFile)
  if( is.null(rowsInOther) || is.na(rowsInOther) ) rowsInOther <- 0

  if(rowsInOther > 0){

    for(col in addTheseColumns){
      if(col!="isInPathway" && col!="gi" && col != "pathway"){
        otherFile[col] <- rep(0,nrow(otherFile))
      }
    }

    #########################################
    #remove uncovered and unconnected genes##
    #########################################
    modelingFile <- rbind(modelingFile,otherFile) # combine connected and rest
    modelingFile <- modelingFile[sample(nrow(modelingFile)),] #shuffle

  }

  return (list(modelingFile, connnectivityInfoFile))
}

#' @export
getGeneIndices <- function(geneNames){
  initIndices <- which(V(g)$label %in% geneNames)
  return (initIndices)
}

#' @export
getConnectedGeneIndices <- function(geneNames){
  initIndices <- which(V(g)$label %in% geneNames)
  subgraph <- induced.subgraph(g, which(V(g) %in% initIndices))
  returnIndices <- V(subgraph)$indexInG[degree(subgraph)>0]
  return (returnIndices)
}

#' @export
getDisconnectedGeneIndices <- function(geneNames){
  initIndices <- which(V(g)$label %in% geneNames)
  subgraph <- induced.subgraph(g, which(V(g) %in% initIndices))
  returnIndices <- V(subgraph)$indexInG[degree(subgraph)==0]
  return (returnIndices)
}

#' @export
computeStatsForAllConnectedNodes <- function(inSymbols,usePathwaySize,multFactor,sampleSize, normalizeMetrics){

  # inSymbols=pathwayGeneSymbols;usePathwaySize = FALSE;multFactor = 0;sampleSize=contextSampleSize;normalizeMetrics = normMetrics

  allIndices <- getGeneIndices(inSymbols)
  connectedIndices <- getConnectedGeneIndices(inSymbols)
  disconnectedIndices <- getDisconnectedGeneIndices(inSymbols)

  neighbors <- neighborhood(g, 1, nodes=connectedIndices) #cant ignore mode mode=c("all", "out", "in") if edges are directed
  neighbors <- unique(unlist(neighbors))

  #context <- neighbors[ ! neighbors %in% allIndices ]
  contextSub <- c()
  #contextSize <- length(context)

  combinedIndicies <- union(allIndices,neighbors)
  subgraphContext <- induced.subgraph(g, which(V(g) %in% combinedIndicies))
  giDegreeOneSymbols <- V(subgraphContext)$label[degree(subgraphContext)>1] # context and pathwya >1
  gn <- getGeneIndices(giDegreeOneSymbols) #get indices
  context <- gn[!gn%in%allIndices] # remove pathway indices to get only context (added !)
  contextSize <- length(context) # reset context size to only those with degree>1

  if(usePathwaySize){
    contextSub <- sample(context, length(allIndices)*multFactor)
  }else{
    if ( contextSize > sampleSize ){
      contextSub <- sample(context, sampleSize)
    }else{
      contextSub <- context
    }
  }

  #connected indices and those that are in context (but not unconnected pathway nodes)
  finalList <- c(contextSub, connectedIndices)

  #remove pathway genes with giDegree==1
  #finalList <- finalList[-which(finalList%in%gn)]

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

    if(normalizeMetrics){
      min.degree <- min(tempdegree); max.degree <- max(tempdegree)
      if(min.degree==max.degree){
        tempdegree <- rep(1, length(tempdegree))
      }else{
        tempdegree <- (tempdegree - min.degree) /  (max.degree - min.degree)
      }
    }

    giDegree <- tempdegree[indexOfgi]

    #tempec <<- NULL
    tryCatch(suppressWarnings(tempec <- evcent(temppg)$vector),
             error = function(e){ tempec <- rep(0, vcount(temppg))    })

    #tempec <- evcent(temppg)$vector
    giEigenv <- tempec[indexOfgi]

    if(giDegree==0) giEigenv <- 0

    tempwdegree <- graph.strength(temppg,loops = FALSE)

    if(normalizeMetrics){
      min.degree <- min(tempwdegree); max.degree <- max(tempwdegree)
      if(min.degree==max.degree){
        tempwdegree <- rep(1, length(tempwdegree))
      }else{
        tempwdegree <- (tempwdegree - min.degree) /  (max.degree - min.degree)
      }
    }

    giWDegree <- tempwdegree[indexOfgi]

    closeness <- closeness(temppg, normalized = normalizeMetrics) # change to normalized = TRUE
    giCloseness <- closeness[indexOfgi]

    betweenness <- betweenness(temppg,normalized = normalizeMetrics)  # change to normalized = TRUE
    giBetweenness <- betweenness[indexOfgi]

    lcc <- transitivity(temppg, type="local") # local clustering coefficient
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

#' @export
splitTrainingAndHoldout <- function(file, holdoutPerc, byPathway){

  if(byPathway){
    unique.pathways <- unique(file$pathway)
    training.pathways <- sample(unique.pathways, (1-holdoutPerc)*length(unique.pathways))
    file$holdout <- ifelse(file$pathway %in% training.pathways,1,0)
    file.holdout <- subset(file,holdout==1)
    file.training <- subset(file,holdout==0)
  }else{
    file$holdout <- ifelse(runif(nrow(file))>holdoutPerc,1,0)
    file.holdout <- subset(file,holdout==1)
    file.training <- subset(file,holdout==0)
  }

  return (c(list(file.holdout), list(file.training)))
}

#' @export
getVariableImportance <- function(model,typeparm){

  ave_oob_error <- importance(model,type=typeparm)
  vars <- ave_oob_error[,0]
  ave_oob_error <- data.frame(ave_oob_error)
  ave_oob_error$vars <- vars
  names(ave_oob_error)[1] <- paste("score")
  order.ave_oob_error <- ave_oob_error[order(-ave_oob_error$score),]

  return (order.ave_oob_error)
}

#' @export
saveVariableImportance <- function(location, name,predList, indata){

  sink(paste(location, name,".txt",sep=""))
  cat(name)
  cat("\n")
  for(i in 1:length(predList)){
    cat(paste(predList[i], " ",round(indata$score[i],2),sep=""))
    cat("\n")
  }
  sink()
}

#' @export
GetQuackClassifier <- function(trainingData, numberOfTrees=500, withReplacement,fileName, saveResults, varImportDir,geneModDir, num_cores){
#  QuackGeneModel <- randomForest(y=response,
#                                 x=predictors,
#                                 importance = TRUE,
#                                 replace = withReplacement,
#                                 na.action = na.omit,
#                                 ntree = numberOfTrees)


#  order.ave_oob_error <- getVariableImportance(QuackGeneModel,1)
#  predList <- row.names(order.ave_oob_error)


#  return (c(list(QuackGeneModel), list(order.ave_oob_error)))

  QuackGeneModel <- ranger(formula      = isInPathway ~ .,
                           data         = trainingData,
                           importance   = 'permutation',
                           replace      = withReplacement,
                           num.trees    = numberOfTrees,
                           num.threads  = num_cores,
                           write.forest = T,
                           probability  = T)

  order.ave_oob_error <- QuackGeneModel$variable.importance
#  predList <- row.names(order.ave_oob_error)

#  if(saveResults){
#    saveVariableImportance(varImportDir, paste("VarImportance_", fileName,sep=""),predList, order.ave_oob_error)
#    save(QuackGeneModel,file = paste(geneModDir, "QuackGeneModel_",fileName,".RData",sep=""),ascii=TRUE)
#  }

  return (c(list(QuackGeneModel), list(order.ave_oob_error)))
}

#' @export
#bootstrapAUC_tm <- function(aucFile, stepSize, bootSampleSize,samplePerc){
bootstrapAUC_tm <- function(aucFile, bootSampleSize, samplePerc, usr_output_loc){

#  aucFile <- holdoutsForAUCs[[1]]
#  bootSampleSize <- 250
#  samplePerc <- 0.1

  numObs <- nrow(aucFile)
  indexList <- c(1:numObs)

#    r <- foreach(i = start_i:end_i,
#                 .export=ls(envir=globalenv()),
#                 .packages=c("randomForest","igraph","pROC")) %dopar% {

#  print("--Start bootstrapping")

  pb <- txtProgressBar(max = bootSampleSize, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  r <- foreach(i = 1:bootSampleSize,
#               .combine = 'rbind',
               .export=ls(envir=globalenv()),
               .packages=c("pROC"),
               .options.snow = opts) %dopar% {
                   this.values <- sample(indexList, round(samplePerc*numObs), replace = TRUE)
                   this.DS     <- aucFile[this.values,]
                   #this.AUC <- NA
                   #try(thisAUC <- as.numeric(auc(thisDS$isInPathway, thisDS$QuackP)),silent = TRUE)
                   this.AUC    <- roc(this.DS$isInPathway, this.DS$QuackP)
                }
#  print("--Done bootstrapping")


  #now append results to data frame
  finalAUCs <- c()
  sensitivityList <- list()
  specificityList <- list()
  for (i in 1:length(r)){
    finalAUCs <- c(finalAUCs, auc(r[[i]]))
    sensitivityList[[length(sensitivityList)+1]] <- r[[i]]$sensitivities
    specificityList[[length(specificityList)+1]] <- r[[i]]$specificities
  }


  remove(r)
  cleanMem()

  finalAUCs <- finalAUCs[! is.na(finalAUCs)]
  sensitivityList <- sensitivityList[! is.na(finalAUCs)]
  specificityList <- specificityList[! is.na(finalAUCs)]
  numObs <- length(finalAUCs)
  if(numObs > 30){
    quackAUCCI <- quantile(finalAUCs, c(0.025,0.5,0.975))
  }else{
    quackAUCCI <- c(NA,NA,NA)
    names(quackAUCCI)<-c("2.5%","50%","97.5%")
  }
  quackAUCCI <- quantile(finalAUCs, c(0.025,0.5,0.975))


  pathways <- unique(aucFile$pathway)
  pb <- txtProgressBar(max = length(pathways), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  r <- foreach(i = 1:length(pathways),
               .combine = 'rbind',
               .export=ls(envir=globalenv()),
               .packages=c("pROC"),
               .options.snow = opts) %dopar% {
                 this.DS  <- aucFile[aucFile$pathway == pathways[i], c("isInPathway", "QuackP")]
                 this.AUC <- auc(roc(this.DS$isInPathway, this.DS$QuackP))[[1]]
                 df <- data.frame(pathway=pathways[i], auc=this.AUC)
              }
  write.table(r, paste(usr_output_loc, modelSpecs, '_pathways.tsv', sep=""), quote=F, col.names=T, row.names=F)

#  specificityOfInterest <- specificityList[[which.quantile(finalAUCs,0.5)]]
#  sensitivityOfInterest <- sensitivityList[[which.quantile(finalAUCs,0.5)]]
#  remove(specificityList)
#  remove(sensitivityList)
#  save(specificityOfInterest,sensitivityOfInterest,finalAUCs,
#       file = paste(SaveQuackGeneModelHere, "QuackAUCROC_",networkKey,".RData",sep=""),
#       ascii=TRUE)

  return(quackAUCCI)
}

#' @export
saveAUCCI <- function(ci, location,name){
  sink(paste(location, name,sep=""))
  cat(paste( "AUC: ",round(ci[[2]],2), " C.I. (", round(ci[[1]],2),",", round(ci[[3]],2),")",sep=""))
  sink()
}

#' @export
splitHoldoutForAUC <- function(holdoutFile, model, response, predList){
#  holdoutFile <- modelingFile.holdout
#  model <- QuackGeneModel
#  response <- "isInPathway"
#  predList <- predictorList

#  tempP <- predict(model,newdata = holdoutFile[,c(predList,response)],type = "prob")
  tempP <- predict(model, data = holdoutFile[,c(predList,response)])

#  probIsPathway <- tempP[,2]
  probIsPathway <- tempP$predictions[,2]

  holdoutFile$QuackP <- probIsPathway

  #correct
  holdoutFile$QuackP <- ifelse(holdoutFile$gi < 1,0,holdoutFile$QuackP) #adjust unconnected

  con <- subset(holdoutFile, gi>0)
  conAndUnc <- subset(holdoutFile, gi>=0)
  conAndUncAndUncv <- subset(holdoutFile, gi>=-1)

  # connected only
  con <- con[,c("isInPathway", "QuackP", "pathway")]
  con <- con[c('QuackP','isInPathway', 'pathway')]
  con <- con[order(con$QuackP,decreasing=TRUE),] #sort

  conAndUnc <- conAndUnc[,c("isInPathway", "QuackP", "pathway")]
  conAndUnc <- conAndUnc[c('QuackP','isInPathway', 'pathway')]
  conAndUnc <- conAndUnc[order(conAndUnc$QuackP,decreasing=TRUE),] #sort

  conAndUncAndUncv <- conAndUncAndUncv[,c("isInPathway", "QuackP", "pathway")]
  conAndUncAndUncv <- conAndUncAndUncv[c('QuackP','isInPathway', 'pathway')]
  conAndUncAndUncv <- conAndUncAndUncv[order(conAndUncAndUncv$QuackP,decreasing=TRUE),] #sort

  return (c(list(con), list(conAndUnc), list(conAndUncAndUncv)))
}

#' @export
saveAUCVariantsToCSV <- function(con,andunc,anduncv, location,name){

  fileName <- paste(location, name,sep="")

  df <- data.frame()
  df <- rbind(df, data.frame(method="con", estimate=round(con[[2]],2), ci= paste("(", round(con[[1]],2),",", round(con[[3]],2),")",sep="")))
  df <- rbind(df, data.frame(method="andunc", estimate=round(andunc[[2]],2), ci= paste("(", round(andunc[[1]],2),",", round(andunc[[3]],2),")",sep="")))
  df <- rbind(df, data.frame(method="anduncv", estimate=round(anduncv[[2]],2), ci= paste("(", round(anduncv[[1]],2),",", round(anduncv[[3]],2),")",sep="")))

  colnames(df) <- c("estimate", "ci")

  write.csv(df, file=fileName)
}

#' @export
returnMinSec <- function(x) {
  result <- paste(floor(x/60)," Min and ", round(x%%60)," Sec",sep="")
  return (result)
}

#' @export
cleanMem <- function(n=10) { for (i in 1:n) gc() }
