#' Train model function
#' 
#' Using user uploaded network and known pathways, train model using Quack algorithm.
#' @param usr_modeling_file,usr_pathway_file,usr_output_loc Characters referring to the location and file name of user's network modeling file, path and name of pathway files, which can be found in data/ or user uploaded, and location to save output.
#' @return rdata file, csv and txt file
#' @import hash igraph randomForest pROC doMC dplyr boot stringr scales lattice MASS mclust changepoint data.table entropy poweRlaw doParallel
#' @export
#' @examples 
#' trainModel("data/InWeb3.RData", "data/stringent835Pathways.RData", "~")
trainModel <- function(trained_network_file, pathway_file, save_location){
  
  ### Register multicore parallel backend 
  registerDoMC(4)
  ### Register doParallel parallel backend
  cl <- makeCluster(4)
  registerDoParallel(cl)
  
  load(file=(usr_modeling_file))
  load(usr_pathway_file)
  
  pathwayName <- strsplit(basename(usr_pathway_file), "\\.")[[1]][1]
  pathwayClassKey <- "General"
  
  quackVersion <- "Quackv1.3"
  predictorList <- c("degreeRatio","wdegreeRatio","ecRatio", "btwnRatio", "closRatio", "lccRatio","giEigenv","giDegree","giWDegree","giCloseness", 
                     "giBetweenness","giLCC","giDegreeInG","giWDegreeInG","giECInG","giClosenessInG","giBetweenessInG","giLCCInG")

  networkKeys <- keys(usrNetworkGraphs)
  useThisPathwayHash <- canonicalPathwayHash
  
  maxFileSize <- 6e5
  numberOfTrees <- 500
  monitor <- TRUE
  NormalizeMetrics <- FALSE
  DoRFSamplingWithReplacement <- FALSE
  SplitTrainingAndHoldout.ByPathway <- TRUE
  
  for(networkKey in networkKeys){
  
    g <<- usrNetworkGraphs[[networkKey]]
    ### Export 'g' to each node
    clusterExport(cl, "g")
    
    pl <- keys(canonicalPathwayHash)
    
    modelSpecs <<- paste(quackVersion, "-",networkKey,"-",pathwayClassKey,pathwayName,sep="")
    
    if(monitor) print(paste("MODEL: ", modelSpecs, sep=""))
    
    Sample.StartTime <- proc.time()
    samplePathways <- getCustomSampleForTrainingPathways(useThisPathwayHash,
                                                         pl,
                                                         restrictionType = "restrictTo",
                                                         sampleSize = maxFileSize,
                                                         doDownSample = TRUE, 
                                                         downSampleSize = 1000,
                                                         pathwayBounds = c(10,100),
                                                         sortByConnected =  FALSE)
    Sample.RunTime <- returnMinSec((proc.time() - Sample.StartTime)[3])
    
    ModFile.StartTime <- proc.time()
    modelingAndConnInfoFile <- GetQuackModelingFile(genesetHash = useThisPathwayHash,
                                                    inSamplePathways = samplePathways,
                                                    logOn = TRUE,
                                                    stepSize = 50,
                                                    normMetrics = NormalizeMetrics)
    ModFile.RunTime <- returnMinSec((proc.time() - ModFile.StartTime)[3])
    
    if(monitor) print("--Completed modelingAndConnInfoFile")
    
    modelingFile <- modelingAndConnInfoFile[[1]]
    connectivityInfoFile <- modelingAndConnInfoFile[[2]]
    
    write.csv(connectivityInfoFile, file=paste(usr_output_loc, modelSpecs,".csv",sep=""))
    
    ### Training / holdouts
    modelingFile$modFileKey <- seq(1:nrow(modelingFile))
    splitFile <- splitTrainingAndHoldout(modelingFile, 0.7, byPathway = SplitTrainingAndHoldout.ByPathway)
    modelingFile.holdout <- splitFile[[1]]
    modelingFile.training <<- splitFile[[2]]
    holdout.keys <- modelingFile.holdout$modFileKey
    modelingFile$isHoldout <- ifelse(modelingFile$modFileKey %in% holdout.keys, 1, 0)
    
    save(modelingFile.holdout,file = paste(usr_output_loc, "Holdout_",modelSpecs,".RData",sep=""),ascii=TRUE)
    
    ModelResults <- GetQuackClassifier(response      = modelingFile.training[,c("isInPathway")],
                                       predictors    = modelingFile.training[,predictorList], 
                                       numberOfTrees = numberOfTrees,                 ### Changed from 500 trees 
                                       withReplacement = DoRFSamplingWithReplacement, ### Do sampling with or without replacement?
                                       fileName      = modelSpecs, 
                                       saveResults   = TRUE,
                                       #varImportDir  = SaveVariableImportanceHere,
                                       #geneModDir    = SaveQuackGeneModelHere,
                                       varImportDir  = usr_output_loc,
                                       geneModDir    = usr_output_loc)
    
    QuackGeneModel <- ModelResults[[1]]
    remove(ModelResults)
    
    save(QuackGeneModel, paste0(usr_output_loc, modelSpecs, '.RData', sep=""), ascii = T)
    #save(QuackGeneModel, paste(usr_output_loc, modelSpecs, '.RData', sep=""), ascii = T)
    save(QuackGeneModel, file="R/sysdata.rda")
    
    if(monitor) print("--Completed QuackGeneModel")
    
    remove(modelingFile.training)
    
    #####################
    #       AUCs        #
    #####################
    
    holdoutsForAUCs <- splitHoldoutForAUC(modelingFile.holdout, QuackGeneModel, "isInPathway", predictorList)
    
    aucCI <- bootstrapAUC_tm(aucFile <- holdoutsForAUCs[[1]],
                          stepSize <- 10,
                          bootSampleSize <- 250,
                          samplePerc <- 0.7)
    
    saveAUCCI(aucCI, usr_output_loc, paste(modelSpecs,"-previous.txt",sep=""))
    
    remove(QuackGeneModel)
    
    auc.ci.c <- ci.auc(holdoutsForAUCs[[1]][,2],holdoutsForAUCs[[1]][,1],levels = c(0,1), direction="<")
    auc.ci.cu <- ci.auc(holdoutsForAUCs[[2]][,2],holdoutsForAUCs[[2]][,1],levels = c(0,1), direction="<")
    auc.ci.cuv <- ci.auc(holdoutsForAUCs[[3]][,2],holdoutsForAUCs[[3]][,1],levels = c(0,1), direction="<")
    
    saveAUCVariantsToCSV(auc.ci.c,auc.ci.cu,auc.ci.cuv, usr_output_loc, paste("AUC-CIs-",modelSpecs,".csv",sep=""))
    
    save(modelingFile,file = paste(usr_output_loc, "ModelingFile_", modelSpecs, ".RData", sep=""), ascii=TRUE)
    
    remove(modelingFile)
    
    sink(paste(usr_output_loc,modelSpecs, "-DONE",sep=""))
    cat(":-)")
    sink()
    
    remove(g)
  }
}
