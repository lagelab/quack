#' Make your graphs function
#' 
#' This function allows user to upload network file and create graph function compatible with Quack.
#' @param usr_network_file,usr_network_name,usr_save_loc Characters referring to the location and file name of user's network file, name of network, and location to save output.
#' @return rdata file
#' @import hash igraph 
#' @export
#' @examples 
#' makeGraphc("data/InWeb3_network.txt", "InWeb3", "data/")
makeGraphs <- function(usr_network_file, usr_network_name, usr_save_loc) {
  
  file.name <- usr_network_file
  network.name <- usr_network_name
  con <- file(file.name, 'r') 
  InWebFile <- readLines(con)
  close(con)
  
  getEdgeKey <- function(n1,n2){
    key <- ""
    if (n1 > n2){
      key <- paste(n1, "!", n2, sep="")
    } else{
      key <- paste(n2, "!", n1, sep="")
    }
    return (key)
  }
  
  InWebHash <- hash()
  i<-1
  for(line in InWebFile){
    theSplit <- unlist(strsplit(line, "\t"))
    
    src <- gsub("(^ +)|( +$)", "",theSplit[1])
    trg <- gsub("(^ +)|( +$)", "",theSplit[2])
    score <- as.numeric(gsub("(^ +)|( +$)", "",theSplit[3]))
    if(is.na(score)){
      print(paste("Source: ", src, "  Target: ",trg, "  Score: ", score, sep=""))
    }
    if(src != "gene1" && (src!=trg)){
      key <- getEdgeKey(src,trg)
      .set(InWebHash, key, score) 
    }  
    
    if(i<10){
      print(paste("Source: ", src, "  Target: ",trg, "  Score: ", score, sep=""))
    }
    i<-i+1
  }
  
  splitEdgeKey <- function(key){
    theSplit <- unlist(strsplit(key, "!"))
    return (theSplit)
  }
  
  makeGraphFromScoreHash <- function(inHash){
    allEdges <- keys(inHash)
    
    nodeHash <- hash()
    for(e in allEdges){
      theSplit <- splitEdgeKey(e)
      src <- theSplit[1]
      trg <- theSplit[2]
      .set( nodeHash, src, TRUE )
      .set( nodeHash, trg, TRUE )  
    }
    
    print("Completed nodeHash")
    allGenes <- keys(nodeHash)
    
    geneToIndexHash <- hash()
    indexToGeneHash <- hash()
    index <- 1
    for(gene in allGenes){
      .set(geneToIndexHash, gene,index)
      .set(indexToGeneHash, as.character(index),gene)
      index <- index + 1
    }
    
    numMatrixRows <- length(allEdges)
    scoresGraph <- matrix(nrow = numMatrixRows, ncol = 1, byrow = TRUE)
    edgeListGraph <- matrix(nrow = numMatrixRows, ncol = 2, byrow = TRUE)
    
    index <- 1
    for (e in allEdges){
      nodes <- splitEdgeKey(e)
      n1 <- nodes[1]
      n2 <- nodes[2]
      edgeListGraph[index,1] <- geneToIndexHash[[n1]]
      edgeListGraph[index,2] <- geneToIndexHash[[n2]]
      scoresGraph[index,1] <- inHash[[e]]
      index <- index + 1 
    }
    print("Completed scoresGraph")
    
    ### Initialize the graph
    g <- graph.edgelist( edgeListGraph , directed=FALSE)
    
    V(g)$label <- allGenes
    V(g)$indexInG <- seq(1,length(allGenes))
    
    ### Assign weights
    E(g)$weight <- scoresGraph
    
    print("starting global properties for nodes")
    
    V(g)$degree <-  degree(g,loops = FALSE)
    V(g)$wDegree <- graph.strength(g,loops = FALSE)
    V(g)$ec <- evcent(g)$vector
    V(g)$closeness <- closeness(g, normalized = FALSE)
    V(g)$betweenness <- betweenness(g, normalized = FALSE)
    V(g)$lcc <- transitivity(g, type="local")
    V(g)$lcc[is.na(V(g)$lcc)] <- 0  #correct NA's with 0
    
    return (g)
  }
  
  Network <- makeGraphFromScoreHash(InWebHash)
  
  usrNetworkGraphs <- hash()
  .set(usrNetworkGraphs,network.name, Network)
  save(usrNetworkGraphs, file = paste(usr_save_loc, network.name, ".RData",sep=""))
}
