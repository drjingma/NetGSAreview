library(ppcor)
library(gdata)
library(graphite)
library(igraph)
library(graph)

# May want to write a wrapper function for each method

#This function drops genes based on whether gene.id is TRUE or FALSE
drop.genes = function(dataset, gene.id = NULL){
  if(is.null(gene.id)){print("no gene selected")}
  else{
    dataset[["gene_info"]] = dataset[["gene_info"]][!gene.id,]
    rownames(dataset[["gene_info"]]) <- NULL
    dataset[["dat"]] = dataset[["dat"]][!gene.id,]
  }
  return(dataset)
}

drop.metabs = function(dataset, metab.id = NULL){
  if(is.null(metab.id)){print("no metabolite selected")}
  else{
    dataset[["metab_info"]] = dataset[["metab_info"]][!metab.id,]
    rownames(dataset[["metab_info"]]) <- NULL
    dataset[["dat"]] = dataset[["dat"]][!metab.id,]
  }
  return(dataset)
}


drop.samples = function(dataset, sample.id = NULL){
  if(is.null(sample.id)){print("no sample selected")}
  else{
    dataset[["sample_info"]] = (dataset[["sample_info"]][!sample.id,])
    dataset[["dat"]] = (dataset[["dat"]][,!sample.id])
  }
  return(dataset)
}

#' My function for comparing nominal pvalues, fdr and family-wise error rate adjusted pvalues
#' with the significance level alpha
adjust.pvalues <- function(pvals,alpha=0.05){
  apply(data.frame('nominal'=pvals, 'fdr'=p.adjust(pvals, 'BH'), 'fwer'=p.adjust(pvals, 'holm')), 2, function(a) (a<alpha))
}

## Community pathway deregulation
## For a given pathway, find communities and choose a subset of communities such that
## genes in these communities are affected with the desired mean changes. 
community.dereg <- function(pathway, genes, DC=0.5, maxIter=10, seed=1){
  # pathway: a pathway object
  # genes: the genes that are used in the enrichment analysis. This can be learned from the B matrix.
  # DC: detection call, a percentage between 0% and 100% indicating the proportiontion of genes to be affected
  #  In this design, DC is fixed to be about 50% to make it easier for control.
  set.seed(seed)
  size <- length(genes)
  
  el <- graph::edges(pathway)
  el$src <- paste0("ENTREZID:", el$src)
  el$dest <- paste0("ENTREZID:", el$dest)
  
  gr <- igraph::graph_from_edgelist(cbind(el$src, el$dest), directed = F)
  
  # get ride of duplicated edges
  adj <- as.matrix(igraph::get.adjacency(gr, type="both"))
  adj[(adj>0)] = 1;   diag(adj) <- 0;
  gr <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
  #subgraph based on the variables in genes
  gr <- igraph::induced_subgraph(gr, genes)
  
  lec <- igraph::cluster_edge_betweenness(gr)
  tmp <- as.numeric(table(lec$membership))
  proportion <- tmp/size
  names(proportion) <- names(table(lec$membership))
  proportion <- sort(proportion, decreasing = T)

  if (proportion[1] > DC){
    # if the first community has more than half the nodes, then sample from this community
    mem <- names(proportion)[1]
    subset <- lec$names[which(lec$membership==mem)]
    subset <- sample(subset, floor(size*DC))
  } else{
    ## need to choose more than one community
    # We first check if there exists a subset of communities that meet the DC criterion
    cumsum_prop <- cumsum(tmp)/size
    tmp_order <- 1:max(lec$membership)
    cnt <- 0
    ## can use a random permutation of tmp to find the cumsum proportion
    while(min(abs(cumsum_prop-DC))>0.05 && cnt < maxIter){
      tmp_order <- sample(1:max(lec$membership))
      tmp <- tmp[tmp_order]
      cumsum_prop <- cumsum(tmp)/size
      cnt <- cnt + 1
    }
    id <- which.min(abs(cumsum_prop-DC))
    names_order <- names(table(lec$membership))[tmp_order]
    mem <- names_order[1:id]
    subset <- lec$names[(lec$membership %in% mem)]
    
    ## check realDC
    if (abs(length(subset)/size-DC)<0.05){
      subset = subset
    } else{
      ## The above procedure doesn't produce a good subset.
      ## Hence we choose a random sample
      cumsum_prop <- cumsum(proportion)
      id <- which.min(cumsum_prop<DC)
      mem <- names(cumsum_prop)[1:id]
      subset <- lec$names[(lec$membership %in% mem)]
      subset <- sample(subset, floor(size*DC))
    }
  }
  res <- list( GeneID=subset, realDC = length(subset)/size, lec=lec$membership) 
  
  return(res)
}
  

betweenness.dereg <- function(pathway, genes, DC=0.5, seed=1){
  set.seed(seed)
  size <- length(genes)
  el <- graph::edges(pathway)
  el$src <- paste0("ENTREZID:", el$src)
  el$dest <- paste0("ENTREZID:", el$dest)
  
  gr <- igraph::graph_from_edgelist(cbind(el$src, el$dest), directed = F)
  
  # get ride of duplicated edges
  adj <- as.matrix(igraph::get.adjacency(gr, type="both"))
  adj[(adj>0)] = 1;   diag(adj) <- 0;
  gr <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
  #subgraph based on the variables in genes
  gr <- igraph::induced_subgraph(gr, genes)
  
  ## calculate edge betweenness based on gr
  vb <- igraph::betweenness(gr, directed = F)

  threshold <- quantile(vb, 1-DC)
  subset <- names(vb)[which(vb>=threshold)] # this subset consists of at least the desired DC, and could be larger
  
  if (abs(length(subset)/size-DC)<0.05){
    subset = subset
  } else{
    ## The above procedure doesn't produce a good subset.
    ## Hence we choose a random sample
    subset <- sample(subset, floor(size*DC))
  }
  return(list(GeneID=subset, realDC=length(subset)/size))
}


## nei: neighborhood distance
## Choose the node that has the largest number of neighbors, together with its neighborhood
neighborhood.dereg <- function(pathway, genes, DC=0.5, nei=2, seed=1){
  set.seed(seed)
  size <- length(genes)
  
  el <- graph::edges(pathway)
  el$src <- paste0("ENTREZID:", el$src)
  el$dest <- paste0("ENTREZID:", el$dest)
  gr <- igraph::graph_from_edgelist(cbind(el$src, el$dest), directed = F)
  
  # get ride of duplicated edges
  adj <- as.matrix(igraph::get.adjacency(gr, type="both"))
  adj[(adj>0)] = 1;   diag(adj) <- 0;
  gr <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
  #subgraph based on the variables in genes
  gr <- igraph::induced_subgraph(gr, genes)
  
  ## calculate edge betweenness based on gr
  ngbh <- igraph::ego(gr, nei)
  ngbh_len <- sapply(ngbh, length)
  names(ngbh_len) <- V(gr)$name
  id <- which.max(ngbh_len)
  subset <- ngbh[[id]]$name
  if (abs(length(subset)/size-DC)<0.05){
    subset = subset
  } else if (length(subset)/size-DC>0.05){
    ## The above procedure doesn't produce a good subset.
    ## Hence we choose a random sample
    subset <- sample(subset, floor(size*DC))
  }
  return(list(GeneID=subset, realDC=length(subset)/size))
}

DAGOnly<-function(pathways){
  pathways<-Filter(function(p) igraph::is_dag(p),pathways)
  return(pathways)
}

FewEdges<-function (pathways, minEdges) {
  pathways <- Filter(function(p) nrow(edges(p)) > minEdges , pathways)
  return(pathways)
}

BigPaths<-function (pathways, maxNodes) {
  pathways <- Filter(function(p) length(nodes(p)) <= maxNodes, pathways)
  return(pathways)
}

CommonGenes<-function (pathways, genes, threshold) {
  pathways<-Filter(function(p) length(intersect(nodes(p), genes)) >= threshold, pathways)
  return(pathways)
}

symAdj <- function(A){
  A <- A + t(A)
  A[A>0] <- 1
  A
}

#' @param pathways A list of pathways
#' @param method Which enrichment method
#' @param g A undirected graph with edges representing metabolic interactions
#' @param dataset A list that consists of metabolomic data, feature information and sample information
prepMetablicPathways<-function(pathways, method, g, metabolites){
  g <- igraph.from.graphNEL(g)
  g <- as.undirected(g, 'collapse')
  if (length(setdiff(metabolites$CID, V(g)$name))>0){
    # add vertices to graph
    g <- add_vertices(g, 2, name=setdiff(metabolites$CID,V(g)$name))
  }
  graList <- lapply(pathways, function(a)  induced_subgraph(g, intersect(V(g)$name, a)))
  
  N <- length(pathways)
  path_membership <- NULL
  adj <- NULL
  
  if (method=="camera"){
    pMembership <- lapply(pathways, function(a) 1*(match(metabolites$CID, a, nomatch = 0)>0))
    pMembership <- lapply(pMembership, function(a) {names(a)=metabolites$CID;a})
    path_membership <- lapply(pMembership, function(a) which(a==1))
    names(path_membership) <- names(pathways)
  }
  if (method=="DEGraph"){
    pathways <- lapply(graList, as_graphnel)
  }
  if (method=="netgsa"){
    # pathway indicator matrix
    pMembership <- lapply(pathways, function(a) 1*(match(metabolites$CID, a, nomatch = 0)>0))
    pMembership <- lapply(pMembership, function(a) {names(a)=metabolites$CID;a})
    path_membership <- do.call(rbind, pMembership)
    
    # network information
    g1 <- induced_subgraph(g, intersect(V(g)$name, metabolites$CID))
    adj <- as.matrix(get.adjacency(g1,type="both"))
    adj[(adj>0)] <- 1
    adj <- adj[match(metabolites$CID, rownames(adj), nomatch = 0),  match(metabolites$CID, rownames(adj), nomatch = 0)]
  }
  
  if (method=="PathNet"){
    # this is a data frame with c("src", "dest", "title")
    edgeList <- lapply(graList, as_edgelist)
    pathways4PathNet <- lapply(1:N, function(j)
      if (length(edgeList[[j]])==0){
        NULL
      } else {
        data.frame(edgeList[[j]], "title"=names(pathways)[j], stringsAsFactors = F)
      })
    pathways4PathNet <- do.call(rbind, pathways4PathNet)
    colnames(pathways4PathNet)[1:2] <- c("src","dest")

    # remember the id needs to be numeric not character
    pathways4PathNet$src <- as.numeric(gsub('C','',pathways4PathNet$src))
    pathways4PathNet$dest <- as.numeric(gsub('C','',pathways4PathNet$dest))
    pathways <- pathways4PathNet
    
    g1 <- induced_subgraph(g, intersect(V(g)$name, metabolites$CID))
    adj <- as.matrix(get.adjacency(g1, type="both"))
    adj[(adj>0)] <- 1
    adj <- adj[match(metabolites$CID, rownames(adj), nomatch = 0),  match(metabolites$CID, rownames(adj), nomatch = 0)]
    rownames(adj) <- as.numeric(gsub('C','',rownames(adj)))
    colnames(adj) <- as.numeric(gsub('C','',colnames(adj)))
  }
  
  if (method=="CePa"){
    # Here the path.catalogue only catelogues interactions, and pathways are listed with interaction id
    edgeList <- lapply(graList, as_edgelist)
    interactionList <- lapply(edgeList, function(a) {colnames(a)=c("input","output");a})
    interactionList <- lapply(1:length(interactionList), function(j) {
      if (length(interactionList[[j]])==0){
        NULL
      } else {
        data.frame('interaction.id'=paste(j, seq(1, nrow(edgeList[[j]])), sep='000'),interactionList[[j]], stringsAsFactors=FALSE)
      }
    })
    pathList <- lapply(interactionList, function(a) a$interaction.id)
    names(pathList) <- names(pathways)
    
    interactionList <- do.call(rbind, interactionList)
    mapping <- lapply(1:N, function(a) data.frame('node.id'=pathways[[a]], 'symbol'=metabolites$CID[match(pathways[[a]], metabolites$CID)], stringsAsFactors = FALSE))
    
    pathways <- set.pathway.catalogue(pathList = pathList,interactionList = interactionList, mapping = do.call(rbind, mapping))
    
  }

  return(list(pathways=pathways, membership=path_membership, adjacency=adj))
  
}
#' @param genes A data frame that consist of genes in the data, named by entrez id and symbol
#' @param pathways A pathwayList object obtained from graphite
#' @param maxNodes
#' @param minEdges
#' @param commonTh Number of shared genes between pathway and genes in the data
prepPathways<-function(pathways, method, genes, maxNodes=NULL, minEdges=NULL, commonTh=NULL, DAGonly=FALSE){
  
  if (!is.list(pathways) & !any(class(pathways)=="PathwayList")) stop("'pathways' must be a list or a PathwayList, found ", class(pathways))

  N<-length(pathways)
  
  if (!is.null(maxNodes)) pathways<-BigPaths(pathways, maxNodes)
  if (!is.null(commonTh)) pathways<-CommonGenes(pathways, genes$EntrezID, commonTh)
  if (!is.null(minEdges)) pathways<-FewEdges(pathways, minEdges)
  if (DAGonly && method %in% "topologyGSA") pathways<-DAGOnly(pathways)
  
  message(paste(N-length(pathways), " pathways were filtered out"))
  
  N <- length(pathways)
  path_membership <- NULL
  adj <- NULL
  
  if (method=="pe"){
    pathways <- lapply(pathways, pathwayGraph)
  }
  
  if (method=="netgsa"){
    ## individual pathway membership
    path_membership <- lapply(pathways, function(a) genes$symbol[which(match(genes$EntrezID, graph::nodes(a), nomatch = 0)>0)])
    names(path_membership) <- names(pathways)
    
    # file_path <- "../../../Michailidis_G/NetGSAextensions/BioGridNet/"
    # file_path <- "BioGrid/"
    el <- read.csv("BioGrid/BioGrid_human_March2019.csv")
    el$src <- paste0("SYMBOL:", el$src)
    el$dest <- paste0("SYMBOL:", el$dest)
    el_genenames <- unique(c(el$src, el$dest))
    g <- graph_from_edgelist(as.matrix(el), directed = FALSE)
    
    ## Remove genes that do not have matched measurements in data
    adj_single <- function(path_j){
      if (length(intersect(el_genenames, path_j))==0){
        warning(paste0('No shared nodes between the ', names(path_j), ' and network information'))
        netgsa_adj <- matrix(0, length(path_j), length(path_j))
        rownames(netgsa_adj) <- colnames(netgsa_adj) <- path_j
      } else if (length(setdiff(el_genenames, path_j))>0){
        g1 <- induced_subgraph(g, intersect(el_genenames, path_j))
        genes2add <- setdiff(path_j, igraph::V(g1)$name)
        if (length(genes2add)>0){
          g2 <- add_vertices(g1, length(genes2add))
          igraph::V(g2)$name[seq(1+length(igraph::V(g1)), length(igraph::V(g2)))] <- genes2add
          g1 <- g2
        }
        ## We redefine the adjacency matrix as some edges have multiplicity greater than 1.
        netgsa_adj <- as.matrix(as_adjacency_matrix(g1, type="both"))
        netgsa_adj[(netgsa_adj>0)] <- 1
        diag(netgsa_adj) <- 0
        
        ## We need to re-order the variables to match the order in the data. 
        new.order <- match(path_j, rownames(netgsa_adj))
        netgsa_adj <- netgsa_adj[new.order, new.order]
        # identical(colnames(netgsa_adj), path_j)

      }
      return(netgsa_adj)
    }
    adj <- lapply(path_membership, adj_single)
  }
  
  if (method=="topologyGSA"){
    pathways <- convertIdentifiers(pathways, "symbol")
    pathways <- lapply(pathways, pathwayGraph)
  }
  
  if (method=="camera"){
    pMembership <- lapply(pathways, function(a) 1*(match(genes$EntrezID, graph::nodes(a), nomatch = 0)>0))
    pMembership <- lapply(pMembership, function(a) {names(a)=genes$EntrezID;a})
    path_membership <- lapply(pMembership, function(a) which(a==1))
    names(path_membership) <- names(pathways)
  }
  if (method=="CePa"){
    # Here the path.catalogue only catelogues interactions,
    # pathways are listed with interaction id
    edgeList <- lapply(pathways, graph::edges)
    names(edgeList) <- NULL
    interactionList <- lapply(edgeList, function(a) data.frame('input'=paste0("ENTREZID:", a$src), 'output'=paste0("ENTREZID:",a$dest), stringsAsFactors=FALSE))
    interactionList <- lapply(1:length(interactionList), function(j) 
      data.frame('interaction.id'=paste(j, seq(1, nrow(edgeList[[j]])), sep='000'),interactionList[[j]], stringsAsFactors=FALSE))
    pathList <- lapply(interactionList, function(a) a$interaction.id)
    names(pathList) <- names(pathways)
    
    interactionList <- do.call(rbind, interactionList)

    paths <- try(sapply(pathways, convertIdentifiers, "symbol"), silent=TRUE)
    if (class(paths)=="try-error") stop ("pathway identifiers could not be converted")
    mapping <- lapply(1:N, function(a) data.frame('node.id'=graph::nodes(pathways[[a]]), 'symbol'=graph::nodes(paths[[a]]), stringsAsFactors = FALSE))
    
    pathways <- set.pathway.catalogue(pathList = pathList,interactionList = interactionList, mapping = do.call(rbind, mapping))

  }
  if (method=="DEGraph"){
    pe_kpg <- lapply(pathways, graphite::pathwayGraph)
    pathways <- lapply(1:N, function(j) graph::subGraph(intersect(genes$EntrezID, graph::nodes(pathways[[j]])), pe_kpg[[j]]))
  }
  
  if (method=="PathNet"){
    entrez <- gsub("ENTREZID:","",genes$EntrezID)
    edgeList <- lapply(pathways, graph::edges)
    pathways4PathNet <- lapply(1:N, function(j) data.frame(edgeList[[j]], "title"=names(pathways)[j]))
    pathways4PathNet <- do.call(rbind, pathways4PathNet)

    pathways <- pathways4PathNet[,c("src","dest","title")]

    ## The adjacency matrix is generated by pooling all interactions (directed or undirected) among used pathways, as described in the manual of PathNet. 
    g <- graph_from_edgelist(cbind(pathways$src, pathways$dest))
    adj <- as.matrix(get.adjacency(g, type="both"))
    adj[(adj>0)] <- 1
    adj <- adj[match(entrez, rownames(adj), nomatch = 0),  match(entrez, rownames(adj), nomatch = 0)]
  }
  
  return(list(pathways=pathways, membership=path_membership, adjacency=adj))
}

### Purpose: to compute the sensitivity and specificity of pathway enrichment.
#' @param  Ahat: a vector of 0-1 that indicates whether a pathway is selected for enrichment.
#' @param Amat: the corresponding true indicator of pathway enrichment.
#' @return A list with
#' \item{FPrate}{}
#' \item{FNrate}{} 
#' \item{sensitivity}{}
StructDiff <- function(Ahat, Amat, eps = 1e-06) {
  TP <- sum((abs(Ahat) > eps) * (abs(Amat) > eps), na.rm = T)
  TN <- sum((abs(Ahat) <= eps) * (abs(Amat) <= eps), na.rm = T)
  FP <- sum((abs(Ahat) > eps) * (abs(Amat) <= eps), na.rm = T)
  FN <- sum((abs(Ahat) <= eps) * (abs(Amat) > eps), na.rm = T)
  
  P <- TP + FN
  N <- TN + FP
  TPrate <- TP / (P + eps) ## Recall, or sensitivity
  TNrate <- TN / (N + eps) ## specificity
  FPrate <- FP / (N + eps)
  FNrate <- FN / (P + eps)
  
  Re <- TPrate
  Pr <- TP / (TP + FP + eps) ## Precision
  
  F1 <- (2 * Pr * Re) / (Pr + Re + eps)
  
  dev <- list(
    sensitivity = TPrate,
    specificity = TNrate,
    FPrate = FPrate,
    FNrate = FNrate,
    Precision = Pr,   
    TP = TP,
    FP = FP,
    TN = TN,
    FN = FN,
    F1 = F1
  )
  return(dev)
}

#' ----Functions for Results Summary-----
#' Function to retrieve power from multiple replications
# retrieve_power_cepa_GSA <- function(sigInd){
#   pvals <- lapply(sigInd, function(a) apply(a[,11:16],1,min))
#   pvals.fdr <- lapply(pvals, function(a) (p.adjust(a, "BH")<0.05))
#   base::Reduce("+",pvals.fdr)/length(pvals.fdr)
# }

retrieve_power_cepa <- function(df){
  pvals <- lapply(df, function(a) apply(a,1,min))
  pvals.fdr <- lapply(pvals, function(a) (p.adjust(a, "BH")<0.05))
  base::Reduce("+",pvals.fdr)/length(pvals.fdr)
}

# retrieve_power_cepa_ORA <- function(sigInd){
#   pvals <- lapply(sigInd, function(a) apply(a[,5:10],1,min))
#   pvals.fdr <- lapply(pvals, function(a) (p.adjust(a, "BH")<0.05))
#   base::Reduce("+",pvals.fdr)/length(pvals.fdr)
# }

retrieve_power <- function(expr, dL){
  j <- which(colnames(dL[[1]])==expr)
  pvals.fdr <- lapply(dL, function(a) (p.adjust(a[,j], "BH")<0.05)*1)
  base::Reduce("+",pvals.fdr)/length(pvals.fdr)
}

#' @param id indicator for which replications to include
retrieve_power_ora <- function(expr, dL, id){
  j <- which(colnames(dL[[1]])==expr)
  pvals.fdr <- lapply(dL[id==1], function(a) (p.adjust(a[,j], "BH")<0.05)*1)
  base::Reduce("+",pvals.fdr)/length(pvals.fdr)
}

print_power_w_filenames <- function(path, filenames, silent=TRUE){
  res <- lapply(paste0(path,filenames),function(a) {
    load(a); res
  }
  )
  
  pe.status <- unlist(lapply(res, function(a) a$pe))
  spia.status <- unlist(lapply(res, function(a) a$spia))
  prs.status <- unlist(lapply(res, function(a) a$prs))
  cepa.status <- unlist(lapply(res, function(a) a$cepa))
  
  if (class(res[[1]]$sigInd) == 'matrix'){
    sigInd <- lapply(res, function(a) a$sigInd)
  } else {# lists
    sigInd <- unlist(lapply(res, function(a) a$sigInd), recursive = FALSE) 
  }
  
  # If the total length is over 200, we only take the first 200 replications.
  if (length(sigInd)>200){
    sigInd <- sigInd[1:200]
    pe.status <- pe.status[1:200]
    spia.status <- spia.status[1:200]
    prs.status <- prs.status[1:200]
    cepa.status <- cepa.status[1:200]
  }
  if (!silent){
    cat('Total replications for ORA is...', sum(spia.status), '/', length(spia.status), '..\n')
  }
  ## check col names
  col23 <- sapply(sigInd, function(a) colnames(a)[23])
  if (col23[1]=="netgsa-joint"){
    sigInd <- lapply(sigInd, function(a) {colnames(a)[23]="netgsa.joint"; a})
  }
  
  methods <- c("netgsa", "netgsa.joint","DEGraph","camera","PathNet","topoGSA.mean","PE.ALL", "CePa.GSA", "CePa.ORA")
  df4power <- matrix(NA, nrow(sigInd[[1]]), length(methods))
  colnames(df4power) <- methods
  rownames(df4power) <- rownames(sigInd[[1]])
  df4power <- as.data.frame(df4power, stringsAsFactors=FALSE)
  for (a in methods[1:7]){
    df4power[,which(colnames(df4power) == a)] <- retrieve_power(a, sigInd)
  }
  
  df4power$CePa.GSA <- retrieve_power_cepa(lapply(sigInd, function(a) a[,grep("cepa.GSA",colnames(sigInd[[1]]))]))
  
  ## Return powers for all ORA methods if the status is nonzero
  df4power$PE <- rep(NA, length(df4power$netgsa))
  df4power$SPIA <- rep(NA, length(df4power$netgsa))
  df4power$PRS <- rep(NA, length(df4power$netgsa))
  
  if (sum(pe.status)>length(spia.status)/2){
    df4power$PE <- retrieve_power_ora("PE",sigInd, pe.status)
  }
  if (sum(spia.status)>length(spia.status)/2){
    df4power$SPIA <- retrieve_power_ora("spia",sigInd, pe.status)
  }
  if (sum(prs.status)>length(spia.status)/2){
    df4power$PRS <- retrieve_power_ora("PRS",sigInd, pe.status)
  }
  if (sum(cepa.status)>length(spia.status)/2){
    df4power$CePa.ORA <- retrieve_power_cepa(lapply(sigInd, function(a) a[,grep("cepa.ORA",colnames(sigInd[[1]]))]))
  }
  colnames(df4power)[which(colnames(df4power)=='netgsa')] <- "NetGSA2"
  colnames(df4power)[which(colnames(df4power)=="netgsa.joint")] <- "NetGSA"
  colnames(df4power)[which(colnames(df4power)=="camera")] <- "CAMERA"
  colnames(df4power)[which(colnames(df4power)=="topoGSA.mean")] <- "topologyGSA"
  colnames(df4power)[which(colnames(df4power)=="CePa.ORA")] <- "CePa"
  colnames(df4power)[which(colnames(df4power)=="PE")] <- "PE.Cut"
  colnames(df4power)[which(colnames(df4power)=="PE.ALL")] <- "PE.noCut"
  
  out <- list()
  out$power <- df4power
  
  df4power[is.na(df4power)] <- 0
  rank4power <- t(apply(df4power, 1, foo))
  colnames(rank4power) <- colnames(df4power)
  out$ranking <- rank4power
  
  return(out)
}


print_power_w_filenames_metabs <- function(path, filenames, silent=TRUE){
  res <- lapply(paste0(path,filenames),function(a) {
    load(a); res
  }
  )
  
  cepa.status <- unlist(lapply(res, function(a) a$cepa))
  
  if (class(res[[1]]$sigInd) == 'matrix'){
    sigInd <- lapply(res, function(a) a$sigInd)
  } else {# lists
    sigInd <- unlist(lapply(res, function(a) a$sigInd), recursive = FALSE) 
  }
  
  # If the total length is over 200, we only take the first 200 replications.
  # if (length(sigInd)>200){
  #   sigInd <- sigInd[1:200]
  #   cepa.status <- cepa.status[1:200]
  # }
  if (!silent){
    cat('Total replications for ORA is...', sum(cepa.status), '/', length(cepa.status), '..\n')
    }
  
  methods <- c("NetGSA", "NetGSA-joint","DEGraph","camera","PathNet", "CePa.GSA", "CePa.ORA")
  df4power <- matrix(NA, nrow(sigInd[[1]]), length(methods))
  colnames(df4power) <- methods
  rownames(df4power) <- rownames(sigInd[[1]])
  
  df4power[,which(colnames(df4power) %in% methods[1:5])] <- do.call(cbind,lapply(methods[1:5], function(a) retrieve_power(a, sigInd)))
  df4power <- as.data.frame(df4power, stringsAsFactors=FALSE)
  df4power$CePa.GSA <- retrieve_power_cepa(lapply(sigInd, function(a) a[,grep("cepa.GSA",colnames(sigInd[[1]]))]))
  
  ## Return powers for all ORA methods if the status is nonzero
  if (sum(cepa.status)>length(cepa.status)/2){
    df4power$CePa.ORA <- retrieve_power_cepa(lapply(sigInd, function(a) a[,grep("cepa.ORA",colnames(sigInd[[1]]))]))
  }
  colnames(df4power)[which(colnames(df4power)=="NetGSA")] <- "NetGSA2"
  colnames(df4power)[which(colnames(df4power)=="NetGSA-joint")] <- "NetGSA"
  colnames(df4power)[which(colnames(df4power)=="camera")] <- "CAMERA"
  colnames(df4power)[which(colnames(df4power)=="CePa.ORA")] <- "CePa"
  
  out <- list()
  out$power <- df4power
  
  df4power[is.na(df4power)] <- 0
  rank4power <- t(apply(df4power, 1, foo))
  colnames(rank4power) <- colnames(df4power)
  out$ranking <- rank4power
  
  return(out)
}

foo <- function(x){
  match(x, sort(unique(x), decreasing = T))
}

geometric.mean <- function (x, na.rm = TRUE) 
{
  if (is.null(nrow(x))) {
    exp(mean(log(x), na.rm = TRUE))
  }
  else {
    exp(apply(log(x), 2, mean, na.rm = na.rm))
  }
}

#' @param dfList A list of data matrices to reshape
reshapeDF <- function(dfList, var.name="power"){
  K <- length(dfList)
  
  outDF <- data.frame('mu'=0.1, reshape2::melt(dfList[['mu1']], id.vars = 1:3, value.name = "power"), stringsAsFactors = F)
  if (K>1){
    for (meanLevel in paste0("mu",seq(2,K))){
      outDF <- rbind(outDF, data.frame('mu'=as.numeric(gsub("mu","",meanLevel))/10, reshape2::melt(dfList[[meanLevel]], id.vars = 1:3, value.name = "power"), stringsAsFactors = F))
    } 
    outDF$name_f <- paste0(outDF$name, " (size=",outDF$size,", DC=", round(outDF$DC,2),")")
    DF.levels <- levels(factor(outDF$name_f))
    df <- list()
    df$size <- as.numeric(sapply(sapply(DF.levels, function(a) strsplit(a,"=")[[1]][2]), function(a) gsub(", DC","",a)))
    df$DC <- as.numeric(sapply(sapply(DF.levels, function(a) strsplit(a,"=")[[1]][3]), function(a) gsub(")","",a)))
    outDF$name_g <- factor(outDF$name_f, levels=DF.levels[with(df, order(DC, size))])
    outDF <- outDF[,-which(colnames(outDF)=="name_f")]
  }
  colnames(outDF)[which(colnames(outDF)=="variable")] <- "method"
  if (var.name=="rank"){
    colnames(outDF)[which(colnames(outDF)=="power")] <- "rank"
  }
  outDF
}

#' @param DEBUG which model
#' @param perm whether there is permutation of sample labels
power.plot <- function(DEBUG, cancerType="BCA", silent=TRUE){

  if (!silent){
    cat("Model is ... ", DEBUG, "...\n")
    cat("Cancer type is ... ", cancerType, "...\n")
  }
  
  cancerName <- ifelse(cancerType=="BCA", "breastcancer2012", "prostatecancer2015")
  load(paste0("../6_Results_New/",cancerType,"/", cancerName,"_ready.rda"))
  size <- sapply(netgsa_paths$membership, length)
  
  pathway_group_propAlt_indicator <- rep(1, length(size))
  empirical_isEnrichment <- rep(0, length(size))
  empirical_propAlt <- rep(0, length(size))
  for (j in 1:length(size)){
    empirical_propAlt[j] <- sum(graph::nodes(deGraph_paths$pathways[[j]]) %in% base::eval(as.name(paste0("genes2affect_",DEBUG)))) / size[j]
    empirical_isEnrichment[j] <- ifelse(empirical_propAlt[j]>0, 1, 0)
  }
  
  path <- paste("../6_Results_New", cancerType, "results/mu0/", sep='/')
  filenames <- list.files(path)
  t1error <- print_power_w_filenames(path, filenames)$power
  
  sigInd <- list()
  sigInd_perm <- list()
  rankMat <- list()
  rankMat_perm <- list()
  
  for (meanLevel in paste0("mu",seq(1,5))){  
    if (!silent){    
      cat("mean level is ... ", meanLevel, "...\n")
    }
    
    path <- paste0("../6_Results_New/", cancerType, "/results/", meanLevel, "/")
    filenames <- list.files(path)
    
    tmp <- grep(paste(cancerName, DEBUG, meanLevel, sep="_"), filenames)
    obj <- print_power_w_filenames(path, filenames[tmp])
    sigInd[[meanLevel]] <- obj$power
    rankMat[[meanLevel]] <- obj$ranking
    
    tmp <- grep(paste(cancerName, DEBUG,  'perm', meanLevel, sep="_"), filenames)
    obj <- print_power_w_filenames(path, filenames[tmp])
    sigInd_perm[[meanLevel]] <- obj$power
    rankMat_perm[[meanLevel]] <- obj$ranking
  }
  
  ## data frame for t1error
  t1error <- data.frame('name'=rownames(t1error), 'size'=size, "DC"=empirical_propAlt, t1error, stringsAsFactors=FALSE)
  rownames(t1error) <- NULL
  df4t1error <- data.frame('mu'=0, reshape2::melt(t1error, id.vars = 1:3, value.name = "power"), stringsAsFactors = F)
  colnames(df4t1error)[which(colnames(df4t1error)=="variable")] <- "method"
  
  out <- list()
  out$t1error <- df4t1error
  out$power <- sigInd
  out$perm_power <- sigInd_perm
  
  for (meanLevel in paste0("mu",seq(1,5))){
    rankMat[[meanLevel]] <- data.frame('name'=rownames(rankMat[[meanLevel]]),
                                       'size'=size,  "DC"=empirical_propAlt,
                                       rankMat[[meanLevel]],
                                       stringsAsFactors=FALSE)
    rownames(rankMat[[meanLevel]]) <- NULL
    
    sigInd[[meanLevel]] <- data.frame('name'=rownames(sigInd[[meanLevel]]),
                                       'size'=size,  "DC"=empirical_propAlt,
                                       sigInd[[meanLevel]],
                                       stringsAsFactors=FALSE)
    rownames(sigInd[[meanLevel]]) <- NULL
    
    
    rankMat_perm[[meanLevel]] <- data.frame('name'=rownames(rankMat_perm[[meanLevel]]),
                                            'size'=size,  "DC"=empirical_propAlt,
                                            rankMat_perm[[meanLevel]],
                                            stringsAsFactors=FALSE)
    rownames(rankMat_perm[[meanLevel]]) <- NULL
    
    sigInd_perm[[meanLevel]] <- data.frame('name'=rownames(sigInd_perm[[meanLevel]]),
                                            'size'=size,  "DC"=empirical_propAlt,
                                            sigInd_perm[[meanLevel]],
                                            stringsAsFactors=FALSE)
    rownames(sigInd_perm[[meanLevel]]) <- NULL
  }
  
  BigDF.p0 <- reshapeDF(rankMat, "rank")
  BigDF.p1 <- reshapeDF(rankMat_perm, "rank")
  BigDF.final <- data.frame("perm"="FALSE", BigDF.p0)
  BigDF.final <- rbind(BigDF.final, data.frame("perm"="TRUE", BigDF.p1))

  out$rank <- rankMat
  out$perm_rank <- rankMat_perm
  out$DF4rank <- BigDF.final
  
  BigDF.p0 <- reshapeDF(sigInd, "power")
  BigDF.p1 <- reshapeDF(sigInd_perm, "power")
  BigDF.final <- data.frame("perm"="FALSE", BigDF.p0)
  BigDF.final <- rbind(BigDF.final, data.frame("perm"="TRUE", BigDF.p1))
  
  out$DF4power <- BigDF.final
  
  return(out)
}

