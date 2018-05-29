library(ppcor)
library(gdata)
library(graphite)
library(graph)


#This function drops metabolites based on whether metab.id is TRUE or FALSE
drop.metabs = function(dataset, metab.id = NULL){
  if(is.null(metab.id)){print("no metabolite selected")}
  else{
    dataset[["metab_info"]] = dataset[["metab_info"]][!metab.id,]
    rownames(dataset[["metab_info"]]) <- NULL
    dataset[["dat"]] = dataset[["dat"]][!metab.id,]
  }
  return(dataset)
}

drop.genes = function(dataset, gene.id = NULL){
  if(is.null(gene.id)){print("no gene selected")}
  else{
    dataset[["gene_info"]] = dataset[["gene_info"]][!gene.id,]
    rownames(dataset[["gene_info"]]) <- NULL
    dataset[["dat"]] = dataset[["dat"]][!gene.id,]
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


## Community pathway deregulation
## For a given pathway, find communities and choose a subset of communities such that
## genes in these communities are affected with the desired mean changes. 
community.dereg <- function(pathway, genes, DC=0.5, maxIter=10, seed=1){
  # pathway: a pathway object
  # genes: the genes that are used in the enrichment analysis. This can be learned from the B matrix.
  # DC: detection call, a percentage between 0% and 100% indicating the proportiontion of genes to be affected
  #  In this design, DC is fixed to be about 50% to make it easier for control.
  set.seed(seed)
  el <- edges(pathway)
  gr <- graph_from_edgelist(cbind(el$src, el$dest), directed = F)
  
  # get ride of duplicated edges
  adj <- as.matrix(get.adjacency(gr, type="both"))
  adj[(adj>0)] = 1;  diag(adj) <- 0;
  gr <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
  gr <- induced_subgraph(gr, which(match(V(gr)$name,genes,nomatch = 0)>0))

  size <- length(genes)
  
  lec <- cluster_edge_betweenness(gr)
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
  res <- list(lec=lec$membership, GeneID=subset, realDC = length(subset)/size) 
  
  return(res)
}
  

betweenness.dereg <- function(pathway, genes, DC=0.5, seed=1){
  set.seed(seed)
  size <- length(genes)
  el <- edges(pathway)
  gr <- igraph::graph_from_edgelist(cbind(el$src, el$dest), directed = F)
  
  # get ride of duplicated edges
  adj <- as.matrix(get.adjacency(gr, type="both"))
  adj[(adj>0)] = 1;   diag(adj) <- 0;
  gr <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
  
  #subgraph based on the variables in genes
  gr <- induced_subgraph(gr, which(match(V(gr)$name,genes,nomatch = 0)>0))
  
  ## calculate edge betweenness based on gr
  vb <- betweenness(gr, directed = F)

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
  el <- edges(pathway)
  gr <- igraph::graph_from_edgelist(cbind(el$src, el$dest), directed = F)
  
  # get ride of duplicated edges
  adj <- as.matrix(get.adjacency(gr, type="both"))
  adj[(adj>0)] = 1
  diag(adj) <- 0
  gr <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
  
  #subgraph based on the variables in genes
  gr <- induced_subgraph(gr, which(match(V(gr)$name,genes,nomatch = 0)>0))
  ## calculate edge betweenness based on gr
  ngbh <- ego(gr, nei)
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
