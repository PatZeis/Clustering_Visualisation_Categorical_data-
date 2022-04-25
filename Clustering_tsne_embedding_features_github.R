library(cluster)
library(fpc)
library(tsne)
library(Rtsne)
library(stringr)
library(plyr) 
library(umap)
library(feather)

## class definition 
SCseq <- setClass("SCseq", slots = c(surveydata = "data.frame", distances = "matrix", umap = "data.frame", tsne = "data.frame", cluster = "list",  cpart = "vector"))

setValidity("SCseq",
            function(object) {
              msg <- NULL
              if ( ! is.data.frame(object@surveydata) ){
                msg <- c(msg, "input data must be a data.frame")
              }else if ( nrow(object@surveydata) < 2 ){
                msg <- c(msg, "input data must have more than one row")
              }else if ( ncol(object@surveydata) < 2 ){
                msg <- c(msg, "input data must have more than one column")
              }
              if (is.null(msg)) TRUE
              else msg
            }
)

setMethod("initialize",
          signature = "SCseq",
          definition = function(.Object, surveydata ){
            .Object@surveydata <- surveydata
            validObject(.Object)
            return(.Object)
          }
)

## required functions for clustering and dimensional reduction

clusGapExt <- function (x, FUNcluster, K.max, B = 100, verbose = TRUE, method="euclidean",random=TRUE,diss=FALSE,
                       ...) 
{
  stopifnot(is.function(FUNcluster), length(dim(x)) == 2, K.max >= 
              2, (n <- nrow(x)) >= 1, (p <- ncol(x)) >= 1)
  if (B != (B. <- as.integer(B)) || (B <- B.) <= 0) 
    stop("'B' has to be a positive integer")
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  ii <- seq_len(n)
  W.k <- function(X, kk) {
    clus <- if (kk > 1) 
      FUNcluster(X, kk, ...)$cluster
    else rep.int(1L, nrow(X))
    0.5 * sum(vapply(split(ii, clus), function(I) {
      if ( diss ){
        xs <- X[I,I, drop = FALSE]
        sum(xs/nrow(xs))
      }else{
        xs <- X[I, , drop = FALSE]
        d <- dist.gen(xs,method=method)
        sum(d/nrow(xs))
      }
    }, 0))
  }
  logW <- E.logW <- SE.sim <- numeric(K.max)
  if (verbose) 
    cat("Clustering k = 1,2,..., K.max (= ", K.max, "): .. \n", 
        sep = "")
  for (k in 1:K.max){
    if (verbose) cat("k =",k,"\r")
    logW[k] <- log(W.k(x, k))
  }
  if (verbose){ 
    cat("\n")
    cat("done.\n")
  }
  if (random){
    xs <- scale(x, center = TRUE, scale = FALSE)
    m.x <- rep(attr(xs, "scaled:center"), each = n)
    V.sx <- svd(xs)$v
    rng.x1 <- apply(xs %*% V.sx, 2, range)
    logWks <- matrix(0, B, K.max)
    
    if (verbose) 
      cat("Bootstrapping, b = 1,2,..., B (= ", B, ")  [one \".\" per sample]:\n", 
          sep = "")
    for (b in 1:B) {
      z1 <- apply(rng.x1, 2, function(M, nn) runif(nn, min = M[1], 
                                                   max = M[2]), nn = n)
      z <- tcrossprod(z1, V.sx) + m.x
      ##z <- apply(x,2,function(m) runif(length(m),min=min(m),max=max(m)))
      ##z <- apply(x,2,function(m) sample(m))
      for (k in 1:K.max) {
        logWks[b, k] <- log(W.k(z, k))
      }
      if (verbose) 
        cat(".", if (b%%50 == 0) 
          paste(b, "\n"))
    }
    if (verbose && (B%%50 != 0)) 
      cat("", B, "\n")
    E.logW <- colMeans(logWks)
    SE.sim <- sqrt((1 + 1/B) * apply(logWks, 2, var))
  }else{
    E.logW <- rep(NA,K.max)
    SE.sim <- rep(NA,K.max)
  }
  structure(class = "clusGap", list(Tab = cbind(logW, E.logW, 
                                                gap = E.logW - logW, SE.sim), n = n, B = B, FUNcluster = FUNcluster))
}
  

pamkdCBI <- function (data, krange = 2:10, k = NULL, criterion = "asw", usepam = TRUE, 
                      scaling = TRUE, diss = inherits(data, "dist"), ...) 
{
  if (!is.null(k)) 
    krange <- k
  c1 <- pamk(as.dist(data), krange = krange, criterion = criterion, 
             usepam = usepam, scaling = scaling, diss = diss, ...)
  partition <- c1$pamobject$clustering
  cl <- list()
  nc <- c1$nc
  
  for (i in 1:nc) cl[[i]] <- partition == i
  out <- list(result = c1, nc = nc, clusterlist = cl, partition = partition, 
              clustermethod = "pam/estimated k", criterion = criterion)
  out
}


## clustering and dimensional reduction 
cluster_survey <- function(object, cln = NULL, tsne_perplexity=30, umap.pars = umap.defaults) {
  dismat <- daisy(object@surveydata, metric = "gower", type = list("ordratio"=ints))
  object@distances <- as.matrix(dismat)
  gpr <- clusGapExt(object@distances, FUN = function(x,k) pam(as.dist(x),k), K.max = 30, B = 50, random=FALSE, method="gower",diss=TRUE)
  g <- gpr$Tab[,1] ## within cluster dispersion
  y <- g[-length(g)] - g[-1] ### change of within cluster disperion from 1 cluster to the other
  mm <- numeric(length(y))  
  nn <- numeric(length(y))
  for ( i in 1:length(y)){
    mm[i] <- mean(y[i:length(y)]) ## average wthin cluser dispersion by increasing cluster size
    nn[i] <- sqrt(var(y[i:length(y)])) ## standart deviation of within cluster dispersion by increasing cluster size
  }
  if ( is.null(cln)) {
  cln <- max(min(which( y - (mm + nn) < 0 )),1) ## selection of cluster number based on saturation of within cluster dispersion
  }
  if ( cln < 2 ) cln <- 2
  else { cln <- cln}
  clb <- clusterboot(object@distances,B=50,bootmethod="boot",clustermethod=pamkdCBI,scaling=FALSE,diss=TRUE,k=cln,multipleboot=FALSE,bscompare=TRUE,seed=17000)
  ts <-  Rtsne(as.dist(object@distances),dims=2,initial_config=cmdscale(as.dist(object@distances),k=2),perplexity=tsne_perplexity)$Y
  object@tsne <- as.data.frame(ts)
  rownames(object@tsne) <- colnames(sc@distances)
  umap.pars$input <- "dist"
  
  object@umap <- as.data.frame(umap(object@distances, config = umap.pars)$layout)
  rownames(object@umap) <- colnames(object@distances)
  object@cluster   <- list( kpart=clb$result$partition, jaccard=clb$bootmean, gap=gpr, clb=clb, features=colnames(object@surveydata))
  return(object)
}

# tSNE map visualisation 
plotclustertsne <- function(object, um=F) {
  if (um==T){
    d <- object@umap
  }
  else { d <- object@tsne}
  plot(d,xlab="",ylab="",pch=20,cex=1.5,col=adjustcolor("lightgrey",1),axes=FALSE)
  part <- object@cluster$kpart
  fcol <- rainbow(length(unique(part)))
  if ( um==F) {
  for ( i in 1:max(part) ){
    if ( sum(part == i) > 0 ) text(object@tsne[part == i,1],object@tsne[part == i,2],i,col=adjustcolor(fcol[i],1),cex=.75,font=4)
  }}
  else {
  for ( i in 1:max(part)) {  
    if (sum(part == i) > 0) text(object@umap[part == i,1], object@umap[part == i, 2], i, col=adjustcolor(fcol[i], 1), cex = .75, font = 4)
  }}
}

### z-scores for QA pairs per cluster
getzscore_cluster <- function(object, zsc=F, ad=T, trh = 0.05 ) {
  zscores <- list()
  tables <- list()
  for ( i in 1:length(unique(object@cluster$kpart))) {
    cat(paste(i, "\n"))
    answers <- NULL
    cl <- object@surveydata[object@cluster$kpart == i,]  
    answers <- apply(cl, 2, function(x) {
      y <- table(x)/nrow(cl)    ### normalize number of answers by number person in cluster
      
    })
    for ( n in 1:length(answers)) { ## collapse question and answers 
      names(answers[[n]]) <- sapply(names(answers[[n]]), function(x) {
        if ( length(ints) == 0) 
          stop("define integers indices")
        if ( n %in% ints){
        y <- paste(names(answers)[n], as.numeric(x), sep  = "_") 
        }
        else {y <- paste(names(answers)[n], x, sep  = "_")}
      })
    }
    clust <- unlist(answers)    ## vector of answers
    #if ( i == 1) {
    QA <- numeric()  ## name of answers 
    for ( n in 1:length(answers)) {
      QA <- append(QA, names(answers[[n]]))
    } 
    #} 
    clust <- as.numeric(clust) 
    names(clust) <- QA ## naming of answers 
    #f <- which(clust > trh)
    #clust <- clust[f]
    cat(paste(i, "\n"))
    tables[[i]] <- data.frame(names(clust), as.numeric(clust)) ## list with answer tables for cluster
    colnames(tables[[i]]) <- c("answers", paste(i, "clus", sep = "_"))
  }
  #return(tables)
  ## merging of answers from clusters
  for ( i in 1:length(tables)) data2 <- if ( i == 1) tables[[i]] else merge(data2, tables[[i]], by="answers", all = T)
  rownames(data2) <- data2[,1]
  data2 <- data2[,-1]
  for ( i in 1:ncol(data2)) {
    data2[is.na(data2[,i]),i] <- 0
  }
  return(data2)
  #data2 <- data2[apply(data2,1,max) >= trh, ]
  #return(data2)
  data3 <- t(apply(data2, 1, function(x){
    mu <- mean(x)
    sd <- sd(x)
    x <- (x-mu)/sd
  }))
  #return(data3)
  zscores <- alply(data3, 2, function(x) {
    x[order(x, decreasing = T)]
  })
  agree_disagree <- lapply(zscores, function(x) {
    y <- data.frame(rank=paste("rank", 1:100, sep = "_"), UPinCluster=head(names(x), 100), DOWNinCluster=tail(names(x), 100)[order(tail(x, 100))])
    #rownames(y) <- paste("rank", 1:nrow(y), sep = "_")
    return(y)
    })
  if ( ad == T) {
  return(agree_disagree)
  }
  if (zsc ==T) {
    return(zscores)
  }
}

## significant enrichment of p-values for QA pairs for clusters
pop <- length(sc@cluster$kpart)
clustersize <- as.numeric(table(sc@cluster$kpart))


for ( i in 1:length(unique(sc@cluster$kpart))){
  a <- apply(agreedisagree_absolute, 1, function(x) {
    b <- sum(x)
    fisher <- fisher.test(matrix(c(pop-clustersize[i], clustersize[i], b - x[i], x[i]), ncol = 2), alternative = "greater")
    pv <- fisher$p.value
    return(pv)
  })
  if( i == 1) {
    data_pval <- a
  }
  else {
    data_pval <- cbind(data_pval, a)
  }
}

pvalues <- list()

for ( i in 1:ncol(data_pval)) {
  pvalues[[i]] <- data_pval[,i][order(as.numeric(data_pval[,i]))]
}
pvalues_sig <- lapply(pvalues, function(x) {
  y <- x[as.numeric(x) <= 0.05]
  return(y)
})


for ( i in 1:length(pvalues_sig)) {
  write.table(data.frame(pvalues_sig[[i]]), file = paste("enrichment_pvalues_QA_pairs_in_cluster", i, ".txt"), row.names = T, col.names = T, sep = "\t", quote = F)
}

### Selection of questions and tSNE visualisation 
plotsymbolsmap <- function (object, types, subset = NULL, samples_col = NULL, cex = 0.5, 
           um = FALSE, leg = TRUE, map = TRUE, ordering=F, orders=NULL) 
{
  require(RColorBrewer)
  paired <- brewer.pal(9, "Set1")
  if (length(object@tsne) == 0  & length(object@umap) == 
      0) 
    stop("run comptsne/compumap before plotlabelsmap")
  if (!is.logical(um)) 
    stop("um has to be TRUE or FALSE")
  if ( um == FALSE & dim(object@tsne)[1] == 0) {
     if (dim(object@umap)[1] != 0) {
      um <- TRUE
    }
  }
  main <- types
  types <- as.character(object@surveydata[,types])
  if (is.null(subset)) 
    subset <- unique(types)
  h <- sort(unique(types)) %in% subset
  if (!is.null(subset)) {
    fp <- rep(FALSE, length(types))
    fp[types %in% subset] <- TRUE
  }
  if (is.null(samples_col)) {
    if ( !(length(unique(types[fp])) > 9))
    {
    samples_col <- paired[1:length(unique(types[fp]))]}
    else {
      qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
      col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      samples_col=sample(col_vector, length(unique(types[fp])))
    }
  }
  else {
    samples_col <- samples_col[h]
  }
  if (um) {
    d <- object@umap
  }
  else {
    d <- object@tsne
  }
  if (map) {
    plot(d, xlab = "", ylab = "", axes = FALSE, cex = cex, 
         pch = 20, col = "grey", main = main)
    for (i in 1:length(unique(types[fp]))) {
      f <- types == sort(unique(types[fp]))[orders][i]
      points(d[f, 1], d[f, 2], col = samples_col[i], pch = 20, 
             cex = cex)
    }
  }
  else {
    plot(d, xlab = "", ylab = "", axes = FALSE, cex = 0, 
         pch = 20, col = "grey", xlim = c(min(d[, 1]), max(d[, 
                                                             1])), ylim = c(min(d[, 2]), max(d[, 2])))
  }
  if (leg) 
    legend("topleft", legend = sort(unique(types[fp]))[orders], col = samples_col, 
           pch = 20, cex = 0.75, bty = "n")
}


