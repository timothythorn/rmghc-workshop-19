read_feature_counts <- function (countDir, sampleID, suffix=".fCounts") {
  
  file_list <- list.files(countDir)
  
  for (i in 1:length(sampleID)) {
    filename <- file_list[grep(sampleID[i], file_list)]
    filename <- filename[grep(paste0("*",suffix,"$"), filename)]
    filename=paste0(countDir,"/",filename)
    print(filename)
    temp=read.table(filename,header=T,stringsAsFactors = F)
    if (i==1) {
      annot=temp[,c(1,5,6)]
      fCount=data.frame(temp[,7])
    } else {
      fCount=data.frame(fCount,temp[,7])
    }
  }
  colnames(annot)=c("geneID","strand","length")
  annot$strand=substr(annot$strand,1,1)
  
  colnames(fCount)=sampleID
  rownames(fCount)=annot$geneID
  return(list(count=fCount,annot=annot))
}

convert_counts <- function(fCount,return="TPM") {
  if (return=="FPKM") {
    out=fCount$count/(fCount$annot$length/1000)
    out=t(t(out)/(colSums(fCount$count)/10^6))
  } else if (return=="TPM") {
    out=fCount$count/(fCount$annot$length/1000)
    out=t(t(out)/(colSums(out)/10^6))
  }
  return(out)
}


extract_attr <- function(x, attr) {
  attr_list <- unlist(strsplit(x, ";"))
  attr_index <- grep(attr, attr_list)
  attr_val <- NA
  if(length(attr_index) != 0) {
    attr_val <- attr_list[attr_index]
    attr_val <- gsub(paste0(attr," "), "", attr_val)
    attr_val <- gsub("\"", "", attr_val)
    attr_val <- gsub(" ", "", attr_val)
  }
  return(attr_val)
}

get_gene_name <- function(gene_id) {
  gen[which(gen$gene_id == gene_id),"gene_name"]
}

get_gene_id <- function(gene_name) {
  gen[which(gen$gene_name == gene_name),"gene_id"]
}


plotPCA <- function (expr, sampleSheet, colorVar=NA, shapeVar=NA, pointSize=5) {
  require("ggplot2")
  pcaObj<-prcomp(t(expr[which(apply(expr,1,var)>0),]), center=T, scale=T)
  pcaRes=data.frame(sampleSheet,PC1=pcaObj$x[,1],PC2=pcaObj$x[,2])
  pcaVarsPct=signif(((pcaObj$sdev)^2)/(sum((pcaObj$sdev)^2)),3)*100
  
  if (!is.na(colorVar) & !is.na(shapeVar)) {
    p=ggplot(pcaRes, aes_string("PC1", "PC2", color=colorVar, shape=shapeVar))
  }  
  
  if (!is.na(colorVar) & is.na(shapeVar)) {
    p=ggplot(pcaRes, aes_string("PC1", "PC2", color=colorVar))
  }  
  
  if (is.na(colorVar) & !is.na(shapeVar)) {
    p=ggplot(pcaRes, aes_string("PC1", "PC2", shape=shapeVar))
  }  
  
  if (is.na(colorVar) & is.na(shapeVar)) {
    p=ggplot(pcaRes, aes_string("PC1", "PC2"))
  }  
  
  
  p + scale_shape_manual(values=seq(1, sampleSheet %>% pull(shapeVar) %>% unique %>% length, by = 1)) + 
    geom_point(size=pointSize) +
    xlab(paste0("PC1: ",pcaVarsPct[1],"% variance")) + ylab(paste0("PC2: ",pcaVarsPct[2],"% variance"))
  
}


kMeansPP <- function(df, k, doPlot = TRUE){
  kCenters <- data.frame(matrix(NA, ncol = ncol(df), nrow = k))
  whichPoints <- rep(NA, k)
  whichPoints[1] <- sample(1:nrow(df), 1)
  kCenters[1, ] <- df[whichPoints[1], ]  # Initial center
  
  for(kk in 2:k){
    distMat <- proxy::dist(df, kCenters[1:(kk-1), ])
    distToNearestCenter <- apply(distMat, 1, min)
    whichPoints[kk] <- sample(1:nrow(df), 1, prob = distToNearestCenter^2)
    kCenters[kk, ] <- df[whichPoints[kk], ]
  }
  
  if(doPlot == TRUE){
    plot(df, col = "GRAY")
    points(kCenters, col = 1:k, pch = 20)
  }
  
  outList <- NULL
  outList$Centers <- kCenters
  outList$whichPoints <- whichPoints
  return(outList) 
}
