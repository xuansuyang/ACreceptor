#activity prediction with the linear regression model(human)
#10/30/24
#Created by Suyang Xuan

####################################################
#Function
####################################################
#count2cpm
cpm=function(expr){
  apply(expr,2,function(x){x/sum(x)*1000000})
}
#another count2cpm function
cpm_normalized=function(counts.matrix,filter=T,log_normalize=T){
  loggeomeans <- rowMeans(log(counts.matrix))
  sf=apply(counts.matrix, 2, function(cnts) {
    exp(stats::median((log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0]))
  })
  normalized.matrix=t( t( counts.matrix ) / sf )
  normalized.cpm=apply(normalized.matrix,2,function(x){x/sum(x)*1000000})
  if(filter){normalized.cpm <- normalized.cpm[rowMeans(normalized.cpm) > 0.25, ]}
  if(log_normalize){normalized.cpm=log(normalized.cpm+0.001)}
  normalized.cpm
}

activity<-function(expr,model,scale=TRUE){
  common_genes <- intersect(rownames(expr), rownames(model))
  result <- t(expr[common_genes,,drop=FALSE]) %*% as.matrix(model[common_genes,,drop=FALSE])
  if (nrow(result) > 1 & scale==TRUE) {
    rn <- rownames(result)
    result <- apply(result, 2, scale)
    rownames(result) <- rn
  }
  return(t(result))
}
####################################################
#Input and output
####################################################
#inout:"output/model.rds","counts.txt"
#output:path activity
#load model
model=readRDS("data/human_model.rds")

#an example
expr=read.table("input/counts/GSE143989_counts.txt")
expr=cpm(expr)
expr=cpm_normalized(expr)
ac=activity(expr,model,scale = F)
ac=activity(expr,model)
