#calculate counts to Zscore for each perturbation
#10/30/24
#Created by Dekang Lv

#attach libraries


####################################################
#Function
####################################################

#use 2sd to reduce DEG of low expression
cpm2zscore<-function(cpms){
  control_group=str_detect(colnames(cpms),"control")
  perturbed_group=str_detect(colnames(cpms),"perturbed")
  control=cpms[,control_group,drop = FALSE]
  print(dim(control))
  perturbed=cpms[,perturbed_group,drop = FALSE]
  print(dim(control))
  mean_control <- apply(control, 1, mean, na.rm = TRUE)
  mean_perturbed <- apply(perturbed, 1, mean, na.rm = TRUE)
  logFC <- mean_perturbed - mean_control
  if(ncol(control)>1 & ncol(perturbed)>1){
    sd_control <- apply(control, 1, sd, na.rm = TRUE)
    model_control <- loess(sd_control ~ mean_control)
    sd_control = predict(model_control,mean_control)
    sd_perturbed <- apply(perturbed, 1, sd, na.rm = TRUE)
    model_perturbed <- loess(sd_perturbed ~ mean_perturbed)
    sd_perturbed = predict(model_perturbed,mean_perturbed)
    sd = sd_control+sd_perturbed
  }else{
    mean_all=apply(cpms, 1, mean, na.rm = TRUE)
    sd_all=apply(cpms, 1, sd, na.rm = TRUE)
    model_all=loess(sd_all ~ mean_all)
    sd=predict(model_all,mean_all)
  }
  zscore <- logFC / sd
  zscore
}

#method get better roc
cpm2zscore_v2<-function(cpms){
  control_group=str_detect(colnames(cpms),"control")
  perturbed_group=str_detect(colnames(cpms),"perturbed")
  control=cpms[,control_group,drop = FALSE]
  print(dim(control))
  perturbed=cpms[,perturbed_group,drop = FALSE]
  print(dim(perturbed))
  mean_control <- apply(control, 1, mean, na.rm = TRUE)
  mean_perturbed <- apply(perturbed, 1, mean, na.rm = TRUE)
  logFC <- mean_perturbed - mean_control

    sd_all=apply(cpms, 1, sd, na.rm = TRUE)
    model_all=loess(sd_all ~ mean_control)
    sd=predict(model_all,mean_perturbed)

  zscore <- logFC / sd
  zscore
}

#old method used by progeny
cpm2zscore_v1<-function(cpms){
  control_group=str_detect(colnames(cpms),"control")
  perturbed_group=str_detect(colnames(cpms),"perturbed")
  control=cpms[,control_group,drop = FALSE]
  print(dim(control))
  perturbed=cpms[,perturbed_group,drop = FALSE]
  print(dim(control))
  mean_control <- apply(control, 1, mean, na.rm = TRUE)
  mean_perturbed <- apply(perturbed, 1, mean, na.rm = TRUE)
  logFC <- mean_perturbed - mean_control
  if(ncol(control)>1 & ncol(perturbed)>1){
    sd_control <- apply(control, 1, sd, na.rm = TRUE)
    model_control <- loess(sd_control ~ mean_control)
    sd = predict(model_control,mean_perturbed)
    #sd_control = predict(model_control,mean_control)
    #sd_perturbed <- apply(perturbed, 1, sd, na.rm = TRUE)
    #model_perturbed <- loess(sd_perturbed ~ mean_perturbed)
    #sd_perturbed = predict(model_perturbed,mean_perturbed)
    #sd = sd_control+sd_perturbed
  }else{
    #mean_all=apply(cpms, 1, mean, na.rm = TRUE)
    sd_all=apply(cpms, 1, sd, na.rm = TRUE)
    model_all=loess(sd_all ~ mean_control)
    sd=predict(model_all,mean_perturbed)
  }
  zscore <- logFC / sd
  zscore
}

# count2zscore_v1<-function(counts){
#   control_group=str_detect(colnames(counts),"control")
#   perturbed_group=str_detect(colnames(counts),"perturbed")
#   group=rep(NA,ncol(counts))
#   group[control_group]="control"
#   group[perturbed_group]="perturbed"
#   group=factor(group)
#   group
#   keep.exprs<- edgeR::filterByExpr(counts,min.total.count = 10,group=group)
#   sum(keep.exprs)
#   counts <- counts[keep.exprs, ]
#   counts = log1p(counts)
#   control=counts[,control_group,drop = FALSE]
#   print(dim(control))
#   perturbed=counts[,perturbed_group,drop = FALSE]
#   print(dim(perturbed))
#
#   mean_control <- apply(control, 1, mean, na.rm = TRUE)
#   mean_perturbed <- apply(perturbed, 1, mean, na.rm = TRUE)
#   logFC <- mean_perturbed - mean_control
#   if(ncol(control)>1 & ncol(perturbed)>1){
#     sd_control <- apply(control, 1, sd, na.rm = TRUE)
#     model_control <- loess(sd_control ~ mean_control)
#     sd = predict(model_control,mean_perturbed)
#     #sd_control = predict(model_control,mean_control)
#     #sd_perturbed <- apply(perturbed, 1, sd, na.rm = TRUE)
#     #model_perturbed <- loess(sd_perturbed ~ mean_perturbed)
#     #sd_perturbed = predict(model_perturbed,mean_perturbed)
#     #sd = sd_control+sd_perturbed
#   }else{
#     #mean_all=apply(cpms, 1, mean, na.rm = TRUE)
#     sd_all=apply(counts, 1, sd, na.rm = TRUE)
#     model_all=loess(sd_all ~ mean_control)
#     sd=predict(model_all,mean_perturbed)
#   }
#   zscore <- logFC / sd
#   zscore
# }
# c2zscore.list=lapply(counts.list,count2zscore_v1)

#count2cpm function
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

####################################################
#Input and output
####################################################
#inout:"input/newtag.csv","input/counts/*_counts.txt"
#output:"output/counts_cpms_zscores.rda"


#setwd("/Users/dekanglv/0project/0receptor_activation_prodiction/")
#1.read relation of experiments,perturbed_genes and directions
newtag=read.csv("input/tag_final.csv",row.names = 1,stringsAsFactors = F)
newtag=newtag[newtag$SeriesID!="GSE133749",]
rownames(newtag)=paste(newtag$gene,newtag$SeriesID,sep="_")
head(newtag)
newtag$pert=ifelse(newtag$tag=="over",1,-1)
table(newtag$pert)
newtag.row=nrow(newtag)

#read counts
paste0(newtag$SeriesID,"_counts.txt") %in% list.files("input/counts/")
counts.list=list()
for(rn in rownames(newtag)){
  print(rn)
  SeriesID=newtag[rn,]$SeriesID
  print(SeriesID)
  counts.list[[rn]]=read.table(paste0("input/counts/",SeriesID,"_counts.txt"),
                               header = T,row.names = 1,stringsAsFactors = F)
}
length(counts.list)
str(counts.list[1])

#caculate cpm for each dataset
cpms.list=lapply(counts.list,cpm_normalized)
length(cpms.list)
head(cpms.list[[1]])
#zscores.list=lapply(cpms.list,cpm2zscore)
#zscores.list2=lapply(cpms.list,cpm2zscore_v1)
zscore.list=lapply(cpms.list,cpm2zscore_v2)
length(zscore.list)
str(zscore.list[1])

save(counts.list,cpms.list,zscore.list,file = "output/counts_cpms_zscore.rda")



