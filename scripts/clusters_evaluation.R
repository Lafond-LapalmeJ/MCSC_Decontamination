
args <- commandArgs(TRUE)

tsvFile = args[1]

tab = read.table(file=tsvFile,header=TRUE,sep="\t",dec=".")

rownames(tab) = tab[,1]

path = dirname(tsvFile)

name = sub(pattern="(.*)_cluster_[0-9]*.*",replacement="\\1",x=tab[,1])[1]
tab[,1]= sub(pattern=".*(cluster_[0-9]*).*",replacement="\\1",x=tab[,1])
tab[,1]=gsub(pattern="[^0-9]",replacement="",x=tab[,1])


group = vector(length=30)
group = sapply(tab$Cluster,FUN=function(x){return(nchar(x))})
tab = cbind(tab,group)

#postscript(file = "rosto_cluster.eps",width =6 ,height =5 ,colormodel = "rgb",horizontal = FALSE,bg = "white")
bitmap(file = paste(path,"/",name,"_clusters_evaluation",".jpeg",sep=""),type="jpeg",res = 400 )
plot(rep(0,length(unique(group)))+0.8,xaxt="n",lwd=2,xlab="Level of clustering",col=0,ylab="WR",type="l",ylim=c(0,max(tab$WR)))#c(min(tab[,n]),max(tab[,n])))
axis(side = 1,at = 1:5,labels = c(2,4,8,16,32))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
abline(h = seq(0,to=1,by=0.2),col=0,lwd=0.5, v = 1:length(unique(group)))

rec_cluster_tree = function(data, cl, index){
  limite = max(nchar(data$Cluster))
  k = which(colnames(data)==index)
  if(data$group[cl]==limite){
    return()
  }
  
  orig = c(data$group[cl],data[cl,k])

  #if(data$Nb_contigs[i] + data$Nb_contigs[j] == somme || data$Nb_contigs[i] + data$Nb_contigs[j] == somme - 1 ){
  
  i = which(data$Cluster==paste(data$Cluster[cl],"0",sep=""))
  j = which(data$Cluster==paste(data$Cluster[cl],"1",sep=""))
  
  lines(rbind(orig,c(data$group[i],data[i,k])),type="b",pch=16)
  lines(rbind(orig,c(data$group[j],data[j,k])),type="b",pch=16)
  rec_cluster_tree(data,i,index)
  rec_cluster_tree(data,j,index)
  return()
}

rec_cluster_tree(tab,which(tab$Cluster=="0"),index = "WR")
rec_cluster_tree(tab,which(tab$Cluster=="1"), index = "WR")

dev.off()
