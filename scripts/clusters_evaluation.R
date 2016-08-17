
args <- commandArgs(TRUE)

tsvFile = args[1]

tab = read.table(file=tsvFile,header=TRUE,sep="\t",dec=".")

rownames(tab) = tab[,1]

path = dirname(tsvFile)
tab = tab[which(grepl(pattern = "cluster",x = tab[,1])),]
name = sub(pattern="(.*)_cluster_[0-9]*.*",replacement="\\1",x=tab[,1])[1]
tab[,1]= sub(pattern=".*(cluster_[0-9]*).*",replacement="\\1",x=tab[,1])
tab[,1]=gsub(pattern="[^0-9]",replacement="",x=tab[,1])


group = vector(length=30)
group = sapply(tab$Cluster,FUN=function(x){return(nchar(x))})
tab = cbind(tab,group)

group

#postscript(file = "rosto_cluster.eps",width =6 ,height =5 ,colormodel = "rgb",horizontal = FALSE,bg = "white")
bitmap(file = paste(path,"/",name,"_clusters_evaluation",".jpeg",sep=""),type="jpeg",res = 400 )
#par(mar=c(4,4,1,1))
plot(rep(0,length(unique(group)))+0.8,xaxt="n",lwd=2,xlab="Level of clustering",col=0,ylab="WR",type="l",ylim=c(0,1))#c(min(tab[,n]),max(tab[,n])))
axis(side = 1,at = 1:7,labels = c(2,4,8,16,32,64,128))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray",)
abline(h = seq(0,to=1,by=0.2),col="white",lwd=0.5, v = 1:length(unique(group)))

rec_cluster_tree = function(data, cl, index){
  size = sum(data$Nb_hits)
  
  limite = max(nchar(data$Cluster))
  k = which(colnames(data)==index)
  if(data$group[cl]==limite){
    return()
  }
  
  orig = c(data$group[cl],data[cl,k])

  #if(data$Nb_contigs[i] + data$Nb_contigs[j] == somme || data$Nb_contigs[i] + data$Nb_contigs[j] == somme - 1 ){
  
  i = which(data$Cluster==paste(data$Cluster[cl],"0",sep=""))
  j = which(data$Cluster==paste(data$Cluster[cl],"1",sep=""))
  
  lines(rbind(orig,c(data$group[i],data[i,k])),type="b",pch=16)#,cex = (data$Nb_Nemata[i]/data$Nb_hits[i])*5)
  lines(rbind(orig,c(data$group[j],data[j,k])),type="b",pch=16)#,cex = (data$Nb_Nemata[j]/data$Nb_hits[j])*5)
  rec_cluster_tree(data,i,index)
  rec_cluster_tree(data,j,index)
  return()
}

rec_cluster_tree(tab,which(tab$Cluster=="0"),index = "WR")
rec_cluster_tree(tab,which(tab$Cluster=="1"), index = "WR")

dev.off()

ssr = rep(0,times = length(unique(tab$group)[order(unique(tab$group))]))


for( i in unique(tab$group)[order(unique(tab$group))]){
  if(i!=1){
    bet = kmeans(tab$WR[which(tab$group == i)],centers = 2)[[6]]
    tot = kmeans(tab$WR[which(tab$group == i)],centers = 2)[[3]]
    ssr[i] = bet/tot
  }
}


mx = max.col(t(ssr))

km = kmeans(tab$WR[which(tab$group == mx)],centers=2)

cl_max = max.col(t(km$centers))
cl = km$cluster

data.frame(rownames(tab)[which(tab$group == mx)][which(cl==cl_max)])

