Biodv=function(file_name,datos,nall,distk,mayorque,menorque,missval,typedata,ht1,ht2,ht3){

#checkpack <- "pbkrtest"%in%rownames(installed.packages())
#pack<-paste(dirApp,"/local/pbkrtest_0.4-6.zip",sep="")
#userdir=.libPaths()[1]
#if(!checkpack) suppressWarnings(install.packages(pack,lib=userdir))

suppressWarnings(library(Hmisc))
suppressWarnings(library(plotly))
resulist=list()
#######################################################
id=as.data.frame(datos[,1])
colnames(id)=c("Markers")
datos=datos[,2:ncol(datos)]
#Add for markers and genotypes replicated
bymark=as.data.frame(apply(datos,1,paste,collapse = "" ))
bygen=as.data.frame(apply(datos,2,paste,collapse = "" ))

bymark[,1]=as.character(bymark[,1])
bygen[,1]=as.character(bygen[,1])

dupm = duplicated(bymark) | duplicated(bymark, fromLast = T)
dupm=split(which(dupm), bymark[dupm,])
dupm=lapply(dupm,function(x)id[x,])
if (length(dupm)!=0){
	max.length <- max(sapply(dupm, length))
	dupm <- lapply(dupm, function(v) { c(v, rep(NA, max.length-length(v)))})
	reptmark=data.frame(do.call(rbind, dupm))
	resulist[[1]]=reptmark
	#write.csv(reptmark,paste("MarkersRep_",file_name,".csv",sep=""),row.names=F)
}else {resulist[[1]]=NULL}

idG=as.data.frame(colnames(datos))
colnames(idG)=c("Genotypes")
#save(bymark,id,idG, bygen,file="chech.RData")

dupg = duplicated(bygen) | duplicated(bygen, fromLast = T)
dupg=split(which(dupg), bygen[dupg,])
dupg=lapply(dupg,function(x)idG[x,])
if (length(dupg)!=0){
	max.lengthg <- max(sapply(dupg, length))
	dupg <- lapply(dupg, function(v) { c(v, rep(NA, max.lengthg-length(v)))})
	reptgen=data.frame(do.call(rbind, dupg))
	resulist[[2]]=reptgen
	#write.csv(reptgen,paste("GenotypesRep_",file_name,".csv",sep=""),row.names=F)
}else {resulist[[2]]=NULL}

#######################################################
## do the filters for missing values in genotypes
if (missval!=0){
missgen=apply(datos,2,function(y) 2*sum(is.na(y))/(nall*nrow(datos)))
datos=datos[,which(missgen<missval)]
}
nacc=ncol(datos)
nalle=nall*nrow(datos)
nmark=nrow(datos)
#######################################################
datos$pmiss=apply(datos,1,function(y) sum(is.na(y))/nacc)
datos$pmiss1=apply((1-datos),1,function(y) sum(is.na(y))/nacc)
datos$pest=apply(datos[,1:nacc],1,mean,na.rm=T)
datos$pest1=apply((1-datos[,1:nacc]),1,mean,na.rm=T)
id2=data.frame(id,datos$pmiss,datos$pmiss1,datos$pest,datos$pest1)
colnames(id2)=c(colnames(id), "pmiss", "pmiss1" ,"pest","pest1")
##do the filter for polymorphism
if(mayorque!=0 & menorque!=0){
exaid=id2[which(id2$pest1<=mayorque & id2$pest<=mayorque & id2$pest>=menorque & id2$pest1>=menorque),]
datos=datos[which(datos$pest<=mayorque & datos$pest1<=mayorque & datos$pest>=menorque & datos$pest1>=menorque),]
}else{
exaid=id2
datos=datos
}
x2=as.data.frame(t(datos[,1:(ncol(datos)-4)]))
idg=as.data.frame(rownames(x2))
colnames(idg)=c("Genotypes")
## with this matrix calculate the genotype he,Ae,ho,etc
tmpmono=which(id2$pest1<=mayorque & id2$pest<=mayorque & id2$pest>=menorque & id2$pest1>=menorque)
if(length(tmpmono)!=nmark) resulist[[3]]=id2[-tmpmono,]#write.csv(id2[-tmpmono,],paste("MarkOutFilterPoly_",file_name,".csv",sep=""),row.names=F)
######################################################
rm(id2)
rm(tmpmono)
## the ID information useful
id1=exaid[,1]
######################################################
#Calculate specificity
cond1=which(datos$pest==0||is.na(datos$pest)==TRUE)
cond2=which(datos$pest1==0||is.na(datos$pest1)==TRUE)
longt=apply(datos[,1:nacc],1, length)-apply(datos[,1:nacc],1, function(y) length(which(is.na(y)==TRUE)))
adif01<--(-1)*(apply((datos[,1:nacc]/datos$pest),1, function (y) sum(y*log(y,2),na.rm=T))/longt)
adif02<--(-1)*(apply(((1-datos[,1:nacc])/datos$pest1),1, function (y) sum(y*log(y,2),na.rm=T))/longt)
if (length(cond1)!=0) adif01[cond1]=0
if (length(cond2)!=0) adif02[cond2]=0
##Calculate rareness
rar1=apply(datos[,1:nacc]*adif01,2,sum,na.rm=T)
rar2=apply((1-datos[,1:nacc])*adif02,2,sum,na.rm=T)
rareness<-(rar1+rar2)/nacc
rm(cond1,cond2,longt,rar1,rar2)
######################################################
######################################################
##Begins calculating indexes
######################################################
## observed heterozigosity
nhom=apply(datos[,1:nacc], 1, function(y) sum(y==1 | y==0, na.rm=T))

## expected heterozigosity, effective alleles, shannon index, Wright statistics
## prepare new columns for calculate the index PER LOCUS
exaid$pest2=exaid$pest^2
exaid$pest12=exaid$pest1^2
exaid$lpest=log2(exaid$pest)
exaid$lpest1=log2(exaid$pest1)
exaid$shan=-exaid$pest*exaid$lpest
exaid$shan1=-exaid$pest1*exaid$lpest1
exaid$ho=1-nhom/(nacc*(1-(exaid$pmiss)))

## calculus of He, Ho, Ae, shannon index, PER LOCUS
he=1-(exaid$pest2+exaid$pest12)
ae=1/(1-he)
ho=exaid$ho
shannon=exaid$shan+exaid$shan1
inb=he-ho

## To observed distributions indexes (each line correspond and create a plot to each index)
noNA=round(exaid$pmiss,4)
refmark=exaid[,1]
exadiv=data.frame(refmark,he,ho,ae,shannon,adif01,adif02,noNA)
colnames(exadiv)=c("Marker","He","Ho","Ae","Shannon","SpeAllele1","SpeAllele2","%NA")
resulist[[4]]=exadiv#write.csv(exadiv,"CalculusPerLocus.csv",row.names=FALSE,quote=FALSE)
rm(exaid)
#############################################################
He=mean(he)                                										## Expected Heterozygosity(diversidad genetica intrapoblacional)
se_he=sqrt(var(he)/length(he))             										## Standar desviation for Expected Heterozygosity
Ho=mean(ho)                                										## Observed Heterozygosity  (heterocigocidad promedio)
se_ho=sqrt(var(ho)/length(ho))             										## Standar desviation Observed Heterozygosity
Ae=mean(ae)                                										## Number of effective allele (diversidad genetica intrapoblacional)
se_ae=sqrt(var(ae)/length(ae))             										## Standar desviation for Number of effective allele
nm=length(unique(exadiv[,1]))
Shan=sum(shannon)/nm                       										## Shannon diversity index
se_shan=sqrt(var(shannon)/length(shannon)) 										## Standar desviation for Shannon diversity index
Inb=mean(inb)
se_inb=sqrt(var(inb)/length(inb))
PPoli=nm/nmark                                                ## Percentage of polymorphic loci
div=t(t(c(PPoli,He, se_he, Ho, se_ho, Ae, se_ae, Shan, se_shan)))
div=cbind(c("% of polymorphic loci","Expected Heterozygosity","Standar desviation for HE","Observed Heterozygosity", "Standar desviation for HO",
             "Number of effective allele", "Standar desviation for Ae", "Shannon diversity Index",
             "Standar desviation for ShanIn"),div)

rm(Ae,ae,He,he,ho,Ho,Inb,inb,nhom,PPoli,refmark,se_ae,se_he,se_ho,se_inb,se_shan,Shan,shannon)

nacc1=ncol(x2)
nmark1=nrow(x2)
x2$pmiss=apply(x2,1,function(y) sum(is.na(y))/nacc1)
x2$pmiss1=apply((1-x2),1,function(y) sum(is.na(y))/nacc1)
x2$pest=apply(x2[,1:nacc1],1,mean,na.rm=T)
x2$pest1=apply((1-x2[,1:nacc1]),1,mean,na.rm=T)
id2g=data.frame(idg,x2$pmiss,x2$pmiss1,x2$pest,x2$pest1)
colnames(id2g)=c(colnames(idg), "pmiss", "pmiss1" ,"pest","pest1")
exaidg=id2g
rm(id2g)
## the ID information useful
id1g=exaidg[,1]
## observed heterozigosity
nhomg=apply(x2[,1:nacc1], 1, function(y) sum(y==1 | y==0, na.rm=T))

## expected heterozigosity, effective alleles, shannon index, Wright statistics
## prepare new columns for calculate the index PER LOCUS
exaidg$pest2=exaidg$pest^2
exaidg$pest12=exaidg$pest1^2
exaidg$lpest=log2(exaidg$pest)
exaidg$lpest1=log2(exaidg$pest1)
exaidg$shan=-exaidg$pest*exaidg$lpest
exaidg$shan1=-exaidg$pest1*exaidg$lpest1
exaidg$ho=1-nhomg/(nacc1*(1-(exaidg$pmiss)))

heg=1-(exaidg$pest2+exaidg$pest12)
aeg=1/(1-heg)
hog=exaidg$ho
shannong=exaidg$shan+exaidg$shan1
## To observed distributions indexes (each line correspond and create a plot to each index)
noNA=round(x2$pmiss,4)
refGenotype=exaidg[,1]
exadivg=data.frame(refGenotype,heg,hog,aeg,shannong,rareness,noNA)
colnames(exadivg)=c("Genotype","He","Ho","Ae","Shannon","rareness","%NA")
resulist[[5]]=exadivg#write.csv(exadivg,"CalculusPerGenotype.csv",row.names=FALSE,quote=FALSE)
rm(nacc1,nmark1,x2,exaidg,id1g,nhomg,refGenotype,heg,aeg,hog,shannong,noNA,rareness)

library(car)
#out="SummaryDiversityAnalysis.csv"
if("div"%in%ls()==TRUE){
#cat("Diversity","\n","\n",file=out)
#write.table(div, file = out, append = T,quote=F, sep=",",col.names=F,row.names=F)
  resulist[[6]]=div
}

######################################################
##Calculate distances and cluster analysis
######################################################
fr=as.matrix(datos[,1:(ncol(datos)-4)])												## recover the marker information
frn=fr
frn[!is.na(frn)]=1																## no missing values convert to 1
frn[is.na(frn)]=0																  ## missing values convert to 0
N=2*crossprod(frn)																## create square matrix markers information
rm(frn)

if (distk=="Rogers"){
aux=matrix(0,nacc,nacc)
for(i in 1:(nacc-2)){
    aux[i,]=t(t(c(rep(0,i),sqrt((2*apply((fr[,i]-fr[,-c(1:i)])^2,2,function(x) sum(x,na.rm=T)))/N[i,-c(1:i)]))))
}
aux[nacc-1,]=t(t(c(rep(0,nacc-1),sqrt((2*sum((fr[,nacc-1]-fr[,-c(1:nacc-1)])^2,na.rm=T))/N[nacc-1,-c(1:(nacc-1))]))))
mrdMAT=aux+t(aux)
}

if (distk=="Nei"){
  fr[is.na(fr)]=0
  aux <- t(fr)%*%fr
  vec <- sqrt(diag(aux))
  aux <- aux/vec[col(aux)]
  aux <- aux/vec[row(aux)]
  aux <- -log(aux)
  mrdMAT=aux
}

rm(fr)
colnames(mrdMAT)=colnames(N)
rownames(mrdMAT)=rownames(N)
mrdMAT1=mrdMAT
colnames(mrdMAT1)=seq(1:dim(mrdMAT)[1])
resulist[[7]]=cbind(ID=seq(1:dim(mrdMAT)[1]),NAME=rownames(mrdMAT),round(mrdMAT1,5))#write.csv(cbind(ID=seq(1:dim(mrdMAT)[1]),NAME=rownames(mrdMAT),round(mrdMAT1,5)),paste(distk,"Distances.csv",sep=""),quote=FALSE,row.names=FALSE)
#write.table(cbind(NAME=rownames(mrdMAT),round(mrdMAT,5)),"DistancesforGAP.txt",sep="\t",quote=FALSE,row.names=FALSE)
rm(mrdMAT1)
##graphical representation, the MDS (multidimensional scaling) analysis
mds=cmdscale(mrdMAT, k=3,eig=T)
coord=mds$points
perctCP12=c(round(mds$eig[1]/sum(mds$eig)*100,2),round(mds$eig[2]/sum(mds$eig)*100,2),round(mds$eig[3]/sum(mds$eig)*100,2))
colnames(coord)=c("dim1","dim2","dim3")
rm(aux,N,mds)
#########################################################################
#########################################################################
###For do dendograms and MDSgraph
#########################################################################
#########################################################################
clust=cluster::agnes(mrdMAT, method = "ward")

library(dendextend)
coord1=as.data.frame(coord)
names(coord1)=c("Factor1","Factor2","Factor3")
resulist[[8]]=coord1#write.csv(coord1,"MDStable.csv",quote=FALSE)

names(resulist)=c("Repeated_Markers","Repeated_Genotypes","MarkOutFilterPoly","CalculusPerMarker","CalculusPerGenotype","SummaryDiversityAnalysis","Distance_Matrix","MDSTable")

coord2=cbind(gen=rownames(coord),coord)
names(coord2)=c("Gen","Factor1","Factor2","Factor3")
res=list(as.data.frame(div),coord2, getwd(), clust, datos, mrdMAT, perctCP12, resulist)
return(res)
}
