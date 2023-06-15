source("localancestryfns.R")
### Start of sim

## This is a stand-alone script that should work regardless of data sources
set.seed(1)
N=150 # Individuals
NI=50 # Introgressing individuals
L=500 # Amount of genome, in SNPs. We assume cM for a basic sim
Fst=0.2 # Fst between MRCA of two populations
f0=runif(L,0.05,0.95) # Frequencies in the MRCA
f=simf(f0,Fst) # Frequencies in the target population
fintrogress=simf(f0,Fst) # Frequencies in the introgressing population
G=10 # Number of gens after introgression starts
infracseq=c(rep(0.3,2),rep(0,G-2)) # Introgression sequence
genomes=sapply(f,function(p)rbinom(2*N,1,p))
genomesIntrogress=sapply(fintrogress,function(p)rbinom(2*NI,1,p))
founder=matrix(rep(1:N,each=2),N*2,L)
data=list(panel=matrix(1,N*2,L),genomes=genomes,founder=founder)
attr(data,"sexf")=rbinom(N,1,0.5)
dataIntrogress=list(panel=matrix(0,NI*2,L),
                    genomes=genomesIntrogress,
                    founder=matrix(rep(1:NI,each=2),NI*2,L))
attr(dataIntrogress,"sexf")=rbinom(NI,1,0.5)

## Generate the initial panel
data0=breedingIntrogress(data,dataIntrogress,G,infracseq)
data0$founder=matrix(rep(1:N,each=2),N*2,L) ## Re-founder the panel
## Generate the captive population
data<-breedingSel(data0,1,selfn=selRandom,frac=0.25)$data ## Random mating bottleneck
data<-breedingSel(data,1,selfn=selMaxQ,frac=0.25)$data ## Removal of highly introgressed individuals

## Plot of what the distribution of introgresion looks like 
par(mfrow=c(2,2))
hist(L-rowSums(data0[[1]]),breaks=seq(0,L,length=21),
     xlab="Introgressed sites",ylab="haplotype count",main="Number of introgressed sites per haplotype (before captive program)")
hist(2*N-colSums(data0[[1]]),breaks=seq(0,2*N,length=21),
     xlab="Introgressed haplotypes",ylab="SNP count",main="Number of introgressed individuals per site (before captive program)")
hist(L-rowSums(data[[1]]),breaks=seq(0,L,length=21),
     xlab="Introgressed sites",ylab="haplotype count",main="Number of introgressed sites per haplotype (after captive program)")
hist(2*N-colSums(data[[1]]),breaks=seq(0,2*N,length=21),
     xlab="Introgressed haplotypes",ylab="SNP count",main="Number of introgressed individuals per site (after captive program)")
allScores(0,data0)
allScores(0,data)

mean(data0[[1]])
par(mfrow=c(2,2))
image(t(data0[[1]]))
image(t(data0[[2]]))
image(1-t(data[[1]]))
image(t(data[[2]]))

### Breeding step


Gforward=40
dataNull<-breedingSel(data,Gforward,selfn=selRandom,frac=1)
dataRandom<-breedingSel(data,Gforward,selfn=selRandom)
dataQ<-breedingSel(data,Gforward,selfn=selMaxQ)
dataHet<-breedingSel(data,Gforward,selfn=selHet)
## dataHetQ<-breedingSel(data,Gforward,selfn=selHetQ) ## Terrible raw score
dataHetQwt<-breedingSel(data,Gforward,selfn=selHetQwt) 
dataKin<-breedingSel(data,Gforward,selfn=selKin)
dataKinQ<-breedingSel(data,Gforward,selfn=selKinQ)
dataExpectedKin<-breedingSel(data,Gforward,selfn=selExpectedKin)
dataExpectedKinQ<-breedingSel(data,Gforward,selfn=selExpectedKinQ)
## dataKinQwt<-breedingSel(data,Gforward,selfn=selKinQwt) ## Terrible weighted score

dataPairwiseRandom<-breedingSel(data,Gforward,selfn=pairwiseRandom,parentfn=rankedParents)
dataPairwiseQ<-breedingSel(data,Gforward,selfn=pairwiseQ,parentfn=rankedParents,min=FALSE)
dataPairwiseMK<-breedingSel(data,Gforward,selfn=pairwiseKinship,parentfn=rankedParents,min=TRUE)
dataPairwiseMH<-breedingSel(data,Gforward,selfn=pairwiseHeterozygosity,parentfn=rankedParents,min=FALSE)
dataPairwiseMHQ<-breedingSel(data,Gforward,selfn=pairwiseHeterozygosityQ,parentfn=rankedParents,min=FALSE)
dataPairwiseMHQScore<-breedingSel(data,Gforward,selfn=pairwiseHeterozygosityQScore,parentfn=rankedParents,min=FALSE)
dataPairwiseMKQ<-breedingSel(data,Gforward,selfn=pairwiseKinshipQ,parentfn=rankedParents,min=TRUE)
dataPairwiseMKQScore<-breedingSel(data,Gforward,selfn=pairwiseKinshipQScore,parentfn=rankedParents,min=TRUE)

scores=c("NonQ","Het","HetQ","K")
names=c("Domestic Admixture","Heterozygosity","Heterozygosity for target loci","Kinship")
names(names)=scores

save.image("sim.RData")

## pdf("SimulatedBreeding.pdf",height=6,width=12)
par(mfrow=c(2,2))
for(score in scores){
    yrange=c(0,0.35)
    if(score=="EffectiveLoci") yrange=c(0,1)
    plot(range(dataNull$score$g),yrange,xlab="Generation",ylab=score,type="n",main=paste("Evolution of",names[score]))
##    lines(dataNull$score[,"g"],dataNull$score[,score],col="grey")
    lines(dataRandom$score[,"g"],dataNull$score[,score],col=1,lty=1)
    lines(dataQ$score[,"g"],dataQ$score[,score],col=2,lty=1)
    lines(dataHet$score[,"g"],dataHet$score[,score],col=3,lty=1)
    lines(dataHetQwt$score[,"g"],dataHetQwt$score[,score],col=4,lty=1)
    legend("topleft",legend=c("Random","Widcat Genome Fraction","Heterozygosity","Target Heterozygosity"),
           text.col=1:4,lty=1,col=c(1,2,3,4),title="Breeding top 50% by:")
}
##dev.off()

scores=c("NonQ","Het","HetQ","Kin","KinQ")
names=c("Domestic Admixture","Heterozygosity","Heterozygosity for Target loci","Kinship","Kinship for Target loci")
#scores=c("NonQ","Het","HetQ","Kin","KinQ", "ExpectedKinQ")
#names=c("Domestic Admixture","Heterozygosity","Heterozygosity for wildcat loci","Kinship","Kinship for wildcat loci","Expected Kinship for wildcat loci")
names(names)=scores

############
## All sampling results
#pdf("SimulatedBreeding2.pdf",height=6,width=12)
par(mfrow=c(2,3))
for(score in scores){
    yrange=c(0,0.35)
    if(score=="EffectiveLoci") yrange=c(0,1)
    if(score=="Kin") yrange=c(0,0.5)
    if(score=="ExpectedKin") yrange=c(0,0.5)
    if(score=="ExpectedKinQ") yrange=c(0,1)
    if(score=="KinQ") yrange=c(0,1)
    plot(range(dataNull$score$g),yrange,xlab="Generation",ylab=score,type="n",main=paste("Evolution of",names[score]))
##    lines(dataNull$score[,"g"],dataNull$score[,score],col="grey")
    lines(dataRandom$score[,"g"],dataNull$score[,score],col=1,lty=1)
    lines(dataQ$score[,"g"],dataQ$score[,score],col=2,lty=1)
    lines(dataHet$score[,"g"],dataHet$score[,score],col=3,lty=1)
#    lines(dataHetQ$score[,"g"],dataHetQ$score[,score],col=4,lty=1)
    lines(dataHetQwt$score[,"g"],dataHetQwt$score[,score],col=4,lty=2)
    lines(dataKin$score[,"g"],dataKin$score[,score],col=5,lty=1)
    lines(dataKinQ$score[,"g"],dataKinQ$score[,score],col=6,lty=1)
    lines(dataExpectedKin$score[,"g"],dataExpectedKin$score[,score],col=5,lty=2,lwd=2)
    lines(dataExpectedKinQ$score[,"g"],dataExpectedKinQ$score[,score],col=6,lty=2,lwd=2)
##    lines(dataKinQwt$score[,"g"],dataKinQwt$score[,score],col=6,lty=2)
    legend("topleft",legend=c("Random","Widcat Genome Fraction","Heterozygosity",
                              "Target Heterozygosity","Kinship","Target Kinship"),
           text.col=1:6,lty=1,col=1:6,title="Breeding top 50% by:",cex=0.8)
}
#dev.off()

############
## All pairwise results
#pdf("SimulatedBreeding2.pdf",height=6,width=12)
par(mfrow=c(2,3))
for(score in scores){
    yrange=c(0,0.4)
    if(score=="EffectiveLoci") yrange=c(0,1)
    if(score=="Kin") yrange=c(0,0.5)
    if(score=="ExpectedKin") yrange=c(0,1)
    if(score=="ExpectedKinQ") yrange=c(0,1)
    if(score=="KinQ") yrange=c(0,1)
    plot(range(dataNull$score$g),yrange,xlab="Generation",ylab=score,type="n",main=paste("Evolution of",names[score]))
##    lines(dataNull$score[,"g"],dataNull$score[,score],col="grey")
##    lines(dataNull$score[,"g"],dataNull$score[,score],col=1,lty=1)
    lines(dataPairwiseRandom$score[,"g"],dataPairwiseRandom$score[,score],col=1)
    lines(dataPairwiseQ$score[,"g"],dataPairwiseQ$score[,score],col=2)
    lines(dataPairwiseMH$score[,"g"],dataPairwiseMH$score[,score],col=3)
    lines(dataPairwiseMK$score[,"g"],dataPairwiseMK$score[,score],col=5)
    lines(dataPairwiseMHQ$score[,"g"],dataPairwiseMHQ$score[,score],col=4,lty=2)
    lines(dataPairwiseMHQScore$score[,"g"],dataPairwiseMHQScore$score[,score],col=4)
    lines(dataPairwiseMKQ$score[,"g"],dataPairwiseMKQ$score[,score],col=6,lty=2)
    lines(dataPairwiseMKQScore$score[,"g"],dataPairwiseMKQScore$score[,score],col=6)
##    lines(dataKinQwt$score[,"g"],dataKinQwt$score[,score],col=6,lty=2)
    legend("topleft",legend=c("Random","Widcat Genome Fraction","Heterozygosity",
                              "Target Heterozygosity","Kinship","Target Kinship"),
           text.col=1:6,lty=1,col=1:6,title="Breeding top 50% by:",cex=0.8)
}
#dev.off()

dataNull$score
dataRandom$score
dataQ$score
dataHet$score
dataHetQ$score

image(t(dataHetQ$data[[1]]))

## A simulation based approach: not needed
## lastdata=data
## nrep=10
## scoremat<-sapply(1:N,function(i){
##     scores=matrix(NA,nrow=N,ncol=nrep)
##     for(rep in 1:nrep){
##         tp1=crossover(getinddata(lastdata,i))
##         for(j in 1:N){
##             tp2=crossover(getinddata(lastdata,j))
##             scores[j,rep]=scorefn(combineData(tp1,tp2))
##         }
##     }
##     rowMeans(scores)
## })

##############################
## Figure 2
scores=c("NonQ","Het","K")
names=c("Introgressed Admixture","Heterozygosity","Kinship")
symbols=c("                  1-Q      (More Introgressed)",
          "                   H       (More Diverse)",
          "                   K       (More Inbred)")
names(symbols)=names(names)=scores
pdf("Figure2-SimpleModels.pdf",height=4,width=10)
par(mfrow=c(1,3),las=1)
for(scoreon in 1:length(scores)){
    score=scores[scoreon]
    panel=letters[scoreon]
    yrange=c(0,0.5)
    if(score=="K") yrange=c(0,1)
    plot(range(dataNull$score$g),yrange,xlab="Generation",ylab=symbols[score],
         type="n",main="",cex.axis=1.5,cex.lab=1.2)
    mtext(paste0(panel,") Evolution of ",names[score]),adj=-0.2,line=1,cex=0.8)
    lines(dataPairwiseRandom$score[,"g"],dataPairwiseRandom$score[,score],col=1,lwd=2,cex=0.75)
    lines(dataPairwiseQ$score[,"g"],dataPairwiseQ$score[,score],col=2,lwd=2)
    lines(dataPairwiseMH$score[,"g"],dataPairwiseMH$score[,score],col=3,lwd=2)
    lines(dataPairwiseMK$score[,"g"],dataPairwiseMK$score[,score],col=5,lwd=2)
    if(score=="Het") legend("topleft",legend=c("Random R","Min Introgression Q","Max Heterozygosity H","Min Kinship K"),
           text.col=c(1,2,3,5),lty=1,col=c(1,2,3,5),title="Breeding ranking:",cex=1,lwd=2)
}
dev.off()

## Figure 3
scores=c("NonQ","Het","Kin","HetQ","KinQ")
names=c("Introgressed Admixture","Heterozygosity","Kinship","Subpop Heterozygosity","Subpop Kinship")
symbols=c("1-Q","H","K","SH","SK")
names(symbols)=names(names)=scores
pdf("Figure3-FullModels.pdf",height=8,width=10)
layout(matrix(c(1,2,3,6,4,5),nrow=2,byrow=TRUE))
par(las=1)
for(scoreon in 1:length(scores)){
    score=scores[scoreon]
    panel=letters[scoreon]
    yrange=c(0,0.5)
    if(score=="Kin") yrange=c(0,1)
    if(score=="KinQ") yrange=c(0,1)
    plot(range(dataNull$score$g),yrange,xlab="Generation",ylab=symbols[score],
         type="n",main="",cex.axis=1.5,cex.lab=1.2)
    mtext(paste0(panel,") Evolution of ",names[score]),adj=-0.2,line=1,cex=0.8)
    lines(dataPairwiseRandom$score[,"g"],dataPairwiseRandom$score[,score],col=1,lty=3,lwd=2)
    lines(dataPairwiseQ$score[,"g"],dataPairwiseQ$score[,score],col=2,lty=3,lwd=2)
    lines(dataPairwiseMH$score[,"g"],dataPairwiseMH$score[,score],col=3,lty=3,lwd=2)
    lines(dataPairwiseMK$score[,"g"],dataPairwiseMK$score[,score],col=5,lty=3,lwd=2)
    lines(dataPairwiseMHQ$score[,"g"],dataPairwiseMHQ$score[,score],col=3,lty=1,lwd=2)
    lines(dataPairwiseMKQ$score[,"g"],dataPairwiseMKQ$score[,score],col=5,lty=1,lwd=2)
    lines(dataPairwiseMHQScore$score[,"g"],dataPairwiseMHQScore$score[,score],col=4,lwd=2)
    lines(dataPairwiseMKQScore$score[,"g"],dataPairwiseMKQScore$score[,score],col=6,lwd=2)
}
par(mar=c(0,0,0,0))
plot(c(0,1),c(0,1), type="n",xlab="",ylab="",axes=FALSE)
legend("center",legend=c("Random","Min Introgression","Max Heterozygosity","Min Kinship",
                          "Max Subpop Heterozygosity SH",
                                              "Min Subpop Kinship SK",
                                              "Max Weighted Heterozygosity WSH",
                                              "Min Weighted Kinship WSK"),
           text.col=c(1,2,3,5,3,5,4,6),lty=c(3,3,3,3,1,1,1,1),col=c(1,2,3,5,3,5,4,6),title="Breeding ranking:",cex=1.2,lwd=2)
dev.off()
