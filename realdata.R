source("localancestryfns.R")
### Start of real data processing


## This script requires that the following files have been created by convertmosaic/convert_mosaic.R

## From Appendix 1 of SW_PVA_2023_Reportv2 (002).docx:
# Mean brood size **per year**:
meanbrood=sum(c(1,2,3,4,5)*c(23,36,30,10,1)/100)
## Mortality in y 1 = 0.3
meanbrood*(1-0.3)
## 1.61
## Age at first offspring 1, max age 8
## F mortality in adults = 0.1, so total number of offspring per F = 
totalbrood=sum(meanbrood * (1-0.3) * ((1-0.1)^(0:6)))
totalbrood
## 8.4
## 50% adult F breeding
totalbrood*0.5
## 4.2 children on average per female per generation
##
## Also 150 is the size of the population

samples=read.csv("chr16data/samples.csv",row.names=1)
localancestry=1-as.matrix(read.csv("chr16data/localancestry.csv",row.names=1)) # reported as domestic proportions
gtdata=as.matrix(read.csv("chr16data/gtdata.csv",row.names=1))
rates=as.matrix(read.table("chr16data/rates.16",skip=1))

## We are going to have to simulate a rapid growth of the population to get enough individuals for a breeding program.
## The following procedure:
## a) goes from 24 to 150 individuals (by using all haplotypes and then randomly adding additional copies of them)
## b) downsamples to a specified number of SNPs. As long as this is large enough to
##    capture LD, it won't affect the mean dynamics but massively changes runtime

set.seed(2) # Due to the upsampling step, some seeds make a bad starting population
## Make a larger population
Nhaptarget=300 # Number of target haplotypes
Ltarget=5000 # Number of target SNPs
## We could breed based on allele sharing via ibs, but no need for this demo
## ibs=(gtdata %*% t(gtdata))/dim(gtdata)[2] 
myhaps=sample(c(1:dim(localancestry)[1],
                sample(1:dim(localancestry)[1],
                       Nhaptarget-dim(localancestry)[1],
                       replace=TRUE)))

## Construct the data we will use: 
gpos=sort(sample(1:dim(localancestry)[2],Ltarget)) # The random SNPs included in the analysis
wildcatdata=list(panel=localancestry[myhaps,gpos],
                 genomes=gtdata[myhaps,gpos],
                 founder=matrix(myhaps,nrow=Nhaptarget,ncol=Ltarget)) # The upsamples panel of data
gdists=c(0,diff(as.numeric(rates[2,gpos]))) # The genetic distances

N=dim(wildcatdata[[1]])[1]/2 # Individuals
L=dim(wildcatdata[[1]])[2] # Number of SNPs

## Perform forward simulation
Gforward=50
verbose=TRUE
ignore<-function(){
    wildcatsimpleNull<-breedingSel(wildcatdata,Gforward,
                                   selfn=selIdentity,d=gdists,verbose=verbose)
    wildcatsimpleRandom<-breedingSel(wildcatdata,Gforward,
                                     selfn=selRandom,d=gdists,verbose=verbose)
    wildcatsimpleQ<-breedingSel(wildcatdata,Gforward,
                                selfn=selMinQ,d=gdists,verbose=verbose)
    wildcatsimpleHet<-breedingSel(wildcatdata,Gforward,
                                  selfn=selHet,d=gdists,verbose=verbose)
    wildcatsimpleHetQ<-breedingSel(wildcatdata,Gforward,
                                   selfn=selHetQ,d=gdists,verbose=verbose)
}

meanoff=4.2
wildcatNull<-breedingSel(wildcatdata,Gforward,selfn=pairwiseRandom,
                         parentfn=rankedParents,offspringdist=function(N)rpois(N,meanoff))
wildcatPairwiseQ<-breedingSel(wildcatdata,Gforward,selfn=pairwiseQ,
                              parentfn=rankedParents,min=FALSE,
                              offspringdist=function(N)rpois(N,meanoff))
wildcatPairwiseMK<-breedingSel(wildcatdata,Gforward,selfn=pairwiseKinship,
                               parentfn=rankedParents,min=TRUE,
                              offspringdist=function(N)rpois(N,meanoff))
wildcatPairwiseMH<-breedingSel(wildcatdata,Gforward,selfn=pairwiseHeterozygosity,
                               parentfn=rankedParents,min=FALSE,
                              offspringdist=function(N)rpois(N,meanoff))
wildcatPairwiseMHQ<-breedingSel(wildcatdata,Gforward,selfn=pairwiseHeterozygosityQ,
                                parentfn=rankedParents,min=FALSE,
                              offspringdist=function(N)rpois(N,meanoff))
wildcatPairwiseMHQScore<-breedingSel(wildcatdata,Gforward,selfn=pairwiseHeterozygosityQScore,
                                     parentfn=rankedParents,min=FALSE,
                              offspringdist=function(N)rpois(N,meanoff))
wildcatPairwiseMKQ<-breedingSel(wildcatdata,Gforward,selfn=pairwiseKinshipQ,
                                parentfn=rankedParents,min=TRUE,
                              offspringdist=function(N)rpois(N,meanoff))
wildcatPairwiseMKQScore<-breedingSel(wildcatdata,Gforward,selfn=pairwiseKinshipQScore,
                                     parentfn=rankedParents,min=TRUE,
                              offspringdist=function(N)rpois(N,meanoff))

save.image(paste0("SimulatedBreedingWildcat_L",L,".RData"))

## Choices of scores to report
scores=c("NonQ","Het","HetQ")
names=c("Domestic Admixture","Heterozygosity","Heterozygosity for wildcat loci")
names(names)=scores
## Plot
pdf(paste0("SimulatedBreedingWildcat_L",L,".pdf"),height=6,width=12)
par(mfrow=c(1,3))
for(score in scores){
    yrange=c(0,1)
    if(score=="Q") yrange=c(0,0.5)
    if(score=="EffectiveLoci") yrange=c(0,1)
    plot(range(wildcatNull$score$g),yrange,xlab="Generation",ylab=score,type="n",main=paste("Evolution of",names[score]))
##    lines(wildcatdataNull$score[,"g"],wildcatdataNull$score[,score],col="grey")
    lines(wildcatNull$score[,"g"],wildcatNull$score[,score],col=1,lty=1)
    lines(wildcatPairwiseQ$score[,"g"],wildcatPairwiseQ$score[,score],col=2,lty=1)
    lines(wildcatPairwiseMH$score[,"g"],wildcatPairwiseMH$score[,score],col=3,lty=1)
    lines(wildcatPairwiseMHQ$score[,"g"],wildcatPairwiseMHQ$score[,score],col=4,lty=1)
    legend("topleft",legend=c("Random","Widcat Genome Fraction","Heterozygosity","Wildcat Heterozygosity"),
           text.col=1:4,lty=1,col=c(1,2,3,4),title="Breeding top 50% by:")
}
dev.off()


## Figure 3
scores=c("NonQ","Het","Kin","HetQ","KinQ")
names=c("Domestic Admixture","Heterozygosity","Kinship","Wildcat Heterozygosity","Wildcat Kinship")
symbols=c("1-Q","H","K","SH","SK")
names(symbols)=names(names)=scores
pdf("Figure5-Wildcat.pdf",height=8,width=10)
layout(matrix(c(1,2,3,6,4,5),nrow=2,byrow=TRUE))
par(las=1)
for(scoreon in 1:length(scores)){
    score=scores[scoreon]
    panel=letters[scoreon]
    yrange=c(0,0.5)
    if(score=="Kin") yrange=c(0,1)
    if(score=="KinQ") yrange=c(0,1)
    plot(range(wildcatNull$score$g),yrange,xlab="Generation",ylab=symbols[score],
         type="n",main="",cex.axis=1.5,cex.lab=1.2)
    mtext(paste0(panel,") Evolution of ",names[score]),adj=-0.2,line=1,cex=0.8)
    lines(wildcatNull$score[,"g"],wildcatNull$score[,score],col=1,lty=3,lwd=2)
    lines(wildcatPairwiseQ$score[,"g"],wildcatPairwiseQ$score[,score],col=2,lty=3,lwd=2)
    lines(wildcatPairwiseMH$score[,"g"],wildcatPairwiseMH$score[,score],col=3,lty=3,lwd=2)
    lines(wildcatPairwiseMK$score[,"g"],wildcatPairwiseMK$score[,score],col=5,lty=3,lwd=2)
    lines(wildcatPairwiseMHQ$score[,"g"],wildcatPairwiseMHQ$score[,score],col=3,lty=1,lwd=2)
    lines(wildcatPairwiseMKQ$score[,"g"],wildcatPairwiseMKQ$score[,score],col=5,lty=1,lwd=2)
    lines(wildcatPairwiseMHQScore$score[,"g"],wildcatPairwiseMHQScore$score[,score],col=4,lwd=2)
    lines(wildcatPairwiseMKQScore$score[,"g"],wildcatPairwiseMKQScore$score[,score],col=6,lwd=2)
}
par(mar=c(0,0,0,0))
plot(c(0,1),c(0,1), type="n",xlab="",ylab="",axes=FALSE)
legend("center",legend=c("Random","Min Introgression","Max Heterozygosity","Min Kinship",
                          "Max Wildcat Heterozygosity SH",
                                              "Min Wildcat Kinship SK",
                                              "Max Weighted Heterozygosity WSH",
                                              "Min Weighted Kinship WSK"),
           text.col=c(1,2,3,5,3,5,4,6),lty=c(3,3,3,3,1,1,1,1),col=c(1,2,3,5,3,5,4,6),title="Breeding ranking:",cex=1.2,lwd=2)
dev.off()
