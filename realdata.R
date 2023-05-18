source("localancestryfns.R")
### Start of real data processing

## This script requires that the following files have been created by convertmosaic/convert_mosaic.R


samples=read.csv("convertmosaic/chr16data/samples.csv",row.names=1)
localancestry=as.matrix(read.csv("convertmosaic/chr16data/localancestry.csv",row.names=1))
gtdata=as.matrix(read.csv("convertmosaic/chr16data/gtdata.csv",row.names=1))
rates=as.matrix(read.table("convertmosaic/chr16data/rates.16",skip=1))

## We are going to have to simulate a rapid growth of the population to get enough individuals for a breeding program.
## The following procedure:
## a) goes from 24 to 100 individuals (by using all haplotypes and then randomly adding additional copies of them)
## b) downsamples to a specified number of SNPs. As long as this is large enough to
##    capture LD, it won't affect the mean dynamics but massively changes runtime

set.seed(2) # Due to the upsampling step, some seeds make a bad starting population
## Make a larger population
Ntarget=200 # Number of target haplotypes
Ltarget=5000 # Number of target SNPs
## We could breed based on allele sharing via ibs, but no need for this demo
## ibs=(gtdata %*% t(gtdata))/dim(gtdata)[2] 
myhaps=sample(c(1:dim(localancestry)[1],
                sample(1:dim(localancestry)[1],
                       Ntarget-dim(localancestry)[1],
                       replace=TRUE)))

## Construct the data we will use: 
gpos=sort(sample(1:dim(localancestry)[2],Ltarget)) # The random SNPs included in the analysis
data=list(panel=localancestry[myhaps,gpos],genomes=gtdata[myhaps,gpos]) # The upsamples panel of data
gdists=c(0,diff(as.numeric(rates[2,gpos]))) # The genetic distances

N=dim(data[[1]])[1]/2 # Individuals
L=dim(data[[1]])[2] # Number of SNPs

## Perform forward simulation
Gforward=100
verbose=FALSE
dataNull<-breedingSel(data,Gforward,selfn=selIdentity,d=gdists,verbose=verbose)
dataRandom<-breedingSel(data,Gforward,selfn=selRandom,d=gdists,verbose=verbose)
dataQ<-breedingSel(data,Gforward,selfn=selMinQ,d=gdists,verbose=verbose)
dataHet<-breedingSel(data,Gforward,selfn=selHet,d=gdists,verbose=verbose)
dataHetQ<-breedingSel(data,Gforward,selfn=selHetQ,d=gdists,verbose=verbose)

## Choices of scores to report
scores=c("Q","Het","HetQ")
names=c("Domestic Admixture","Heterozygosity","Heterozygosity for wildcat loci")
names(names)=scores
## Plot
pdf(paste0("SimulatedBreedingWildcat_L",L,".pdf"),height=6,width=12)
par(mfrow=c(1,3))
for(score in scores){
    yrange=c(0,0.35)
    if(score=="Q") yrange=c(0,0.5)
    if(score=="EffectiveLoci") yrange=c(0,1)
    plot(range(dataNull$score$g),yrange,xlab="Generation",ylab=score,type="n",main=paste("Evolution of",names[score]))
##    lines(dataNull$score[,"g"],dataNull$score[,score],col="grey")
    lines(dataRandom$score[,"g"],dataNull$score[,score],col=1,lty=1)
    lines(dataQ$score[,"g"],dataQ$score[,score],col=2,lty=1)
    lines(dataHet$score[,"g"],dataHet$score[,score],col=3,lty=1)
    lines(dataHetQ$score[,"g"],dataHetQ$score[,score],col=4,lty=1)
    legend("topleft",legend=c("Random","Widcat Genome Fraction","Heterozygosity","Wildcat Heterozygosity"),
           text.col=1:4,lty=1,col=c(1,2,3,4),title="Breeding top 50% by:")
}
dev.off()
