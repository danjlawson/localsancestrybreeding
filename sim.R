source("localancestryfns.R")
### Start of sim

## This is a stand-alone script that should work regardless of data sources

N=100 # Individuals
NI=50 # Introgressing individuals
p=2 # ploidy
L=500 # Amount of genome, in SNPs. We assume cM for a basic sim
Fst=0.2 # Fst between MRCA of two populations
f0=runif(L,0.05,0.95) # Frequencies in the MRCA
f=simf(f0,Fst) # Frequencies in the target population
fintrogress=simf(f0,Fst) # Frequencies in the introgressing population
G=10 # Number of gens after introgression starts
infracseq=c(rep(0.3,2),rep(0,G-2)) # Introgression sequence
genomes=sapply(f,function(p)rbinom(2*N,1,p))
genomesIntrogress=sapply(fintrogress,function(p)rbinom(2*NI,1,p))
data=list(panel=matrix(0,N*2,L),genomes=genomes)
dataIntrogress=list(panel=matrix(1,NI*2,L),genomes=genomesIntrogress)

## Generate the initial panel
data0=breedingIntrogress(data,dataIntrogress,G,infracseq)
## Generate the captive population
data<-breedingSel(data0,1,selfn=selRandom,frac=0.25)$data ## Random mating bottleneck
data<-breedingSel(data,1,selfn=selMinQ,frac=0.25)$data ## Removal of highly introgressed individuals

## Plot of what the distribution of introgresion looks like 
par(mfrow=c(2,2))
hist(rowSums(data0[[1]]),breaks=seq(0,L,length=21),
     xlab="Introgressed sites",ylab="haplotype count",main="Number of introgressed sites per haplotype (before captive program)")
hist(colSums(data0[[1]]),breaks=seq(0,2*N,length=21),
     xlab="Introgressed haplotypes",ylab="SNP count",main="Number of introgressed individuals per site (before captive program)")
hist(rowSums(data[[1]]),breaks=seq(0,L,length=21),
     xlab="Introgressed sites",ylab="haplotype count",main="Number of introgressed sites per haplotype (after captive program)")
hist(colSums(data[[1]]),breaks=seq(0,2*N,length=21),
     xlab="Introgressed haplotypes",ylab="SNP count",main="Number of introgressed individuals per site (after captive program)")
allScores(0,data0)
allScores(0,data)

mean(data0[[1]])
par(mfrow=c(2,2))
image(t(data0[[1]]))
image(t(data0[[2]]))
image(t(data[[1]]))
image(t(data[[2]]))

### Breeding step


Gforward=10
dataNull<-breedingSel(data,Gforward,selfn=selIdentity)
dataRandom<-breedingSel(data,Gforward,selfn=selRandom)
dataQ<-breedingSel(data,Gforward,selfn=selMinQ)
dataHet<-breedingSel(data,Gforward,selfn=selHet)
dataHetQ<-breedingSel(data,Gforward,selfn=selHetQ)

scores=c("Q","Het","HetQ")
names=c("Domestic Admixture","Heterozygosity","Heterozygosity for wildcat loci")
names(names)=scores

pdf("SimulatedBreeding.pdf",height=6,width=12)
par(mfrow=c(1,3))
for(score in scores){
    yrange=c(0,0.35)
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
