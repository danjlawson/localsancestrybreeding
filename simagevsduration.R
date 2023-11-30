source("localancestryfns.R")
### Start of sim

set.seed(1)
nrep=10 ## How many sims per datapoint?

rangehigh<-function(x)quantile(x,0.75)
rangelow<-function(x)quantile(x,0.25)

########################################
## gens since introgression scan

## genetic distance
genvals=c(2,(2:10)^2)
nrep=10
genall=expand.grid(genback=genvals,genforward=genvals)
indices=rep(1:dim(genall)[1],each=nrep)
genall=genall[indices,]

gendatalist<-lapply(1:dim(genall)[1],function(i){
    simCaptive(N=100,NI=50,L=500,seed=1+(i%%nrep),G=genall[i,1])
})

## Simulate breeding program
gen2MHQScore=lapply(1:dim(genall)[1],function(i){
    print(i)
    Gforward=genall[i,2]
    dataPairwiseMHQScore<-breedingSel(gendatalist[[i]],Gforward,
                                      selfn=pairwiseHeterozygosityQScore,parentfn=rankedParents,
                                      min=FALSE)
    dataPairwiseMHQScore
})
gen2MKQ=lapply(1:dim(genall)[1],function(i){
    print(i)
    Gforward=genall[i,2]
    dataPairwiseMHQScore<-breedingSel(gendatalist[[i]],Gforward,
                                      selfn=pairwiseKinshipQ,parentfn=rankedParents,
                                      min=TRUE)
    dataPairwiseMHQScore
})

## Construct score
gen2scoreres=data.frame(t(sapply(1:length(gen2MHQScore),function(i){
    score=tail(gen2MHQScore[[i]]$score,1)
    score$genback = genall[i,1]
    score$genforward = genall[i,2]
    as.numeric(score)
})))
colnames(gen2scoreres)=c(colnames(gen2MHQScore[[1]]$score),"genback","genforward")

## Construct score
gen2Kscoreres=data.frame(t(sapply(1:length(gen2MKQ),function(i){
    score=tail(gen2MKQ[[i]]$score,1)
    score$genback = genall[i,1]
    score$genforward = genall[i,2]
    as.numeric(score)
})))
colnames(gen2Kscoreres)=c(colnames(gen2MHQScore[[1]]$score),"genback","genforward")

save.image("simagevsduration.RData")
## On the HPC
q(save=FALSE)


gen2Kscorelist=lapply(1:dim(gen2Kscoreres)[2],function(i){
    r=matrix(gen2Kscoreres[,i],nrow=10,ncol=10)
    rownames(r)=colnames(r)=genvals
    r
})
names(gen2Kscorelist)=colnames(gen2Kscoreres)

gen2scorelist=lapply(1:dim(gen2scoreres)[2],function(i){
    r=matrix(gen2scoreres[,i],nrow=10,ncol=10)
    rownames(r)=colnames(r)=genvals
    r
})
names(gen2scorelist)=colnames(gen2scoreres)


genscoremeans=aggregate.data.frame(genscoreres,by=list(groupl=genscoreres$gen),mean)
genscoremins=aggregate.data.frame(genscoreres,by=list(groupl=genscoreres$gen),rangelow)
genscoremaxs=aggregate.data.frame(genscoreres,by=list(groupl=genscoreres$gen),rangehigh)

genKscoremeans=aggregate.data.frame(genKscoreres,by=list(groupl=genKscoreres$gen),mean)
genKscoremins=aggregate.data.frame(genKscoreres,by=list(groupl=genKscoreres$gen),rangelow)
genKscoremaxs=aggregate.data.frame(genKscoreres,by=list(groupl=genKscoreres$gen),rangehigh)


# remotes::install_github("danjlawson/CLARITY/Clarity")
library("Clarity")
## How to check x and y axes are correct
## Clarity_Chart(gen2Kscorelist[["genforward"]],text=T)
par(mfrow=c(1,2))
Clarity_Chart(1-gen2Kscorelist[["Q"]],text=T,xlab="Generations back",ylab="Generations forward",main="Kinship")
Clarity_Chart(1-gen2scorelist[["Q"]],text=T,xlab="Generations back",ylab="Generations forward",main="Heterozygosity")

gen2Kscorelist=lapply(1:dim(gen2Kscoreres)[2],function(i) matrix(gen2Kscoreres[,i],nrow=10,ncol=10))
matrix(gen2Kscoreres[,"genforward"],nrow=10,ncol=10)
