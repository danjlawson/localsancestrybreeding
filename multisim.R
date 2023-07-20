source("localancestryfns.R")

nrep=10 ## How many sims per datapoint?

rangehigh<-function(x)quantile(x,0.75)
rangelow<-function(x)quantile(x,0.25)

########################################
## genetic distance scan

## genetic distance
gLtest=rep(seq(1,10,1)^2,each=nrep)
gdistlist<-lapply(gLtest,function(gL,L=500){
    gdist=diff(seq(0,gL,length.out=L+1))})
## Sim 1: we stick with the same genetic structure; every case gets a different expected admixture fraction
gdatalist<-lapply(1:length(gLtest),function(i){
    x=gLtest[i]
    simCaptive(N=100,NI=50,L=500,gL=x,seed=i)
})
## Sim 2: No change to the genetic structure in the panel, only post-panel
gdatalist<-lapply(1:length(gLtest),function(i){
    x=gLtest[i]
    simCaptive(N=100,NI=50,L=500,seed=1+(i%%nrep))
})


tmp= simCaptive(N=100,NI=50,L=500,seed=1)
galpha0=sapply(gdatalist,function(x)mean(x$panel))

## Simulate breeding program
Gforward=20
gLMHQScore=lapply(1:length(gLtest),function(i){
    print(i)
    dataPairwiseMHQScore<-breedingSel(gdatalist[[i]],Gforward,d=gdistlist[[i]],
                                      selfn=pairwiseHeterozygosityQScore,parentfn=rankedParents,min=FALSE)
    dataPairwiseMHQScore
})
gLMKQ=lapply(1:length(gLtest),function(i){
    print(i)
    dataPairwiseMHQScore<-breedingSel(gdatalist[[i]],Gforward,d=gdistlist[[i]],
                                      selfn=pairwiseKinshipQ,parentfn=rankedParents,min=TRUE)
    dataPairwiseMHQScore
})

## Construct score
gLscoreres=data.frame(t(sapply(1:length(gLMHQScore),function(i){
    gL=gLtest[i]
    score=tail(gLMHQScore[[i]]$score,1)
    score$gL=gL
    as.numeric(score)
})))
colnames(gLscoreres)=c(colnames(gLMHQScore[[1]]$score),"gL")

gLscoremeans=aggregate.data.frame(gLscoreres,by=list(groupgL=gLscoreres$gL),mean)
gLscoremins=aggregate.data.frame(gLscoreres,by=list(groupgL=gLscoreres$gL),rangelow)
gLscoremaxs=aggregate.data.frame(gLscoreres,by=list(groupgL=gLscoreres$gL),rangehigh)

## Construct score for Kinship
gLKscoreres=data.frame(t(sapply(1:length(gLMKQ),function(i){
    gL=gLtest[i]
    score=tail(gLMKQ[[i]]$score,1)
    score$gL=gL
    as.numeric(score)
})))
colnames(gLKscoreres)=c(colnames(gLMKQ[[1]]$score),"gL")

gLKscoremeans=aggregate.data.frame(gLKscoreres,by=list(groupgL=gLKscoreres$gL),mean)
gLKscoremins=aggregate.data.frame(gLKscoreres,by=list(groupgL=gLKscoreres$gL),rangelow)
gLKscoremaxs=aggregate.data.frame(gLKscoreres,by=list(groupgL=gLKscoreres$gL),rangehigh)


########################################
## admixture scan

## admixture sims
atest=rep(seq(0.05,0.5,length.out=10),each=nrep)
adatalist<-lapply(1:length(gLtest),function(i){
    x=gLtest[i]
    simCaptive(N=100,NI=50,L=500,alpha=atest[i],seed=i)
})
aalpha0=sapply(adatalist,function(x)mean(x$panel))

## Simulate breeding program
Gforward=20
aMHQScore=lapply(1:length(atest),function(i){
    print(i)
    dataPairwiseMHQScore<-breedingSel(adatalist[[i]],Gforward,
                                      selfn=pairwiseHeterozygosityQScore,parentfn=rankedParents,min=FALSE)
    dataPairwiseMHQScore
})
aMKQ=lapply(1:length(atest),function(i){
    print(i)
    dataPairwiseMHQScore<-breedingSel(adatalist[[i]],Gforward,
                                      selfn=pairwiseKinshipQ,parentfn=rankedParents,min=TRUE)
    dataPairwiseMHQScore
})

## Construct score
alphascoreres=data.frame(t(sapply(1:length(aMHQScore),function(i){
    alpha=atest[i]
    score=tail(aMHQScore[[i]]$score,1)
    score$alpha=alpha
    as.numeric(score)
})))
colnames(alphascoreres)=c(colnames(aMHQScore[[1]]$score),"alpha")

alphaxvals=aggregate(1-aalpha0,list(groupalpha=alphascoreres$alpha),mean)
alphascoreres$NonQ0=1-aalpha0
alphascoremeans=aggregate.data.frame(alphascoreres,
                                     by=list(groupalpha=alphascoreres$alpha),mean)
alphascoremins=aggregate.data.frame(alphascoreres,
                                    by=list(groupalpha=alphascoreres$alpha),rangelow)
alphascoremaxs=aggregate.data.frame(alphascoreres,
                                    by=list(groupalpha=alphascoreres$alpha),rangehigh)
to=order(alphascoremeans$NonQ0)
alphascoremeans=alphascoremeans[to,]
alphascoremins=alphascoremins[to,]
alphascoremaxs=alphascoremaxs[to,]


## Construct scorefor Kinship
alphaKscoreres=data.frame(t(sapply(1:length(aMKQ),function(i){
    alpha=atest[i]
    score=tail(aMKQ[[i]]$score,1)
    score$alpha=alpha
    as.numeric(score)
})))
colnames(alphaKscoreres)=c(colnames(aMKQ[[1]]$score),"alpha")

alphaKxvals=aggregate(1-aalpha0,list(groupalpha=alphaKscoreres$alpha),mean)
alphaKscoreres$NonQ0=1-aalpha0
alphaKscoremeans=aggregate.data.frame(alphaKscoreres,
                                     by=list(groupalpha=alphaKscoreres$alpha),mean)
alphaKscoremins=aggregate.data.frame(alphaKscoreres,
                                    by=list(groupalpha=alphaKscoreres$alpha),rangelow)
alphaKscoremaxs=aggregate.data.frame(alphaKscoreres,
                                    by=list(groupalpha=alphaKscoreres$alpha),rangehigh)
to=order(alphaKscoremeans$NonQ0)
alphaKscoremeans=alphaKscoremeans[to,]
alphaKscoremins=alphaKscoremins[to,]
alphaKscoremaxs=alphaKscoremaxs[to,]

########################################
## broodsize scan

## genetic distance
lambdatest=rep(seq(2.1,8,length.out=15),each=nrep)
## Sim 2: No change to the genetic structure in the panel, only post-panel
ldatalist<-lapply(1:length(lambdatest),function(i){
    x=gLtest[i]
    simCaptive(N=100,NI=50,L=500,seed=1+(i%%nrep))
})
lfunctionlist=lapply(lambdatest,function(x){
    function(N)rpois(N,x)
})

## Simulate breeding program
Gforward=20
lMHQScore=lapply(1:length(lambdatest),function(i){
    print(i)
    dataPairwiseMHQScore<-breedingSel(ldatalist[[i]],Gforward,
                                      selfn=pairwiseHeterozygosityQScore,parentfn=rankedParents,
                                      offspringdist=lfunctionlist[[i]],min=FALSE)
    dataPairwiseMHQScore
})
lMKQ=lapply(1:length(lambdatest),function(i){
    print(i)
    dataPairwiseMHQScore<-breedingSel(ldatalist[[i]],Gforward,
                                      selfn=pairwiseKinshipQ,parentfn=rankedParents,
                                      offspringdist=lfunctionlist[[i]],min=TRUE)
    dataPairwiseMHQScore
})

## Construct score
lambdascoreres=data.frame(t(sapply(1:length(lMHQScore),function(i){
    lamndafn=gLtest[i]
    score=tail(lMHQScore[[i]]$score,1)
    score$lambda = lambdatest[i]
    as.numeric(score)
})))
colnames(lambdascoreres)=c(colnames(lMHQScore[[1]]$score),"lambda")
lambdascoreres=na.omit(lambdascoreres)

lambdascoremeans=aggregate.data.frame(lambdascoreres,by=list(groupl=lambdascoreres$lambda),mean)
lambdascoremins=aggregate.data.frame(lambdascoreres,by=list(groupl=lambdascoreres$lambda),rangelow)
lambdascoremaxs=aggregate.data.frame(lambdascoreres,by=list(groupl=lambdascoreres$lambda),rangehigh)

## Construct score
lambdaKscoreres=data.frame(t(sapply(1:length(lMKQ),function(i){
    lamndafn=gLtest[i]
    score=tail(lMKQ[[i]]$score,1)
    score$lambda = lambdatest[i]
    as.numeric(score)
})))
colnames(lambdaKscoreres)=c(colnames(lMHQScore[[1]]$score),"lambda")
lambdaKscoreres=na.omit(lambdaKscoreres)

lambdaKscoremeans=aggregate.data.frame(lambdaKscoreres,by=list(groupl=lambdaKscoreres$lambda),mean)
lambdaKscoremins=aggregate.data.frame(lambdaKscoreres,by=list(groupl=lambdaKscoreres$lambda),rangelow)
lambdaKscoremaxs=aggregate.data.frame(lambdaKscoreres,by=list(groupl=lambdaKscoreres$lambda),rangehigh)


########################################
## gens since introgression scan

## genetic distance
gentest=rep( c(2,(2:10)^2),each=nrep)
## Sim 2: No change to the genetic structure in the panel, only post-panel
gendatalist<-lapply(1:length(gentest),function(i){
    simCaptive(N=100,NI=50,L=500,seed=1+(i%%nrep),G=gentest[i])
})

## Simulate breeding program
Gforward=20
genMHQScore=lapply(1:length(gentest),function(i){
    print(i)
    dataPairwiseMHQScore<-breedingSel(gendatalist[[i]],Gforward,
                                      selfn=pairwiseHeterozygosityQScore,parentfn=rankedParents,
                                      min=FALSE)
    dataPairwiseMHQScore
})
genMKQ=lapply(1:length(gentest),function(i){
    print(i)
    dataPairwiseMHQScore<-breedingSel(gendatalist[[i]],Gforward,
                                      selfn=pairwiseKinshipQ,parentfn=rankedParents,
                                      min=TRUE)
    dataPairwiseMHQScore
})

## Construct score
genscoreres=data.frame(t(sapply(1:length(genMHQScore),function(i){
    score=tail(genMHQScore[[i]]$score,1)
    score$gen = gentest[i]
    as.numeric(score)
})))
colnames(genscoreres)=c(colnames(genMHQScore[[1]]$score),"gen")
genscoreres=na.omit(genscoreres)

genscoremeans=aggregate.data.frame(genscoreres,by=list(groupl=genscoreres$gen),mean)
genscoremins=aggregate.data.frame(genscoreres,by=list(groupl=genscoreres$gen),rangelow)
genscoremaxs=aggregate.data.frame(genscoreres,by=list(groupl=genscoreres$gen),rangehigh)

## Construct score
genKscoreres=data.frame(t(sapply(1:length(genMKQ),function(i){
    score=tail(genMKQ[[i]]$score,1)
    score$gen = gentest[i]
    as.numeric(score)
})))
colnames(genKscoreres)=c(colnames(genMHQScore[[1]]$score),"gen")
genKscoreres=na.omit(genKscoreres)

genKscoremeans=aggregate.data.frame(genKscoreres,by=list(groupl=genKscoreres$gen),mean)
genKscoremins=aggregate.data.frame(genKscoreres,by=list(groupl=genKscoreres$gen),rangelow)
genKscoremaxs=aggregate.data.frame(genKscoreres,by=list(groupl=genKscoreres$gen),rangehigh)


save.image("multisim.RData")

load("multisim.RData")

############################################
## plot
myline<-function(x,ymean,ymin,ymax,col,...){
    pcol=paste0(col,"55")
    polygon(c(x,rev(x)),c(ymin,rev(ymax)),border=NA,col=pcol,...)
    lines(x,ymean,col=col,lwd=2,...)
}

########################
## Figure 4
pdf("Figure4-ChangingSim-version0.pdf",height=4,width=10)
par(las=1)
layout(matrix(c(1,2,3),nrow=1))
plot(alphascoremeans$NonQ0,alphascoremeans$NonQ,type="n",col=2,ylim=c(0,0.5),xlab="Founder introgression proportion",ylab="Score after 20 Generations",cex.axis=1.5,cex.lab=1.2)
abline(v=0.2,lty=3)
myline(alphascoremeans$NonQ0,alphascoremeans$NonQ,alphascoremins$NonQ,alphascoremaxs$NonQ,col="#FF0000")
myline(alphascoremeans$NonQ0,alphascoremeans$Het,alphascoremins$Het,alphascoremaxs$Het,col="#00FF00")
myline(alphascoremeans$NonQ0,alphascoremeans$Kin,alphascoremins$Kin,alphascoremaxs$Kin,col="#00FFFF")
myline(alphascoremeans$NonQ0,alphascoremeans$HetQ,alphascoremins$HetQ,alphascoremaxs$HetQ,col="#0000FF")
myline(alphascoremeans$NonQ0,alphascoremeans$KinQ,alphascoremins$KinQ,alphascoremaxs$KinQ,col="#FF00FF")
mtext("a) Impact of Introgression Magnitude",adj=-0.2,line=1,cex=0.8)
##
plot(gLscoremeans$gL,gLscoremeans$NonQ,type="n",col=2,ylim=c(0,0.5),xlab="Genome Length (Morgans)",ylab="Score after 20 Generations",cex.axis=1.5,cex.lab=1.2)
abline(v=5,lty=3)
myline(gLscoremeans$gL,gLscoremeans$NonQ,gLscoremins$NonQ,gLscoremaxs$NonQ,col="#FF0000")
myline(gLscoremeans$gL,gLscoremeans$Het,gLscoremins$Het,gLscoremaxs$Het,col="#00FF00")
myline(gLscoremeans$gL,gLscoremeans$Kin,gLscoremins$Kin,gLscoremaxs$Kin,col="#00FFFF")
myline(gLscoremeans$gL,gLscoremeans$HetQ,gLscoremins$HetQ,gLscoremaxs$HetQ,col="#0000FF")
myline(gLscoremeans$gL,gLscoremeans$KinQ,gLscoremins$KinQ,gLscoremaxs$KinQ,col="#FF00FF")
mtext("b) Impact of Genome Length",adj=-0.2,line=1,cex=0.8)
plot(lambdascoremeans$lambda,lambdascoremeans$NonQ,type="n",col=2,ylim=c(0,0.5),xlab="Mean offspring per pair surviving to adulthood",ylab="Score after 20 Generations",cex.axis=1.5,cex.lab=1.2)
abline(v=3,lty=3)
myline(lambdascoremeans$lambda,lambdascoremeans$NonQ,lambdascoremins$NonQ,lambdascoremaxs$NonQ,col="#FF0000")
myline(lambdascoremeans$lambda,lambdascoremeans$Het,lambdascoremins$Het,lambdascoremaxs$Het,col="#00FF00")
myline(lambdascoremeans$lambda,lambdascoremeans$Kin,lambdascoremins$Kin,lambdascoremaxs$Kin,col="#00FFFF")
myline(lambdascoremeans$lambda,lambdascoremeans$HetQ,lambdascoremins$HetQ,lambdascoremaxs$HetQ,col="#0000FF")
myline(lambdascoremeans$lambda,lambdascoremeans$KinQ,lambdascoremins$KinQ,lambdascoremaxs$KinQ,col="#FF00FF")
mtext("c) Impact of Mean offspring",adj=-0.2,line=1,cex=0.8)
##
## par(mar=c(0,0,0,0))
## plot(c(0,1),c(0,1), type="n",xlab="",ylab="",axes=FALSE)
legend("topright",legend=c("Introgression","Heterozygosity","Kinship",
                          "Population Heterozygosity PH",
                         "Population Kinship PK"),
       text.col=c(2,3,5,4,6),lty=1,col=c(2,3,5,4,6),title="Score:",title.col=1,cex=1,lwd=2,bg="white")
##
dev.off()


########################
## Figure 4
pdf("Figure4-ChangingSim-version0.1.pdf",height=8,width=10)
par(las=1)
layout(matrix(c(1:6),nrow=2,ncol=3,byrow=T))
plot(alphascoremeans$NonQ0,alphascoremeans$NonQ,type="n",col=2,ylim=c(0,0.5),xlab="Founder introgression proportion",ylab="Score after 20 Generations",cex.axis=1.5,cex.lab=1.2)
abline(v=0.2,lty=3)
myline(alphascoremeans$NonQ0,alphascoremeans$NonQ,alphascoremins$NonQ,alphascoremaxs$NonQ,col="#FF0000")
myline(alphascoremeans$NonQ0,alphascoremeans$Het,alphascoremins$Het,alphascoremaxs$Het,col="#00FF00")
myline(alphascoremeans$NonQ0,alphascoremeans$Kin,alphascoremins$Kin,alphascoremaxs$Kin,col="#00FFFF")
myline(alphascoremeans$NonQ0,alphascoremeans$HetQ,alphascoremins$HetQ,alphascoremaxs$HetQ,col="#0000FF")
myline(alphascoremeans$NonQ0,alphascoremeans$KinQ,alphascoremins$KinQ,alphascoremaxs$KinQ,col="#FF00FF")
mtext("a) Introgression (Selecting WPH)",adj=0,line=1,cex=0.8)
##
plot(gLscoremeans$gL,gLscoremeans$NonQ,type="n",col=2,ylim=c(0,0.5),xlab="Genome Length (Morgans)",ylab="Score after 20 Generations",cex.axis=1.5,cex.lab=1.2)
abline(v=5,lty=3)
myline(gLscoremeans$gL,gLscoremeans$NonQ,gLscoremins$NonQ,gLscoremaxs$NonQ,col="#FF0000")
myline(gLscoremeans$gL,gLscoremeans$Het,gLscoremins$Het,gLscoremaxs$Het,col="#00FF00")
myline(gLscoremeans$gL,gLscoremeans$Kin,gLscoremins$Kin,gLscoremaxs$Kin,col="#00FFFF")
myline(gLscoremeans$gL,gLscoremeans$HetQ,gLscoremins$HetQ,gLscoremaxs$HetQ,col="#0000FF")
myline(gLscoremeans$gL,gLscoremeans$KinQ,gLscoremins$KinQ,gLscoremaxs$KinQ,col="#FF00FF")
mtext("b) Genome Length (Selecting WPH)",adj=-0.2,line=1,cex=0.8)
plot(lambdascoremeans$lambda,lambdascoremeans$NonQ,type="n",col=2,ylim=c(0,0.5),xlim=c(2,6),xlab="Mean offspring per pair surviving to adulthood",ylab="Score after 20 Generations",cex.axis=1.5,cex.lab=1.2)
abline(v=3,lty=3)
myline(lambdascoremeans$lambda,lambdascoremeans$NonQ,lambdascoremins$NonQ,lambdascoremaxs$NonQ,col="#FF0000")
myline(lambdascoremeans$lambda,lambdascoremeans$Het,lambdascoremins$Het,lambdascoremaxs$Het,col="#00FF00")
myline(lambdascoremeans$lambda,lambdascoremeans$Kin,lambdascoremins$Kin,lambdascoremaxs$Kin,col="#00FFFF")
myline(lambdascoremeans$lambda,lambdascoremeans$HetQ,lambdascoremins$HetQ,lambdascoremaxs$HetQ,col="#0000FF")
myline(lambdascoremeans$lambda,lambdascoremeans$KinQ,lambdascoremins$KinQ,lambdascoremaxs$KinQ,col="#FF00FF")
mtext("c) Mean offspring (Selecting WPH)",adj=-0.2,line=1,cex=0.8)
legend("topright",legend=c("Introgression","Heterozygosity","Kinship",
                          "Population Heterozygosity PH",
                         "Population Kinship PK"),
       text.col=c(2,3,5,4,6),lty=1,col=c(2,3,5,4,6),title="Score:",title.col=1,cex=1,lwd=2,bg="white")
##
plot(alphaKscoremeans$NonQ0,alphaKscoremeans$NonQ,type="n",col=2,ylim=c(0,0.5),xlab="Founder introgression proportion",ylab="Score after 20 Generations",cex.axis=1.5,cex.lab=1.2)
abline(v=0.2,lty=3)
myline(alphaKscoremeans$NonQ0,alphaKscoremeans$NonQ,alphaKscoremins$NonQ,alphaKscoremaxs$NonQ,col="#FF0000")
myline(alphaKscoremeans$NonQ0,alphaKscoremeans$Het,alphaKscoremins$Het,alphaKscoremaxs$Het,col="#00FF00")
myline(alphaKscoremeans$NonQ0,alphaKscoremeans$Kin,alphaKscoremins$Kin,alphaKscoremaxs$Kin,col="#00FFFF")
myline(alphaKscoremeans$NonQ0,alphaKscoremeans$HetQ,alphaKscoremins$HetQ,alphaKscoremaxs$HetQ,col="#0000FF")
myline(alphaKscoremeans$NonQ0,alphaKscoremeans$KinQ,alphaKscoremins$KinQ,alphaKscoremaxs$KinQ,col="#FF00FF")
mtext("d) Introgression (Selecting PK)",adj=-0.2,line=1,cex=0.8)
##
plot(gLKscoremeans$gL,gLKscoremeans$NonQ,type="n",col=2,ylim=c(0,0.5),xlab="Genome Length (Morgans)",ylab="Score after 20 Generations",cex.axis=1.5,cex.lab=1.2)
abline(v=5,lty=3)
myline(gLKscoremeans$gL,gLKscoremeans$NonQ,gLKscoremins$NonQ,gLKscoremaxs$NonQ,col="#FF0000")
myline(gLKscoremeans$gL,gLKscoremeans$Het,gLKscoremins$Het,gLKscoremaxs$Het,col="#00FF00")
myline(gLKscoremeans$gL,gLKscoremeans$Kin,gLKscoremins$Kin,gLKscoremaxs$Kin,col="#00FFFF")
myline(gLKscoremeans$gL,gLKscoremeans$HetQ,gLKscoremins$HetQ,gLKscoremaxs$HetQ,col="#0000FF")
myline(gLKscoremeans$gL,gLKscoremeans$KinQ,gLKscoremins$KinQ,gLKscoremaxs$KinQ,col="#FF00FF")
mtext("e) Genome Length (Selecting PK)",adj=-0.2,line=1,cex=0.8)
plot(lambdaKscoremeans$lambda,lambdaKscoremeans$NonQ,type="n",col=2,ylim=c(0,0.5),xlim=c(2,6),xlab="Mean offspring per pair surviving to adulthood",ylab="Score after 20 Generations",cex.axis=1.5,cex.lab=1.2)
abline(v=3,lty=3)
myline(lambdaKscoremeans$lambda,lambdaKscoremeans$NonQ,lambdaKscoremins$NonQ,lambdaKscoremaxs$NonQ,col="#FF0000")
myline(lambdaKscoremeans$lambda,lambdaKscoremeans$Het,lambdaKscoremins$Het,lambdaKscoremaxs$Het,col="#00FF00")
myline(lambdaKscoremeans$lambda,lambdaKscoremeans$Kin,lambdaKscoremins$Kin,lambdaKscoremaxs$Kin,col="#00FFFF")
myline(lambdaKscoremeans$lambda,lambdaKscoremeans$HetQ,lambdaKscoremins$HetQ,lambdaKscoremaxs$HetQ,col="#0000FF")
myline(lambdaKscoremeans$lambda,lambdaKscoremeans$KinQ,lambdaKscoremins$KinQ,lambdaKscoremaxs$KinQ,col="#FF00FF")
mtext("f) Mean Offspring (Selecting PK)",adj=-0.2,line=1,cex=0.8)
##
## par(mar=c(0,0,0,0))
## plot(c(0,1),c(0,1), type="n",xlab="",ylab="",axes=FALSE)
legend("topright",legend=c("Introgression","Heterozygosity","Kinship",
                          "Population Heterozygosity PH",
                         "Population Kinship PK"),
text.col=c(2,3,5,4,6),lty=1,col=c(2,3,5,4,6),title="Score:",title.col=1,cex=1,lwd=2,bg="white")
##
dev.off()



###############################
cex.lab=1.5
cex.axis=1.5
cex.mtext=0.95
pdf("Figure4-ChangingSim.pdf",height=7,width=12)
par(las=1,cex=0.8,mar=c(5,5,3,1))
layout(matrix(c(1:8),nrow=2,ncol=4,byrow=T))
plot(alphascoremeans$NonQ0,alphascoremeans$NonQ,type="n",col=2,ylim=c(0,0.5),xlab="Founder introgression proportion",ylab="Score after 20 Generations",cex.axis=cex.axis,cex.lab=cex.lab)
abline(v=0.2,lty=3)
myline(alphascoremeans$NonQ0,alphascoremeans$NonQ,alphascoremins$NonQ,alphascoremaxs$NonQ,col="#FF0000")
myline(alphascoremeans$NonQ0,alphascoremeans$Het,alphascoremins$Het,alphascoremaxs$Het,col="#00FF00")
myline(alphascoremeans$NonQ0,alphascoremeans$Kin,alphascoremins$Kin,alphascoremaxs$Kin,col="#00FFFF")
myline(alphascoremeans$NonQ0,alphascoremeans$KinQ,alphascoremins$KinQ,alphascoremaxs$KinQ,col="#0000FF")
myline(alphascoremeans$NonQ0,alphascoremeans$HetQ,alphascoremins$HetQ,alphascoremaxs$HetQ,col="#FF00FF")
mtext("a) Introgression (Selecting WPH)",adj=0.6,line=1,cex=cex.mtext)
plot(gLscoremeans$gL,gLscoremeans$NonQ,type="n",col=2,ylim=c(0,0.5),xlab="Genome Length (Morgans)",ylab="Score after 20 Generations",cex.axis=cex.axis,cex.lab=cex.lab,log="x")
abline(v=5,lty=3)
myline(gLscoremeans$gL,gLscoremeans$NonQ,gLscoremins$NonQ,gLscoremaxs$NonQ,col="#FF0000")
myline(gLscoremeans$gL,gLscoremeans$Het,gLscoremins$Het,gLscoremaxs$Het,col="#00FF00")
myline(gLscoremeans$gL,gLscoremeans$Kin,gLscoremins$Kin,gLscoremaxs$Kin,col="#00FFFF")
myline(gLscoremeans$gL,gLscoremeans$KinQ,gLscoremins$KinQ,gLscoremaxs$KinQ,col="#0000FF")
myline(gLscoremeans$gL,gLscoremeans$HetQ,gLscoremins$HetQ,gLscoremaxs$HetQ,col="#FF00FF")
mtext("b) Genome Length (Selecting WPH)",adj=0.6,line=1,cex=cex.mtext)
plot(lambdascoremeans$lambda,lambdascoremeans$NonQ,type="n",col=2,ylim=c(0,0.5),xlim=c(2,6),xlab="Mean breeding offspring per pair",ylab="Score after 20 Generations",cex.axis=1.5,cex.lab=cex.lab)
abline(v=4.2,lty=3)
myline(lambdascoremeans$lambda,lambdascoremeans$NonQ,lambdascoremins$NonQ,lambdascoremaxs$NonQ,col="#FF0000")
myline(lambdascoremeans$lambda,lambdascoremeans$Het,lambdascoremins$Het,lambdascoremaxs$Het,col="#00FF00")
myline(lambdascoremeans$lambda,lambdascoremeans$Kin,lambdascoremins$Kin,lambdascoremaxs$Kin,col="#00FFFF")
myline(lambdascoremeans$lambda,lambdascoremeans$KinQ,lambdascoremins$KinQ,lambdascoremaxs$KinQ,col="#0000FF")
myline(lambdascoremeans$lambda,lambdascoremeans$HetQ,lambdascoremins$HetQ,lambdascoremaxs$HetQ,col="#FF00FF")
mtext("c) Mean Offspring (Selecting WPH)",adj=0.6,line=1,cex=cex.mtext)
##
plot(genscoremeans$gen,genscoremeans$NonQ,type="n",col=2,ylim=c(0,0.5),xlab="Gens before breeding program",ylab="Score after 20 Generations",cex.axis=1.5,cex.lab=cex.lab,log="x")
abline(v=10,lty=3)
myline(genscoremeans$gen,genscoremeans$NonQ,genscoremins$NonQ,genscoremaxs$NonQ,col="#FF0000")
myline(genscoremeans$gen,genscoremeans$Het,genscoremins$Het,genscoremaxs$Het,col="#00FF00")
myline(genscoremeans$gen,genscoremeans$Kin,genscoremins$Kin,genscoremaxs$Kin,col="#00FFFF")
myline(genscoremeans$gen,genscoremeans$KinQ,genscoremins$KinQ,genscoremaxs$KinQ,col="#0000FF")
myline(genscoremeans$gen,genscoremeans$HetQ,genscoremins$HetQ,genscoremaxs$HetQ,col="#FF00FF")
mtext("d) Introgression Age (Selecting WPH)",adj=0.6,line=1,cex=cex.mtext)
##
plot(alphaKscoremeans$NonQ0,alphaKscoremeans$NonQ,type="n",col=2,ylim=c(0,0.5),xlab="Founder introgression proportion",ylab="Score after 20 Generations",cex.axis=1.5,cex.lab=cex.lab)
abline(v=0.2,lty=3)
myline(alphaKscoremeans$NonQ0,alphaKscoremeans$NonQ,alphaKscoremins$NonQ,alphaKscoremaxs$NonQ,col="#FF0000")
myline(alphaKscoremeans$NonQ0,alphaKscoremeans$Het,alphaKscoremins$Het,alphaKscoremaxs$Het,col="#00FF00")
myline(alphaKscoremeans$NonQ0,alphaKscoremeans$Kin,alphaKscoremins$Kin,alphaKscoremaxs$Kin,col="#00FFFF")
myline(alphaKscoremeans$NonQ0,alphaKscoremeans$KinQ,alphaKscoremins$KinQ,alphaKscoremaxs$KinQ,col="#0000FF")
myline(alphaKscoremeans$NonQ0,alphaKscoremeans$HetQ,alphaKscoremins$HetQ,alphaKscoremaxs$HetQ,col="#FF00FF")
mtext("e) Introgression (Selecting PK)",adj=0.6,line=1,cex=cex.mtext)
##
legend("topleft",
       legend=c(expression("Introgression"~S^Q),
                expression("Heterozygosity"~S^H),
                expression("Kinship"~S^K),
                expression("Pop. Kinship"~S^{PK}),
                expression("Pop. Heterozygosity"~S^{PH})),
       text.col=c(2,3,5,4,6),lty=1,
       col=c(2,3,5,4,6),title="Score:",title.col=1,cex=1.2,lwd=2,bg="white")
##
##
plot(gLKscoremeans$gL,gLKscoremeans$NonQ,type="n",col=2,ylim=c(0,0.5),xlab="Genome Length (Morgans)",ylab="Score after 20 Generations",cex.axis=1.5,cex.lab=cex.lab,log="x")
abline(v=5,lty=3)
myline(gLKscoremeans$gL,gLKscoremeans$NonQ,gLKscoremins$NonQ,gLKscoremaxs$NonQ,col="#FF0000")
myline(gLKscoremeans$gL,gLKscoremeans$Het,gLKscoremins$Het,gLKscoremaxs$Het,col="#00FF00")
myline(gLKscoremeans$gL,gLKscoremeans$Kin,gLKscoremins$Kin,gLKscoremaxs$Kin,col="#00FFFF")
myline(gLKscoremeans$gL,gLKscoremeans$KinQ,gLKscoremins$KinQ,gLKscoremaxs$KinQ,col="#0000FF")
myline(gLKscoremeans$gL,gLKscoremeans$HetQ,gLKscoremins$HetQ,gLKscoremaxs$HetQ,col="#FF00FF")
mtext("f) Genome Length (Selecting PK)",adj=0.6,line=1,cex=cex.mtext)
plot(lambdaKscoremeans$lambda,lambdaKscoremeans$NonQ,type="n",col=2,ylim=c(0,0.5),xlim=c(2,6),xlab="Mean breeding offspring per pair",ylab="Score after 20 Generations",cex.axis=1.5,cex.lab=cex.lab)
abline(v=4.2,lty=3)
myline(lambdaKscoremeans$lambda,lambdaKscoremeans$NonQ,lambdaKscoremins$NonQ,lambdaKscoremaxs$NonQ,col="#FF0000")
myline(lambdaKscoremeans$lambda,lambdaKscoremeans$Het,lambdaKscoremins$Het,lambdaKscoremaxs$Het,col="#00FF00")
myline(lambdaKscoremeans$lambda,lambdaKscoremeans$Kin,lambdaKscoremins$Kin,lambdaKscoremaxs$Kin,col="#00FFFF")
myline(lambdaKscoremeans$lambda,lambdaKscoremeans$KinQ,lambdaKscoremins$KinQ,lambdaKscoremaxs$KinQ,col="#0000FF")
myline(lambdaKscoremeans$lambda,lambdaKscoremeans$HetQ,lambdaKscoremins$HetQ,lambdaKscoremaxs$HetQ,col="#FF00FF")
mtext("g) Mean Offspring (Selecting PK)",adj=0.6,line=1,cex=cex.mtext)
##
plot(genKscoremeans$gen,genKscoremeans$NonQ,type="n",col=2,ylim=c(0,0.5),xlab="Gens before breeding program",ylab="Score after 20 Generations",cex.axis=1.5,cex.lab=cex.lab,log="x")
abline(v=10,lty=3)
myline(genKscoremeans$gen,genKscoremeans$NonQ,genKscoremins$NonQ,genKscoremaxs$NonQ,col="#FF0000")
myline(genKscoremeans$gen,genKscoremeans$Het,genKscoremins$Het,genKscoremaxs$Het,col="#00FF00")
myline(genKscoremeans$gen,genKscoremeans$Kin,genKscoremins$Kin,genKscoremaxs$Kin,col="#00FFFF")
myline(genKscoremeans$gen,genKscoremeans$KinQ,genKscoremins$KinQ,genKscoremaxs$KinQ,col="#0000FF")
myline(genKscoremeans$gen,genKscoremeans$HetQ,genKscoremins$HetQ,genKscoremaxs$HetQ,col="#FF00FF")
mtext("h) Introgression Age (Selecting PK)",adj=0.6,line=1,cex=cex.mtext)
dev.off()
