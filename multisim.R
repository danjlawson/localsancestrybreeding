



########################################
## genetic distance scan

## genetic distance
nrep=3
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

## Construct score
gLscoreres=data.frame(t(sapply(1:length(gLMHQScore),function(i){
    gL=gLtest[i]
    score=tail(gLMHQScore[[i]]$score,1)
    score$gL=gL
    as.numeric(score)
})))
colnames(gLscoreres)=c(colnames(gLMHQScore[[i]]$score),"gL")

gLscoremeans=aggregate.data.frame(gLscoreres,by=list(groupgL=gLscoreres$gL),mean)
gLscoremins=aggregate.data.frame(gLscoreres,by=list(groupgL=gLscoreres$gL),min)
gLscoremaxs=aggregate.data.frame(gLscoreres,by=list(groupgL=gLscoreres$gL),max)


########################################
## admixture scan

## admixture sims
atest=rep(seq(0.05,0.5,length.out=10),each=3)
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

## Construct score
alphascoreres=data.frame(t(sapply(1:length(aMHQScore),function(i){
    alpha=atest[i]
    score=tail(aMHQScore[[i]]$score,1)
    score$alpha=alpha
    as.numeric(score)
})))
colnames(alphascoreres)=c(colnames(aMHQScore[[i]]$score),"alpha")

alphaxvals=aggregate(1-aalpha0,list(groupalpha=alphascoreres$alpha),mean)
alphascoreres$NonQ0=1-aalpha0
alphascoremeans=aggregate.data.frame(alphascoreres,
                                     by=list(groupalpha=alphascoreres$alpha),mean)
alphascoremins=aggregate.data.frame(alphascoreres,
                                    by=list(groupalpha=alphascoreres$alpha),min)
alphascoremaxs=aggregate.data.frame(alphascoreres,
                                    by=list(groupalpha=alphascoreres$alpha),max)
to=order(alphascoremeans$NonQ0)
alphascoremeans=alphascoremeans[to,]
alphascoremins=alphascoremins[to,]
alphascoremaxs=alphascoremaxs[to,]

############################################
## plot
myline<-function(x,ymean,ymin,ymax,col,...){
    pcol=paste0(col,"55")
    polygon(c(x,rev(x)),c(ymin,rev(ymax)),border=NA,col=pcol,...)
    lines(x,ymean,col=col,lwd=2,...)
}

########################
## Figure 4
pdf("Figure4-ChangingSim.pdf",height=4,width=10)
par(las=1,xpd=NA)
layout(matrix(c(1,3,2),nrow=1),widths=c(1,0.7,1))
plot(alphascoremeans$NonQ0,alphascoremeans$NonQ,type="n",col=2,ylim=c(0,0.5),xlab="Genome Length (Morgans)",ylab="Score after 20 Generations")
myline(alphascoremeans$NonQ0,alphascoremeans$NonQ,alphascoremins$NonQ,alphascoremaxs$NonQ,col="#FF0000")
myline(alphascoremeans$NonQ0,alphascoremeans$Het,alphascoremins$Het,alphascoremaxs$Het,col="#00FF00")
myline(alphascoremeans$NonQ0,alphascoremeans$Kin,alphascoremins$Kin,alphascoremaxs$Kin,col="#00FFFF")
myline(alphascoremeans$NonQ0,alphascoremeans$HetQ,alphascoremins$HetQ,alphascoremaxs$HetQ,col="#0000FF")
myline(alphascoremeans$NonQ0,alphascoremeans$KinQ,alphascoremins$KinQ,alphascoremaxs$KinQ,col="#FF00FF")
mtext("a) Impact of Introgression Magnitude",adj=-0.2,line=1,cex=0.8)
##
plot(gLscoremeans$gL,gLscoremeans$NonQ,type="n",col=2,ylim=c(0,0.5),xlab="Genome Length (Morgans)",ylab="Score after 20 Generations")
myline(gLscoremeans$gL,gLscoremeans$NonQ,gLscoremins$NonQ,gLscoremaxs$NonQ,col="#FF0000")
myline(gLscoremeans$gL,gLscoremeans$Het,gLscoremins$Het,gLscoremaxs$Het,col="#00FF00")
myline(gLscoremeans$gL,gLscoremeans$Kin,gLscoremins$Kin,gLscoremaxs$Kin,col="#00FFFF")
myline(gLscoremeans$gL,gLscoremeans$HetQ,gLscoremins$HetQ,gLscoremaxs$HetQ,col="#0000FF")
myline(gLscoremeans$gL,gLscoremeans$KinQ,gLscoremins$KinQ,gLscoremaxs$KinQ,col="#FF00FF")
mtext("b) Impact of Genome Length",adj=-0.2,line=1,cex=0.8)
##
par(mar=c(0,0,0,0))
plot(c(0,1),c(0,1), type="n",xlab="",ylab="",axes=FALSE)
legend("center",legend=c("Introgression","Heterozygosity","Kinship",
                          "Subpopulation Heterozygosity SH",
                         "Subpopulation Kinship SK"),
       text.col=c(2,3,5,4,6),lty=1,col=c(2,3,5,4,6),title="Score:",title.col=1,cex=1.2,lwd=2)
dev.off()
