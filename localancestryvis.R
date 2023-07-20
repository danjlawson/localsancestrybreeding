source("localancestryfns.R")

visdata<-simCaptive(100,50,500,gL=5,seed=1,usesex=FALSE)

Gforward=20
#dataPairwiseRandom<-breedingSel(data,Gforward,selfn=pairwiseRandom,parentfn=rankedParents)
visdataPairwiseQ<-breedingSel(visdata,Gforward,selfn=pairwiseQ,parentfn=rankedParents,min=FALSE)
visdataPairwiseMK<-breedingSel(visdata,Gforward,selfn=pairwiseKinship,parentfn=rankedParents,min=TRUE)
#dataPairwiseMH<-breedingSel(data,Gforward,selfn=pairwiseHeterozygosity,parentfn=rankedParents,min=FALSE)
#dataPairwiseMHQ<-breedingSel(data,Gforward,selfn=pairwiseHeterozygosityQ,parentfn=rankedParents,min=FALSE)
visdataPairwiseMHQScore<-breedingSel(visdata,Gforward,selfn=pairwiseHeterozygosityQScore,parentfn=rankedParents,min=FALSE)
#dataPairwiseMKQ<-breedingSel(data,Gforward,selfn=pairwiseKinshipQ,parentfn=rankedParents,min=TRUE)
#dataPairwiseMKQScore<-breedingSel(data,Gforward,selfn=pairwiseKinshipQScore,parentfn=rankedParents,min=TRUE)

myimage<-function(data){
    esh=getExpectedSubpopulationHetScore(data)*2
    tdend=as.dendrogram(hclust(dist(esh)))
    tdend=reorder(tdend,rowMeans(esh))
    torder=as.numeric(labels(tdend))
    esh=esh[torder,]
    image(t(esh),col=c("white","cyan","orange","blue","red"))
}

myimage2<-function(data,inds=rep(TRUE,dim(data[[1]])[1]/2),snps=rep(TRUE,dim(data[[1]])[2])){
    ## tdend=as.dendrogram(hclust(dist(data$genomes)))
    ## tdend=reorder(tdend,rowMeans(data$genomes))
    ## torder=as.numeric(labels(tdend))
    ## data$panel=data$panel[torder,]
    ## data$genomes=data$genomes[torder,]
    ## data$founder=data$founder[torder,]
    c1=(data$panel)&(data$genomes)
    c2=(data$panel)&(!data$genomes)
    c3=(!data$panel)&(data$genomes)
    c4=(!data$panel)&(!data$genomes)
    pp= c3 + 2*c1+ 3*c2
    L=sum(snps)
    N=sum(inds)
    image(0:L,0:(2*N),t(pp[rep(inds,each=2),snps]),
          col=c("white","cyan","orange","orange","red"),zlim=c(0,3),
          xlab="",ylab="",axes=FALSE)
    for(i in 1:N){
        rect(0,2*(i-1),L,2*i,border="white")
    }
    invisible(pp[rep(inds,each=2),snps])
}
mergemitosis<-function(p1m,p2m,idx){
    ret=list()
    for(p in 1:length(p1m[[1]])){
        ret[[p]]=rbind(p1m[[idx[1]]][[p]],
                       p2m[[idx[1]]][[p]])
    }
    for(i in idx[-1]){
        for(p in 1:length(p1m[[1]])){
            ret[[p]]=rbind(ret[[p]],
                           p1m[[i]][[p]],
                           p2m[[i]][[p]])
        }
    }
    names(ret)=c("panel","genomes","founder")
    ret
}
applyinds<-function(data,inds){
    for(i in 1:length(data)) data[[i]]=data[[i]][rep(inds,each=2),]
    data
}
applysnps<-function(data,snps){
    for(i in 1:length(data)) data[[i]]=data[[i]][,snps]
    data
}
    
ttab=visdata$panel*(visdata$genome+1)
snps=apply(ttab,2,function(x){ sum(x==2)*sum(x==1) })>1000
snps[1:400]=FALSE
inds=rep(FALSE,100)
set.seed(1)
inds[sample(1:100,8)]<-TRUE
#inds[1:8]<-TRUE

png("selecionpanels.png",height=800,width=1600)
layout(matrix(c(5,6,1,2,3,4),ncol=2))
par(mar=c(1,6,6,6))
pp=myimage2(visdata,inds=inds,snps=snps)
ppQ=myimage2(visdataPairwiseQ$data,inds=inds,snps=snps)
ppK=myimage2(visdataPairwiseMK$data,inds=inds,snps=snps)
ppHQ=myimage2(visdataPairwiseMHQScore$data,inds=inds,snps=snps)
dev.off()

visdatazoomed=applysnps(visdata,snps)
p1=which(inds==1)[1]; h1=2*p1+c(-1,0)
p2=which(inds==1)[6]; h2=2*p2+c(-1,0)
indsel=rep(FALSE,100)
indsel[c(p1,p2)]=TRUE
pardata=applyinds(visdatazoomed,indsel)
for(i in 1:3)for(rep in 1:2) pardata[[i]]<-rbind(NA,pardata[[i]],NA)
p1data=getinddata(visdatazoomed,p1)
p2data=getinddata(visdatazoomed,p2)
set.seed(1)
p1mitosis=lapply(1:1000,function(i)crossover(p1data))
p2mitosis=lapply(1:1000,function(i)crossover(p2data))
p1mitosis=p1mitosis[order(sapply(p1mitosis,function(x)mean(x[[1]])),decreasing=FALSE)]
p2mitosis=p2mitosis[order(sapply(p2mitosis,function(x)mean(x[[1]])),decreasing=FALSE)]

newoffspring=mergemitosis(p1mitosis,p2mitosis,c(1,250,750,1000))

png("mitosis.png",height=400,width=800)
par(mfrow=c(1,2),mar=c(1,3,1,3))
myimage2(pardata)
myimage2(newoffspring)
dev.off()

pp[11:12,]
pp[15:16,]

image(1-t(dataPairwiseMK$data$panel))
image(1-t(dataPairwiseMHQScore$data$panel))


png("selecionpanels2.png",height=800,width=1600)
layout(matrix(c(5,6,1,2,3,4),ncol=2))
par(mar=c(1,6,6,6))
pp=myimage2(visdata,inds=inds,snps=snps)
ppQ=myimage2(visdataPairwiseQ$data,inds=inds,snps=snps)
ppK=myimage2(visdataPairwiseMK$data,inds=inds,snps=snps)
ppHQ=myimage2(visdataPairwiseMHQScore$data,inds=inds,snps=snps)
dev.off()

load("sim.RData")
scores=c("NonQ","KinQ")
names=c("Introgressed Admixture","Subpop Kinship")
symbols=c("Domestic Cat proportion","Local Ancestry Kinship")
names(symbols)=names(names)=scores

png("simplifiedbreedingcurve.png",height=800,width=1600)
mylwd=6
par(mfrow=c(1,2),las=1,cex=3)
for(scoreon in 1:length(scores)){
    score=scores[scoreon]
    panel=c("e","f")[scoreon]
    yrange=c(0,0.5)
    if(score=="KinQ") yrange=c(0,1)
    plot(range(dataNull$score$g),yrange,xlab="Generation",ylab=symbols[score],
         type="n",main="",cex.axis=1.5,cex.lab=1.2)
#    mtext(paste0(panel,") Evolution of ",names[score]),adj=-0.2,line=1,cex=4)
    lines(dataPairwiseRandom$score[,"g"],dataPairwiseRandom$score[,score],col="grey",
          lwd=mylwd,cex=0.75,lty=3)
    lines(dataPairwiseQ$score[,"g"],dataPairwiseQ$score[,score],col="grey",lwd=mylwd,lty=3)
    lines(dataPairwiseMK$score[,"g"],dataPairwiseMK$score[,score],col="grey",lwd=mylwd,lty=3)
    lines(dataPairwiseMKQ$score[,"g"],dataPairwiseMKQ$score[,score],col="red",lty=1,lwd=mylwd)
    if(score=="NonQ") legend("topleft",legend=c("Traditional options","Local ancestry"),
                             text.col=c("black","red"),lty=1,col=c("grey","red"),
                             title="Breeding using:",cex=1,lwd=mylwd)
}
dev.off()
