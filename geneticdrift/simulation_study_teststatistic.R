## Simulation study to test local ancestry anomaly scores via a detailed simulation
## This is a simulation in support of " Genetic swamping of the critically endangered Scottish wildcat was recent and accelerated by disease" by Howard-McCombe, et al
## Author:  Daniel Lawson (dan.lawson@bristol.ac.uk) 2023
## Licence: GPLv3

## Simulation functions: See https://doi.org/10.1101/2023.07.27.550812
## 
source("../localancestryfns.R")
### Start of real data processing

## Obtaining the empirical results from mosaic. External users will not have these so we save the histograms, which is all you need for the analysis presented.
hide <-function(){
    tmp=load("../convertmosaic/localanc_scot_2way_1-36_1-19_130_60_0.99_58.RData")
    ## [1] "localanc"    "final.flips" "g.loc"     
    tmp2=load("../convertmosaic/scot_2way_1-36_1-19_130_60_0.99_58.RData") 
    ## [1] "target"  "logfile" "Mu"      "lambda"  "theta"   "alpha"   "PI"     
    ## [8] "rho"     "A"       "NUMA"    "nchrno"  "chrnos"  "dr"      "NL"     
    ##[15] "kLL"     "acoancs" "all_Fst" "GpcM"   
    alphamosaic=sapply(alpha,function(x)x[2]) # alpha, according to mosaic (domestic,wildcat) pairs
    alpha=mean(alphamosaic) # [1] 0.579633
    ## Processing the empirical results into histogram bins
    localancestrychr5=localanc[[5]][1,,]
    catsmeanlist=lapply(localanc,function(x)apply(x[1,,],2,mean))
    catsmean=do.call("c",catsmeanlist[-5])
    meanchr5=apply(localancestrychr5,2,mean)
    write.table(catsmean,"wildcat_catsmean_exclchr5.txt",quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)
    write.table(meanchr5,"wildcat_chr5mean.txt",quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)
}
## Obtain the data above
alpha=0.579633
catsmean=read.table("wildcat_catsmean_exclchr5.txt")[,1]
meanchr5=read.table("wildcat_chr5mean.txt")[,1]
## Compute histograms
thcats=hist((1-catsmean)*72,breaks=seq(-0.5,72.5,by=1),plot=FALSE)
th5=hist((1-meanchr5)*72,breaks=seq(-0.5,72.5,by=1),plot=FALSE)

## 95% central quantile range
maskrange5=quantile(1-meanchr5,c(0.025,0.975))*72
maskrangecats=quantile(1-catsmean,c(0.025,0.975))*72


########## Simulation setup

## Simulation parameters
## These are chosen as the best-fit parameters from the scan below
N=400 # Individuals
G=20 # Number of gens after introgression starts
L=300 # Amount of genome, in SNPs. We assume 1 SNP/cM = 100 SNP/ Morgan for a basic sim. Doesn't affect the details of the fit, just increases the amount of averaging.
# 1-q/2 = alpha, q=2-2 alpha
infracseq=c(rep(2-2*alpha,1),rep(0,G-1)) # Introgression sequence
# Irrelevant parameters
NI=10 # Introgressing individuals (don't matter as we don't use SNPs)
## Construct the inputs to the simulation, data matrices of ancestry
set.seed(2)
data=list(panel=matrix(1,N*2,L))
dataIntrogress=list(panel=matrix(0,NI*2,L))
## Run the simulation for G generations with the specified introgression sequence
data0=breedingIntrogress(data,dataIntrogress,G,infracseq,usesex=FALSE)
# Extract data
panel=data0$panel
c(alpha=alpha,meansim=mean(panel)) # Check for the mean (it can drift)
## Extract 50 samples of 72 individuals, to match the real data
## Then sum the total ancestry in 72 individuals
cssimcats=sapply(1:50,function(x){
    testdata=panel[sort(sample(1:dim(panel)[1],72)),]
    colSums(testdata)
})
## Make histograme of each fo the 50 samples
thsimcatslist=lapply(1:50,function(i){
    hist(cssimcats[,i],breaks=seq(-0.5,72.5,by=1),plot=FALSE)})
## Annotate these with the 95% range
thsimcatslist=lapply(thsimcatslist,function(x){
    cs=cumsum(x$density)    
    x$maskrange=x$mids[c(max(which(cs<0.025)),min(which(cs>0.975)))]
    x$mask=as.numeric((x$mids>=x$maskrange[1])&(x$mids<=x$maskrange[2]))
    x
})
## Compute the Earth-mover distance and sort by increasing distance from the Chr5 results
thsimcatsemd=sapply(thsimcatslist,function(x)emdL1(x$density,th5$density))
thsimcatslist=thsimcatslist[order(thsimcatsemd)]
thsimcatsemd=thsimcatsemd[order(thsimcatsemd)]
 

########
## Run a range of simulations with different N and G
set.seed(2)
nreps=100
Nseq=c(400,600,800,1000,1200)
Gseq=c(5,10,15,20)
thsimlistind=list()
for(Non in 1:length(Nseq)){
    for(Gon in 1:length(Gseq)){
        G=Gseq[Gon]
        N=Nseq[Non]
        print(paste("Processing N =",N,"(",Non,"of",length(Nseq),")",
              "G =",G,"(",Gon,"of",length(Gseq),")"
              ))
        thsimlistind[[Gon+(Non-1)*length(Gseq)]]=list()
        for(rep in 1:nreps){
            infracseq=c(rep(2-2*alpha,1),rep(0,G-1)) # Introgression sequence
            data=list(panel=matrix(1,N*2,L))
            dataIntrogress=list(panel=matrix(0,NI*2,L))
            ## Generate the initial panel
            data0=breedingIntrogress(data,dataIntrogress,G,infracseq,usesex=FALSE)
            panel=data0$panel
            ## Sample 50 panels of 72 individuals
            cssimcats=sapply(1:50,function(x){
                testdata=panel[sort(sample(1:dim(panel)[1],72)),]
                colSums(testdata)
            })
            ## Here we use the results from all 50 runs to obtain an average
            th=hist(cssimcats,breaks=seq(-0.5,72.5,by=1),plot=FALSE)
            th$emd=emdL1(th$density,thcats$density)
            th$rmse=sum((th$density-thcats$density)^2)
            th$name=paste("N =",N,"G =",G)
            thsimlistind[[Gon+(Non-1)*length(Gseq)]][[rep]]=th
        }
    }
}

######## Average and annotate the histograms by their 95% central quantile range
thsimlist=lapply(thsimlistind,function(x){
    ret=x[[1]]
    ret$counts=rowMeans(sapply(x,function(y)y$counts))
    ret$density=rowMeans(sapply(x,function(y)y$density))
    cs=cumsum(ret$density)    
    ret$maskrange=ret$mids[c(max(which(cs<0.025)),min(which(cs>0.975)))]
    ret$mask=as.numeric((ret$mids>=ret$maskrange[1])&(ret$mids<=ret$maskrange[2]))
    ret
})

### Compute various measures of distance from the GENOME WIDE AVERAGE
rmse=sqrt(sapply(thsimlist,function(x)sum((x$density-thcats$density)^2)))
abse=sapply(thsimlist,function(x)sum(abs(x$density-thcats$density)))
emd=sapply(thsimlist,function(x)emdL1(x$density,thcats$density))
rmsetext=format(rmse,digits=2)
absetext=format(abse,digits=3)
emdtext=format(emd,digits=3)

## Plot the results when varying parameters
#pdf("SimulationModelCheck.pdf",height= 12,width=16)
png("SimulationModelCheck.png",height= 12*72,width=16*72)
par(mar=c(3.5,3.5,2,1),cex=1)
layout(matrix(1:length(thsimlist),nrow=length(Nseq),ncol=length(Tseq),byrow=T))
for(i in 1:length(thsimlist)){
    plot(thsimlist[[i]]$mids,thsimlist[[i]]$density,ylim=c(0,0.08),col="#0000FF55",
         main="",xlab="",ylab="",type="h",log="",lwd=2,cex.axis=1.5)
    mtext("Count of Wildcat alleles",side=1,line=2,cex=1)
    mtext("Density",side=2,cex=1,line=2)
#    mtext(paste0(letters[i],") ",thsimlist[[i]]$name," RMSE = ",rmsetext[i]),cex=0.75)
#    mtext(paste0(letters[i],") ",thsimlist[[i]]$name," |Err| = ",absetext[i]),cex=0.75)
    mtext(paste0(letters[i],") ",thsimlist[[i]]$name," EMD = ",emdtext[i]),cex=1)
    lines(thcats$mids-1/3,thcats$density,col="#00000099",type="h",lwd=2)
    lines(th5$mids+1/3,th5$density,col="#FF000055",type="h",lwd=2)
    lines(thsimlist[[i]]$mids,thsimlist[[i]]$density,col="#0000FFFF")
    lines(thcats$mids,thcats$density,col="#000000FF")
    lines(th5$mids,th5$density,col="#FF0000FF")
    abline(v=thsimlist[[i]]$maskrange,lty=2,col="blue")
    abline(v=maskrange5,lty=2,col="red")
    abline(v=maskrangecats,lty=2,col="black")
    if(i==1){
        legend("topleft",
               legend=c("Genome wide","Chr 5 (inc MHC)","Simulation"),
               lty=1,
               col=c("black","red","blue"),
               text.col=c("black","red","blue"),bg="white")
    }
}
dev.off()

## Plot the results for the best parameters
pdf("SimulationModelCheckFixedPars.pdf",height= 12,width=16)
par(mar=c(3,3,2,1))
layout(matrix(1:length(thsimlist),nrow=length(Tseq),ncol=length(Nseq)))
for(i in 1:length(thsimlist)){
    plot(thsimcatslist[[i]]$mids,thsimcatslist[[i]]$density,ylim=c(0,0.08),col="#0000FF55",
         main="",xlab="",ylab="",type="h",log="",lwd=3)
    mtext("Count of Wildcat alleles",side=1,line=2,cex=0.75)
    mtext("Density",side=2,cex=0.75,line=2)
#    mtext(paste0(letters[i],") ",thsimlist[[i]]$name," RMSE = ",rmsetext[i]),cex=0.75)
#    mtext(paste0(letters[i],") ",thsimlist[[i]]$name," |Err| = ",absetext[i]),cex=0.75)
    mtext(paste0(letters[i],") EMD (Sim,Chr5) = ",format(thsimcatsemd[i],digits=3)),cex=0.75)
    lines(thcats$mids,thcats$density,col="#00000099",type="h",lwd=3)
    lines(th5$mids,th5$density,col="#FF000055",type="h",lwd=3)
    lines(thsimcatslist[[i]]$mids,thsimcatslist[[i]]$density,col="#0000FFFF")
    lines(thcats$mids,thcats$density,col="#000000FF")
    lines(th5$mids,th5$density,col="#FF0000FF")
    abline(v=thsimcatslist[[i]]$maskrange,lty=2,col="blue")
    abline(v=maskrange5,lty=2,col="red")
    abline(v=maskrangecats,lty=2,col="black")
    if(i==1){
        legend("topleft",
               legend=c("Genome wide","Chr 5 (inc MHC)","Simulation"),
               lty=1,
               col=c("black","red","blue"),
               text.col=c("black","red","blue"),bg="white")
    }
}
dev.off()
