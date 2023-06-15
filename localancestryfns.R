## Simulation for breeding and introgression where:
## Genetic and local ancestry data are stored as a list(ancestry,genetic)

simCaptive<-function(N,NI,L,gL=L/100,Fst=0.2,G=10,Galpha=2,alpha=0.3,seed=NA,usesex=TRUE){
    ## Simulate a captive population
    ## N= Individuals
    ## NI = Introgressing individuals
    ## L = Amount of genome, in SNPs. We assume separated by 1cM for a basic sim
    ## This is a stand-alone script that should work regardless of data sources
    ## Fst=0.2 # Fst between MRCA of two populations
    ## G=10 # Number of gens after introgression starts
    ## Galpha=2 # Number of gens of introgression
    ## alpha=0.3 # Introgression alpha during Galpha period
    ## usesex: whether to model sex
    ## seed =NA random seed for replicability
    if(!is.na(seed)) set.seed(seed)
    f0=runif(L,0.05,0.95) # Frequencies in the MRCA
    f=simf(f0,Fst) # Frequencies in the target population
    fintrogress=simf(f0,Fst) # Frequencies in the introgressing population
    gdist=diff(seq(0,gL,length.out=L+1))
    infracseq=c(rep(alpha,Galpha),rep(0,G-Galpha)) # Introgression sequence

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
    data0=breedingIntrogress(data,dataIntrogress,G,infracseq,d=gdist,usesex=usesex)
    data0$founder=matrix(rep(1:N,each=2),N*2,L) ## Re-founder the panel
    ## Generate the captive population
    data<-breedingSel(data0,1,selfn=selRandom,frac=0.25,d=gdist,usesex=usesex)$data ## Random mating bottleneck
    data<-breedingSel(data,1,selfn=selMaxQ,frac=0.25,d=gdist,usesex=usesex)$data ## Removal of highly introgressed individuals
    data
}

crossover=function(data,d=rep(0.01,dim(data[[1]])[2])){
    ## Crosses over data for two individuals, stored as rows of ancestry: data[[1]] and genome: data[[2]]
    cdistance=cumsum(d)
    distance=sum(d)
    ploidy=dim(data[[1]])[1] ## must be 2 for the code to make any sense!
    if(ploidy!=2) stop("ploidy must be 2")
    p=0 # position
    i=1 # recombination index
    rec=numeric() # recombination end points
    while((i<length(d)) &&(p<distance)){ 
        p=p+rexp(1,1) # next recombination location
        while((i<length(d)) &&(p>cdistance[i])){
            i=i+1
        }
        rec=c(rec,i) # report segment
    }
    rec

    ret=list()
    for(v in 1:length(data)){
        ret[[v]]=rep(NA,dim(data[[v]])[2]) # sequence we return
    }
    curr=sample(1:ploidy,1) # current sample
    i=1 ## current start index
    for(j in 1:length(rec)){ ## loop over recombination end points
        for(v in 1:length(data)) ret[[v]][i:rec[j]]=data[[v]][curr,i:rec[j]] ## update current
        if(ploidy==2){curr=3-curr ## to make it polyploid: make it random from the others
        }else if(ploidy>2) {curr=sample((1:ploidy)[-curr],1)
        }else curr=1
        i=rec[j]+1
    }
    
    ret
}

### Introgression step
getinddata<-function(data,i){
    ## Extract panel and genome data
    ret<-list()
    for(v in 1:length(data)){
        ret[[v]]<-data[[v]][(-1:0)+2*i,]
    }
    ret
}
combineData<-function(x,y){
    ## Combine two datasets stored in list format
    for(v in 1:length(x)){
        x[[v]]<-rbind(x[[v]],y[[v]])
    }
    x
}
sampleParentsIntrogress<-function(data,dataIntrogress,infrac,usesex=TRUE){
    ## Sample parents allowing introgression
    N=dim(data[[1]])[1]/2
    NI=dim(dataIntrogress[[1]])[1]/2
    if(!usesex){
        parents=data.frame(
            p1=sample(1:N,N,TRUE),
            p2=sample(1:N,N,TRUE),
            p3=sample(1:NI,N,TRUE),
            introgressed=sample(0:1,N,TRUE,prob=c(1-infrac,infrac)))
    }else{
        parents=data.frame(
            p1=sample((1:N)[attr(data,"sexf")==1],N,TRUE),
            p2=sample((1:N)[attr(data,"sexf")==0],N,TRUE),
            p3=sample((1:NI)[attr(dataIntrogress,"sexf")==0],N,TRUE),
            introgressed=sample(0:1,N,TRUE,prob=c(1-infrac,infrac)))
    }
    return(parents)
}

breedingIntrogress<-function(data,dataIntrogress,G,infracseq=rep(0,G),d=rep(0.01,dim(data[[1]])[2]),usesex=TRUE){
    ## Breed at random, with a particular introgression profile, from a fixed panel of introgressing individuals
    lastdata=data
    N=dim(data[[1]])[1]/2
    for(g in 1:G){
        infrac=infracseq[g]
        parents=sampleParentsIntrogress(data,dataIntrogress,infrac,usesex=usesex)
        for(i in 1:N){
            tp1=crossover(getinddata(lastdata,parents[i,1]),d=d)
            if(parents[i,"introgressed"]==1){
                tp2=crossover(getinddata(dataIntrogress,parents[i,3]),d=d)
            }else tp2=crossover(getinddata(lastdata,parents[i,2]),d=d)
            for(v in 1:length(data)){
                data[[v]][2*i-1,]=tp1[[v]]
                data[[v]][2*i,]=tp2[[v]]
            }            
        }
        sexok=FALSE ## Forbid all one gender for technical reasons
        while(!sexok){
            attr(data,"sexf")=rbinom(N,1,0.5)
            if(sd(attr(data,"sexf"))>0) sexok=TRUE
        }
        lastdata=data
    }
    data
}

simf<-function(f0,Fst){
    ##    https://www.nature.com/articles/s41431-021-00873-2
    ## Frequency change according to genetic drift by an amount equal to Fst
    rbeta(length(f0),f0*(1-Fst)/Fst,(1-f0)*(1-Fst)/Fst)
}

sampleParents<-function(data,selfn,usesex=FALSE,min=TRUE,...){
    ## Sample parents with a simple weighting function
    N=dim(data[[1]])[1]/2
    score=selfn(data,...) ## MOVE INTO parentfn :(
    if(!usesex){
        ## Implement their algorithm without the constraint of sex
        parents=data.frame(p1=sample(1:N,N,replace=TRUE,prob=score),
                           p2=sample(1:N,N,replace=TRUE,prob=score))
    }else{
        parents=data.frame(
            p1=sample((1:N)[attr(data,"sexf")==1],N,TRUE,prob=score[attr(data,"sexf")==1]),
            p2=sample((1:N)[attr(data,"sexf")==0],N,TRUE,prob=score[attr(data,"sexf")==0]))
    }
    return(parents)
}
rankedParents<-function(data,selfn=pairwiseRandom,min=TRUE,usesex=FALSE,offspringdist=function(N)rpois(N,4.2),...){
    ## Sample parents with a ranked weighting procedure
    ## https://academic.oup.com/jhered/article/103/2/186/886373?login=false#81041438
    ## This is Ranked MK Selection if we provide selfn=pairwiseKinship
    ## selfn MUST be a pairwise quantity
    N=dim(data[[1]])[1]/2
    scoremat=selfn(data,...) ## MOVE INTO parentfn :(
    if(!usesex){
        ## Implement their algorithm without the constraint of sex
        presentlist=1:N
        score=rowMeans(scoremat)
        parentlist<-numeric()
        for(i in 1:N){ ## Remove the worst candidate one by one, computing the mean of the scorematrix for the remaining
            scoreorder=order(score,decreasing=min)
            who=presentlist[scoreorder[1]]
            parentlist<-c(parentlist,who)
            scoremat=scoremat[-scoreorder[1],-scoreorder[1],drop=FALSE]
            score=rowMeans(scoremat)
            presentlist=presentlist[-scoreorder[1]]
        }
        nchildren<-offspringdist(N)
        parents=data.frame(p1=numeric(),numeric())
        i=1
        parentlist=rev(parentlist)
        while(dim(parents)[1]<N){
            parents<-rbind(parents,
                           data.frame(p1=rep(parentlist[2*i-1],nchildren[i]),
                                      p2=rep(parentlist[2*i],nchildren[i])))
            i=i+1
        }
        parents=parents[1:N,]
    }else{
        ## Urgh its a faff INCOMPLETE
        ## Implement their algorithm as described
        presentlistf=(1:N)[attr(data,"sexf")==1]
        presentlistm=(1:N)[attr(data,"sexf")==0]
        score=rowMeans(scoremat)
        parentlistf<-numeric()
        parentlistm<-numeric()
        while(length(presentlistf)>0){
            scoreorder=order(score,decreasing=min)
            who=presentlist[scoreorder[1]]
            parentlist<-c(parentlist,who)
            scoremat=scoremat[-scoreorder[1],-scoreorder[1],drop=FALSE]
            score=rowMeans(scoremat)
            presentlist=presentlist[-scoreorder[1]]
        }
        nchildren<-offspringdist(N)
        parents=data.frame(p1=numeric(),numeric())
        i=1
        while(dim(parents)[1]<N){
            parents<-rbind(parents,
                           data.frame(p1=rep(parentlist[2*i-1],nchildren[i]),
                                      p2=rep(parentlist[2*i],nchildren[i])))
            i=i+1
        }
        parents=parents[1:N,]
        parents=data.frame(
            p1=sample((1:N)[attr(data,"sexf")==1],N,TRUE,prob=wt[attr(data,"sexf")==1]),
            p2=sample((1:N)[attr(data,"sexf")==0],N,TRUE,prob=wt[attr(data,"sexf")==0]))
    }
    return(parents)
}

breedingSel<-function(data,G,scorefn=scoreMeanQ,
                      selfn=selIdentity,d=rep(0.01,dim(data[[1]])[2]),
                      parentfn=sampleParents, usesex=FALSE,verbose=TRUE,...){
    ## Breed for G gens, applying a specified selection procedure for choosing parents
    lastdata=data
    score=allScores(0,data)
    N=dim(data[[1]])[1]/2
    for(g in 1:G){
        if(verbose) print(tail(score,1))
        parents=parentfn(data,selfn,usesex=usesex,...)
        for(i in 1:N){
            tp1=crossover(getinddata(lastdata,parents[i,1]),d=d)
            tp2=crossover(getinddata(lastdata,parents[i,2]),d=d)
            for(v in 1:length(data)){
                data[[v]][2*i-1,]=tp1[[v]]
                data[[v]][2*i,]=tp2[[v]]
            }
        }
        sexok=FALSE ## Forbid all one gender for technical reasons
        while(!sexok){
            attr(data,"sexf")=rbinom(N,1,0.5)
            if(sd(attr(data,"sexf"))>0) sexok=TRUE
        }
        score=rbind(score,allScores(g,data))
        lastdata=data
    }
    list(score=score,data=data)
}
getQ<-function(x){
    ## Average admixture fraction per individual
    Qm=matrix(rowMeans(x[[1]]),nrow=2)
    Q=colMeans(Qm)
}
getHet<-function(x){
    ## Return an N by L matrix reporting if each site is a het
    g=x[[2]]
    N<-dim(g)[1]/2
    L<-dim(g)[2]
    het<-matrix(0,N,L)
    for(i in 1:N){
        het[i,]<-g[2*i-1,]!=g[2*i,]
    }
    het
}
getExpectedSubpopulationHetScore<-function(x,delta=0.5){
    ## Return an N by L matrix reporting if each site is a
    a=x[[1]] ## ancestry
    g=x[[2]] ## genotype (0,1)
    N<-dim(g)[1]/2
    L<-dim(g)[2]
    ret<-matrix(0,N,L)
    for(i in 1:N){
        ## score 0 for no wildcat, 0.5 for 1 wildcat hap, 1 for 2 wildcat haps with same variant, 2 for 2 wildcat variants
        ret[i,]<-(1-delta)*(a[2*i,] + a[2*i-1,]) + 2*delta* ((a[2*i,] + a[2*i-1,]) ==2)*abs(g[2*i,]-g[2*i-1,])
    }
    ret
}
getExpectedKinshipCoeficient<-function(x,fthresh=0.05){
    ## Return an N by L matrix reporting if each site is a
    g=x[[2]] ## genotype (0,1)
    N<-dim(g)[1]/2
    L<-dim(g)[2]
    f=apply(g,2,mean)
    f[f<fthresh]=fthresh
    f[f>1-fthresh]=1-fthresh
    ret<-matrix(0,N,L)
    for(i in 1:N){
        ##        ret[i,]<-(1 - abs(g[2*i,]==g[2*i-1,]))/f/(1-f)/2
        ret[i,]<-1-(g[2*i,]!=g[2*i-1,])/f/(1-f)/2
    }
    ret
}
getKinship<-function(x){
    ## Empirical kinship from founder
    k=x[[3]]
    N<-dim(k)[1]/2
    L<-dim(k)[2]
    ret<-matrix(0,N,L)
    for(i in 1:N){
        ret[i,]<- k[2*i,]==k[2*i-1,]
    }
    ret
}
getSubpopulationKinship<-function(x,fthresh=0.05,delta=0.5){
    ## Empirical kinship from founder
    a=x[[1]] ## ancestry
    k=x[[3]]
    N<-dim(k)[1]/2
    L<-dim(k)[2]
    ret<-matrix(0,N,L)
    for(i in 1:N){
        ret[i,]<- ((1-delta)*(2-a[2*i,] - a[2*i-1,])+
                   delta * (k[2*i,]==k[2*i-1,]))
    }
    ret    
}
getExpectedSubpopulationKinshipCoefficient<-function(x,fthresh=0.05,delta=0.5){
    ## Return an N by L matrix reporting if each site is a
    a=x[[1]] ## ancestry
    g=x[[2]] ## genotype (0,1)
    f=apply(g,2,mean)
    f[f<fthresh]=fthresh
    f[f>1-fthresh]=1-fthresh
    N<-dim(g)[1]/2
    L<-dim(g)[2]
    ret<-matrix(0,N,L)
    for(i in 1:N){
        ret[i,]<-1 - ((1-delta)*(a[2*i,] + a[2*i-1,])
            + 2*delta*((a[2*i,] + a[2*i-1,]) ==2)*(g[2*i,]!=g[2*i-1,])/f/(1-f)/2)
    }
    ret
}
getSnpQ<-function(x){
    ## Return an N by L matrix of the average introgression for each site for each individual
    d=x[[1]]
    N<-dim(d)[1]/2
    L<-dim(d)[2]
    snpQ<-matrix(0,N,L)
    for(i in 1:N){
        snpQ[i,]<-(d[2*i-1,]+d[2*i,])/2
    }
    snpQ
}

### Scores
scoreMeanKinship<-function(x){
    K=getKinship(x)
    mean(K)
}
scoreHet<-function(x){
    ## score is the average Heterozygosity
    h<-getHet(x)
    mean(h)
}
scoreHetQ<-function(x){
    ## score is the average Heterozygosity at wildcat sites
    Q=getSnpQ(x)
    h<-getHet(x)
    h[Q!=1]= 0  ## reward Q==1 (maximising)
#    wt<-(Q==1) ## reward Q==1 (maximising)
    mean(h)
}
scoreKin<-function(x){
    ## score is the average Kinship
    K<-getKinship(x)
    mean(K)
}
scoreKinQ<-function(x){
    ## score is the average Kinship at wildcat sites
    Q=getSnpQ(x)
    K<-getKinship(x)
    K[Q!=1]=1  ## reward Q==1 (minimising)
    mean(K)
}
scoreExpectedKin<-function(x){
    ## score is the average Kinship
    K<-getExpectedKinshipCoeficient(x)
    mean(K)
}
scoreExpectedKinQ<-function(x){
    ## score is the average Kinship at wildcat sites
    Q=getSnpQ(x)
    K<-getExpectedKinshipCoeficient(x)
    K[Q!=1] = 1  ## reward Q==1 (minimising)
    mean(K)
}
scoreMeanQ<-function(x){
    ## score is the average admixture
    mean(x[[1]])
}
scoreInfo<-function(x){
    p=colMeans(x[[2]])
    logp=log(p)
    logp[p==0]=0
    q=1-p
    logq=log(q)
    logq[q==0]=0
    -mean(p*logp+q*logq)
}
scoreAe<-function(x){
    ## expected heterozygosity
    p=colMeans(x[[2]])
    mean(1-p^2-(1-p)^2)
}
scoreInfoQ<-function(x){
    ## score is the average Heterozygosity at wildcat sites    
    x[[2]]=2*x[[2]]-1
    tx=x[[2]]*(1-x[[1]])
    nsites=colSums(abs(tx))
    nsites[nsites==0]=1
    np=colSums(tx==1)
    nq=colSums(tx==-1)
    p=np/nsites
    q=nq/nsites
    logp=log(p)
    logp[p==0]=0
    logq=log(q)
    logq[q==0]=0
    -mean(p*logp+q*logq)
}
scoreEffectiveLoci<-function(x){
    ## Effective number of loci, compared to all being frequency 0.5
    p=colMeans(x[[2]])
    q=1-p
    mean(p*q/0.25)
}
scoreEffectiveLociQ<-function(x){
    ## Effective number of wildcat loci, compared to all being frequency 0.5
    x[[2]]=2*x[[2]]-1
    tx=x[[2]]*(1-x[[1]])
    nsites=colSums(abs(tx))
    nsites[nsites==0]=1
    nsites=rep(dim(x[[1]])[1],length(nsites))
    np=colSums(tx==1)
    nq=colSums(tx==-1)
    p=np/nsites
    q=nq/nsites
    mean(p*q/0.25)
}
allScores<-function(g,x){
    ## Report all available scores
    r=data.frame(g=g,
                 Q=scoreMeanQ(x),
                 K=scoreMeanKinship(x),
                 Het=scoreHet(x),
                 HetQ=scoreHetQ(x),
                 Ae=scoreAe(x),
                 Info=scoreInfo(x),
                 InfoQ=scoreInfoQ(x),
                 EffectiveLociQ=scoreEffectiveLociQ(x),
                 EffectiveLoci=scoreEffectiveLoci(x),
                 Kin=scoreKin(x),
                 KinQ=scoreKinQ(x),
                 ExpectedKin=scoreExpectedKin(x),
                 ExpectedKinQ=scoreExpectedKinQ(x)
                 )
    r$NonQ=1-r$Q
    r
}
######################################################################
### Pairwise scores
pairwiseKinship<-function(x){
    ## Expected value of the Pairwise kinship probability
    ## where x is the data list and the result is an N by N matrix
    ## 
    X=x[[3]]
    N<-dim(X)[1]/2
    ret<-matrix(0,N,N)
    for(i in 1:(N)){
        for(j in (i):N){
            ret[j,i]<-ret[i,j]<- mean((X[2*i-1,]==X[2*j-1,])+(X[2*i-1,]==X[2*j,]) + (X[2*i,]==X[2*j-1,])+ (X[2*i,]==X[2*j,]))/4
        }
    }
    ret
}
pairwiseHeterozygosity<-function(x){
    ## Expected value of the Pairwise Heterozygosity
    ## where x is the data list and the result is an N by N matrix
    ## 
    X=x[[2]]
    N<-dim(X)[1]/2
    ret<-matrix(0,N,N)
    for(i in 1:(N)){
        for(j in (i):N){
            ret[j,i]<-ret[i,j]<- mean((X[2*i-1,]!=X[2*j-1,])+(X[2*i-1,]!=X[2*j,]) + (X[2*i,]!=X[2*j-1,])+ (X[2*i,]!=X[2*j,]))/4
        }
    }
    ret
}
pairwiseQ<-function(x){
    ## Expected value of the Pairwise Admixture
    ## where x is the data list and the result is an N by N matrix
    ## 
    Q=x[[1]]
    N<-dim(Q)[1]/2
    ret<-matrix(0,N,N)
    for(i in 1:(N)){
        for(j in (i):N){
            ret[j,i]<-ret[i,j]<- mean(Q[2*i-1,]+Q[2*i,]+Q[2*j-1,]+Q[2*j,])/4
        }
    }
    ret
}
pairwiseRandom<-function(x){
    ## Obtains a random ordering as a pairwise matrix with identical columns, for use in rankedParents
    N=dim(x[[1]])[1]/2
    torder=sample(1:N,N)
    matrix(torder,nrow=N,ncol=N)
}
pairwiseHeterozygosityQ<-function(x){
    ## Expected value of the Pairwise Heterozygosity for wildcat loci
    ## where x is the data list and the result is an N by N matrix
    ##
    Q=x[[1]]
    X=x[[2]]
    N<-dim(X)[1]/2
    ret<-matrix(0,N,N)
    for(i in 1:(N)){
        for(j in (i):N){
            ret[j,i]<-ret[i,j]<- mean(
            ((Q[2*i-1,]+Q[2*j-1,])==2)*(X[2*i-1,]!=X[2*j-1,])+
            ((Q[2*i-1,]+Q[2*j,])==2)  *(X[2*i-1,]!=X[2*j,]) +
            ((Q[2*i,]+Q[2*j-1,])==2)  *(X[2*i,]!=X[2*j-1,])+
            ((Q[2*i,]+Q[2*j,])==2)    *(X[2*i,]!=X[2*j,]))/4
        }
    }
    ret
}
pairwiseHeterozygosityQScore<-function(x,delta=0.5){
    ## Expected value of the Pairwise Heterozygosity
    ## where x is the data list and the result is an N by N matrix
    ##
    wt=getExpectedSubpopulationHetScore(x)
    wt2=1/colMeans(wt)
    wt2=wt2/mean(wt2)
    Q=x[[1]]
    X=x[[2]]
    N<-dim(X)[1]/2
    ret<-matrix(0,N,N)
    for(i in 1:(N)){
        for(j in (i):N){
            ret[j,i]<-ret[i,j]<- mean(
            (1-delta)*(
                Q[2*i-1,]+Q[2*j-1,]+
                Q[2*i-1,]+Q[2*j,]+
                Q[2*i,]+Q[2*j-1,]+
                Q[2*i,]+Q[2*j,] ) +
            2*delta*(
                ((Q[2*i-1,]+Q[2*j-1,])==2)*(X[2*i-1,]!=X[2*j-1,])+
                ((Q[2*i-1,]+Q[2*j,])==2)*(X[2*i-1,]!=X[2*j,]) +
                ((Q[2*i,]+Q[2*j-1,])==2)*(X[2*i,]!=X[2*j-1,])+
                ((Q[2*i,]+Q[2*j,])==2)*(X[2*i,]!=X[2*j,])
            ))/4
        }
    }
    ret
}
pairwiseKinshipQ<-function(x,delta=0.5){
    ## Expected value of the Pairwise Heterozygosity
    ## where x is the data list and the result is an N by N matrix
    ##
    Q=x[[1]]
    K=x[[3]]
    N<-dim(Q)[1]/2
    ret<-matrix(0,N,N)
    for(i in 1:(N)){
        for(j in (i):N){
            ret[j,i]<-ret[i,j]<- 1-mean(
            ((Q[2*i-1,]+Q[2*j-1,])==2)*(K[2*i-1,]!=K[2*j-1,])+
            ((Q[2*i-1,]+Q[2*j,])==2)  *(K[2*i-1,]!=K[2*j,]) +
            ((Q[2*i,]+Q[2*j-1,])==2)  *(K[2*i,]!=K[2*j-1,])+
            ((Q[2*i,]+Q[2*j,])==2)    *(K[2*i,]!=K[2*j,]))/4
        }
    }
    ret
}
pairwiseKinshipQScore<-function(x,delta=0.5){
    ## Expected value of the Pairwise Heterozygosity
    ## where x is the data list and the result is an N by N matrix
    ##
    Q=x[[1]]
    K=x[[3]]
    N<-dim(Q)[1]/2
    ret<-matrix(0,N,N)
    for(i in 1:(N)){
        for(j in (i):N){
            ret[j,i]<-ret[i,j]<- mean(
            (1-delta)*(
                (2-Q[2*i-1,]+Q[2*j-1,])+
                (2-Q[2*i-1,]+Q[2*j,])+
                (2-Q[2*i,]+Q[2*j-1,])+
                (2-Q[2*i,]+Q[2*j,]) ) +
            2*delta*(
                1-((Q[2*i-1,]+Q[2*j-1,])==2)*(K[2*i-1,]!=K[2*j-1,])+
                1-((Q[2*i-1,]+Q[2*j,])==2)  *(K[2*i-1,]!=K[2*j,]) +
                1-((Q[2*i,]+Q[2*j-1,])==2)  *(K[2*i,]!=K[2*j-1,])+
                1-((Q[2*i,]+Q[2*j,])==2)    *(K[2*i,]!=K[2*j,])
            ))/4
            }
    }
    ret
}

######################################################################
### Selection procedures
selIdentity<-function(x,frac=NA){
    ## Return a weighting of 1 for all individuals
    Q=rep(1,dim(x[[1]])[1]/2)
}
selRandom<-function(x,frac=0.5){
    ## Breed only a random fraction of the individuals
    N=dim(x[[1]])[1]/2
    ret=rep(0,N)
    ret[sample(1:N)[1:ceiling(N*frac)]]<-1
    ret
}
selMaxQ<-function(x,frac=0.5){
    ## Select the individuals with the lowest introgression
    Q=getQ(x)
    N=length(Q)
    ret=rep(0,N)
    ret[order(Q,decreasing=TRUE)[1:ceiling(frac*N)]]=1
    ret
}
selHetQwt<-function(x,frac=0.5){
    ## Select the individuals with the lowest introgression and the highest heterozyogisty within wildcat segments
    ## using an inverse weighting scheme
    wt=getExpectedSubpopulationHetScore(x)
    wt2=t(t(wt)/colMeans(wt))
    score=rowMeans(wt2)
    ret=rep(0,N)
    ret[order(score,decreasing=TRUE)[1:ceiling(frac*N)]]=1
    ret
}
selHetQ<-function(x,frac=0.5){ ## Doesn't make much sense
    ## Select the individuals with the lowest introgression and the highest heterozyogisty within wildcat segments
    Q=getSnpQ(x)
    h<-getHet(x)
    wt<-(Q==1)
    rowMeans(h * wt)
}
selHet<-function(x,frac=0.5){
    ## Select the individuals with the highest heterozygosity
    h<-getHet(x)
    score=rowMeans(h)
    ret=rep(0,N)
    ret[order(score,decreasing=TRUE)[1:ceiling(frac*N)]]=1
    ret
}
selKin<-function(x,frac=0.5){
    ## Select the individuals with the lowest kinship
    K<-getKinship(x)
    score=rowMeans(K)
    ret=rep(0,N)
    ret[order(score,decreasing=FALSE)[1:ceiling(frac*N)]]=1
    ret
}
selKinQ<-function(x,frac=0.5){
    ## Select the individuals with the lowest wildcat kinship
    K<-getSubpopulationKinship(x) 
    score=rowMeans(K)
    ret=rep(0,N)
    ret[order(score,decreasing=FALSE)[1:ceiling(frac*N)]]=1
    ret
}
selExpectedKin<-function(x,frac=0.5){
    ## Select the individuals with the lowest kinship
    K<-getExpectedKinshipCoeficient(x)
    score=rowMeans(K)
    ret=rep(0,N)
    ret[order(score,decreasing=FALSE)[1:ceiling(frac*N)]]=1
    ret
}
selExpectedKinQ<-function(x,frac=0.5){
    ## Select the individuals with the lowest wildcat kinship
    K<-getExpectedSubpopulationKinshipCoefficient(x)
    score=rowMeans(K)
    ret=rep(0,N)
    ret[order(score,decreasing=FALSE)[1:ceiling(frac*N)]]=1
    ret
}
selKinQwt<-function(x,frac=0.5){ ## Doesn't make much sense
    ## Select the individuals with the lowest wildcat kinship using an inverse weighting scheme
    wt<-getExpectedSubpopulationKinshipCoefficient(x)
    wt2=t(t(wt)/colMeans(wt))
    score=rowMeans(wt2)
    ret=rep(0,N)
    ret[order(score,decreasing=FALSE)[1:ceiling(frac*N)]]=1
    ret
}
