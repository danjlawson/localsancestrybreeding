
## Simulation for breeding and introgression where:
## Genetic and local ancestry data are stored as a list(ancestry,genetic)


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

breedingIntrogress<-function(data,dataIntrogress,G,infracseq=rep(0,G),d=rep(0.01,dim(data[[1]])[2])){
    ## Breed at random, with a particular introgression profile, from a fixed panel of introgressing individuals
    N=dim(data[[1]])[1]/2
    NI=dim(dataIntrogress[[1]])[1]/2
    lastdata=data
    for(g in 1:G){
        infrac=infracseq[g]
        parents=data.frame(
            p1=sample(1:N,N,TRUE),
            p2=sample(1:N,N,TRUE),
            p3=sample(1:NI,N,TRUE),
            introgressed=sample(0:1,N,TRUE,prob=c(1-infrac,infrac)))
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
        lastdata=data
    }
    data
}

simf<-function(f0,Fst){
    ##    https://www.nature.com/articles/s41431-021-00873-2
    rbeta(length(f0),f0*(1-Fst)/Fst,(1-f0)*(1-Fst)/Fst)
}

breedingSel<-function(data,G,scorefn=scoreMeanQ,selfn=selIdentity,d=rep(0.01,dim(data[[1]])[2]),verbose=TRUE,...){
    ## Breed for G gens, applying a specified selection procedure for choosing parents
    lastdata=data
    score=allScores(0,data)
    N=dim(data[[1]])[1]/2
    for(g in 1:G){
        if(verbose) print(tail(score,1))
        wt=selfn(lastdata,...)
        parents=data.frame(
            p1=sample(1:N,N,TRUE,prob=wt),
            p2=sample(1:N,N,TRUE,prob=wt))
        for(i in 1:N){
            tp1=crossover(getinddata(lastdata,parents[i,1]),d=d)
            tp2=crossover(getinddata(lastdata,parents[i,2]),d=d)
            for(v in 1:length(data)){
                data[[v]][2*i-1,]=tp1[[v]]
                data[[v]][2*i,]=tp2[[v]]
            }
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
getNnative<-function(x){
    ## Return an N by L matrix reporting if each site is a
    a=1-x[[1]] ## admixture
    g=x[[2]]*2-1 ## genotype (-1,1)
    g=g*a ## genotype (-1 for wildcat 0, 0 for domestic, 1 for wildcat 1)
    N<-dim(g)[1]/2
    L<-dim(g)[2]
    ret<-matrix(0,N,L)
    for(i in 1:N){
        ## score 0 for no wildcat, 0.5 for 1 wildcat hap, 1 for 2 wildcat haps with same variant, 2 for 2 wildcat variants
        ret[i,]<-g[2*i,]^2+g[2*i-1,]^2 - abs(g[2*i-1,]+g[2*i,])/2
        ## head(cbind(g[2*i-1,],g[2*i,],ret[i,]))
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
scoreHet<-function(x){
    ## score is the average Heterozygosity
    h<-getHet(x)
    mean(h)
}
scoreHetQ<-function(x){
    ## score is the average Heterozygosity at wildcat sites
    Q=getSnpQ(x)
    h<-getHet(x)
    wt<-(Q==0)
    mean(h * wt)
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
    data.frame(g=g,Q=scoreMeanQ(x),
               Het=scoreHet(x),
               HetQ=scoreHetQ(x),
               Ae=scoreAe(x),
               Info=scoreInfo(x),
               InfoQ=scoreInfoQ(x),
               EffectiveLociQ=scoreEffectiveLociQ(x),
               EffectiveLoci=scoreEffectiveLoci(x))
}

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
selMinQ<-function(x,frac=0.5){
    ## Select the individuals with the lowest introgression
    Q=getQ(x)
    N=length(Q)
    ret=rep(0,N)
    ret[order(Q,decreasing=FALSE)[1:ceiling(frac*N)]]=1
    ret
}
selHetQ<-function(x,frac=0.5){
    ## Select the individuals with the lowest introgression and the highest heterozyogisty within wildcat segments
    ## Q=getSnpQ(x)
    ## h<-getHet(x)
    ## wt<-(Q==0)
    ## score=rowMeans(h * wt)
    wt=getNnative(x)
    wt2=t(t(wt)/colMeans(wt))
    score=rowMeans(wt2)
    ret=rep(0,N)
    ret[order(score,decreasing=TRUE)[1:ceiling(frac*N)]]=1
    ret
}
selHet<-function(x,frac=0.5){
    ## Select the individuals with the highest heterozygosity
    h<-getHet(x)
    score=rowMeans(h)
    ret=rep(0,N)
    ret[order(score,decreasing=TRUE)[1:ceiling(frac*N)]]=1
    ret
}
