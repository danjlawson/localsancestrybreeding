require(MOSAIC)
library(vcfR)

wildcat_data<-read.vcfR("mosaic_scottish_cats_E3.vcf.gz") # VCF

tmp=load("localanc_scot_2way_1-36_1-19_130_60_0.99_58.RData")
## [1] "localanc"    "final.flips" "g.loc"     
tmp2=load("scot_2way_1-36_1-19_130_60_0.99_58.RData") 
## [1] "target"  "logfile" "Mu"      "lambda"  "theta"   "alpha"   "PI"     
## [8] "rho"     "A"       "NUMA"    "nchrno"  "chrnos"  "dr"      "NL"     
##[15] "kLL"     "acoancs" "all_Fst" "GpcM"   
alphamosaic=sapply(alpha,function(x)x[2]) # alpha, according to mosaic (domestic,wildcat) pairs

samples=read.csv("mosaic_scottish_cats_sample_info.csv",header=FALSE,as.is=TRUE)
## Format is name, class, wildcat proportion
dataorder=colnames(wildcat_data@gt)[-1] # extract order of VCF files (and ence mosaic)
rownames(samples)=gsub("WCQ0","WCQ",samples[,1]) # Fix annoying mismatch between labels

samples=samples[dataorder,]


sapply(localanc,dim)


dim(localanc[[16]])

## Compare local ancestry to global ancestry
## Check with chrnos=1:19
chrnos=16 # set this to the chromosome you are working on, e.g., for chromosome 16
local_pos=grid_to_pos(localanc,"INPUT2/",g.loc,chrnos)
## local_pos[[chr]][1,,] is the domestic ancestry
## local_pos[[chr]][2,,] is the wildcat ancestry
allrs=sapply(local_pos,function(x)rowSums(x[2,,]))
hapsums=rowSums(allrs) # Total local ancestry
indsums=colSums(matrix(hapsums,nrow=2))

## Sanity checking
thresh=0.5
par(mfrow=c(1,2))
plot(alphamosaic,indsums,col=1+(samples[,3]>thresh),xlab="Mosaic global ancestry estimate", ylab="Mosaic Chr16 ancestry estimate") ## Perfectly correlated
plot(alphamosaic,samples[,3],col=1+(samples[,3]>thresh),xlab="Mosaic global ancestry estimate", ylab="Prior ancestry estimate") ## Uncorrelated

## Extract data
localancestry=local_pos[[1]][2,,] # proportion wildcat
nsnps<-dim(localancestry)[2] # Number of SNPs
nindv<-dim(localancestry)[1]/2 # Number of individuals (number of haplotypes divided by 2)

## Who do we want?
mysamples=which(samples[,3]>thresh)
myhaps=rep(mysamples,each=2)*2+c(-1,0)
## Get local ancestry data
mylocalancestry=1-(localancestry[myhaps,]>0.8) ## Thresholded; 0 for wildcat, 1 for introgressed
## Get genotype data
gtdata=apply(wildcat_data@gt[,-1],1,function(x)as.numeric(unlist(strsplit(x,"|",fixed=T))))
mygtdata=gtdata[myhaps,]

system("mkdir -p ../chr16data")
write.csv(mygtdata,file="../chr16data/gtdata.csv")
write.csv(mylocalancestry,file="../chr16data/localancestry.csv")
write.csv(samples[mysamples,],file="../chr16data/samples.csv")


