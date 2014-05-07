require(vbglmm)
set.seed(1)
#-------------- generate some synthetic data -------------------------
n.loci=1000 # number of loci
n.samples=200 # number of individuals
alt=list() # list of alternative counts
n=list() # list of total counts
xFull=list() # list of design matrices for the alternative models
xNull=list() # list of design matrices for the null models
environmentVar=rnorm(n.samples) # values of the environment variable, e.g. age
trueBeta=ifelse(runif(n.loci) < 0.02, rgamma(n.loci,shape=2,rate=1), 0.0)
logistic=function(x) 1/(1+exp(-x))
for (i in 1:n.loci){
    maf=runif(1)*.4+.1 # MAF between .1 and .5
    hap1=runif(n.samples)<maf
    hap2=runif(n.samples)<maf
    hets=xor(hap1,hap2)
    numHets=sum(hets)
    print(numHets)
    ones=numeric(numHets)+1.0
    x=environmentVar[hets]
    xFull[[i]]=cbind(x,ones)
    xNull[[i]]=cbind(ones)
    n[[i]]=rpois(numHets,100*rgamma(numHets,shape=2,rate=2)) # overdispersed Poisson (NB?)
    p=logistic(trueBeta[i]*x+.3*rnorm(numHets)) 
    alt[[i]]=rbinom(numHets,n[[i]],p)
    #for (j in 1:numHets)
    #    alt[[i]][j]=min(alt[[i]][j],n[[i]][j]-alt[[i]][j])
}
res=run.all(alt,n,xFull,xNull,max.its=100000,tol=0.0001,debug=F,flips="none",learn.rev=T,traceEvery=1,rev.model="local.regression",null.first=T,storeAllCoeffs=F)
save(file="testNoRerun.RData",res,trueBeta)
