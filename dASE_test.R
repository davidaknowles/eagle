# This is a test script for the Environment ASE GLMM
# In this case the model is used to test for differential ASE rather than association
# with an environment variable. 
# Data sampled from the generative model and then we run inference to attempt to
# recover the true "hits"

require(vbglmm) # load the package
set.seed(1) # for reproducibility
#-------------- generate some synthetic data -------------------------
n.loci=200 # number of loci
n.samples=100 # number of individuals
alt=list() # list of alternative counts
n=list() # list of total counts
xFull=list() # list of design matrices for the alternative models
xNull=list() # list of design matrices for the null models
# make "true" regression coefficients, most of which are 0

trueBeta=ifelse(runif(n.loci) < 0.05, rgamma(n.loci,shape=2,scale=.2)*(1-2*(runif(n.loci)<0.5)), 0.0)

logistic=function(x) 1/(1+exp(-x))
for (i in 1:n.loci){
    maf=runif(1)*.4+.1 # MAF between .1 and .5 for this locus
    # sample which individuals are heterozygous at this locus assuming HWE
    hap1=runif(n.samples)<maf 
    hap2=runif(n.samples)<maf
    hets=xor(hap1,hap2)
    numHets=sum(hets)
    ones=numeric(numHets)+1.0 
    xFull[[i]]=cbind(ones) # design matrix at this locus for alternative ("full") model
    xNull[[i]]=matrix(NA,ncol=0,nrow=numHets) # design matrix at this locus for null model
    n[[i]]=rpois(numHets,100*rgamma(numHets,shape=2,rate=2)) # sample read depth from overdispersed Poisson (NB?)
    p=logistic(trueBeta[i]+.3*rnorm(numHets)) # sample underlying probabilities, with overdispersion
    alt[[i]]=rbinom(numHets,n[[i]],p) # sample alternative counts
    #for (j in 1:numHets) # "min" model: doesn't work well with the way I sample betas here
    #    alt[[i]][j]=min(alt[[i]][j],n[[i]][j]-alt[[i]][j])
}

#---------------- run the model --------------------------
res=run.all(alt,n,xFull,xNull,max.its=100000,tol=0.0001,debug=F,flips="none",learn.rev=T,traceEvery=1,rev.model="local.regression",null.first=T,storeAllCoeffs=F)

cat("p-values for true hits:",res$p.values[trueBeta!=0],"\n")
cat("true hit betas: ",trueBeta[trueBeta!=0],"\n")

hist( res$p.values[ trueBeta==0.0 ], main="p values for null sites")
