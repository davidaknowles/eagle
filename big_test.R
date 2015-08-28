# This is a test script for the EAGLE, the Environment-Ase Generalized Linear modEl
# Data sampled from the generative model and then we run inference to attempt to
# recover the true "hits"

require(eagle) # load the package
set.seed(1) # for reproducibility

#-------------- generate some synthetic data -------------------------
n.loci=300 # number of loci
n.samples=200 # number of individuals
alt=matrix(NA,n.samples,n.loci) # matrix of alternative counts
totalreads=matrix(NA,n.samples,n.loci) # matrix of total counts
het=matrix(NA,n.samples,n.loci) # matrix of which individuals are hets at each locus
environmentVar=runif(n.samples) # values of the environment variable, e.g. age
eqtlGenotypes=matrix(NA,n.samples,n.loci) # genotypes (0,1,2) of the lead eQTL for each gene
# make "true" regression coefficients, most of which are 0
#trueBeta=ifelse(runif(n.loci) < 0.05, rgamma(n.loci,shape=5,rate=5), 0.0) 
trueBeta=ifelse(runif(n.loci) < 0.05, runif(n.loci), 0.0) 
logistic=function(x) 1/(1+exp(-x))
for (i in 1:n.loci){
    maf=runif(1)*.4+.1 # MAF between .1 and .5 for this locus
    # sample which individuals are heterozygous at this locus assuming HWE
    hap1=runif(n.samples)<maf 
    hap2=runif(n.samples)<maf
    het[,i]=xor(hap1,hap2)
    totalreads[,i]=rpois(n.samples,100*rgamma(n.samples,shape=2,rate=2)) # sample read depth from overdispersed Poisson (NB?)
    logOdds=.3+trueBeta[i]*environmentVar+.3*rnorm(n.samples)
    if (runif(1)<.5){ # add eQTL
      eSNPmaf=runif(1)*.4+.1 
      ehap1=runif(n.samples)<eSNPmaf 
      ehap2=runif(n.samples)<eSNPmaf
      eqtlGenotypes[,i]=ehap1+ehap2
      logOdds=logOdds+.3*(eqtlGenotypes[,i]==1)
    }
    p=logistic(logOdds) # sample underlying probabilities, with overdispersion
    alt[,i]=rbinom(n.samples,totalreads[,i],p) # sample alternative counts
}
alt[het!=1]=NA
ref=totalreads-alt

#------------------ filters -------------------
alt.list <- list()
n.list <- list()
x.null <- list()
x.full <- list()
original.index <- list()
count.cutoff=3
prop.cutoff=0.01
prop.mono.cutoff=.5
non.problematic.counter=1
for (snp.index in 1:n.loci) {
  if (snp.index %% 1000 == 0) print(snp.index)
  valid=het[,snp.index]==1 & totalreads[,snp.index]>5
  # check we have at least 20 valid samples
  if (sum(valid)<20)
    next
  a=alt[valid,snp.index]
  r=ref[valid,snp.index]
  heteq= eqtlGenotypes[valid,snp.index]==1
  x=environmentVar[valid]
  # check not too many are in one group (e.g. female, non-smokers)
  if ((length(x)-max(table(x)))<10) 
    next
  n=a+r
  # check there isn't too much mono-allelic expression
  min.a.r=apply(cbind(a,r),1,min)
  is.mono=(min.a.r<count.cutoff)|((as.double(min.a.r)/n)<prop.cutoff)
  if (mean(is.mono)>prop.mono.cutoff)
    next
  alt.list[[non.problematic.counter]]=pmin(a,n-a) # min model
  n.list[[non.problematic.counter]]=n
  original.index[[non.problematic.counter]]=snp.index
  num.samples=length(x)
  x.full[[non.problematic.counter]]=cbind(env=x,intercept=numeric(num.samples)+1.0)
  x.null[[non.problematic.counter]]=cbind(intercept=numeric(num.samples)+1.0)
  # add eQTL heterozygosity
  if (!any(is.na(heteq))) if ((sd(heteq)>0) & (min(table(heteq))>5)){ 
    x.full[[non.problematic.counter]]=cbind(x.full[[non.problematic.counter]],eqtl=heteq)
    x.null[[non.problematic.counter]]=cbind(x.null[[non.problematic.counter]],eqtl=heteq)
    
    # lots of conditions to decide whether to add an interaction terms
    if (min(table(heteq))>=10) { # only if at least 10 hets and 10 non-hets
      e1=x[heteq==1.0]
      e0=x[heteq==0.0]
      if ((length(e1)-max(table(e1)))>=10 & (length(e0)-max(table(e0)))>=10){ 
        if (length(unique(x))>2 | min(table(x,heteq))>=5){ 
          temp=cbind(x.full[[non.problematic.counter]],interaction=x*heteq)
          if (det(t(temp) %*% temp) > 1e-8)
            x.full[[non.problematic.counter]]=temp
        }
      }
    }
    
  }
  non.problematic.counter=non.problematic.counter+1
}
original.index=unlist(original.index)

trueBetaAtTested=trueBeta[original.index]

#---------------- run the model --------------------------
s=eagle.settings()
s$debug=F
#s$rev.model=as.integer(3) # local regression
s$rev.model=2
s$normalised.depth=scale(log10(unlist(lapply(n,mean)))) # only required for rev.model=3
s$max.iterations=10000
s$convergence.tolerance=.001
s$coeff.regulariser=0.0

# learnt parameters from the DGN dataset
s$rep.global.shape=1.0
s$rep.global.rate=0.0033

s$traceEvery=1

system.time( res <- eagle.helper(alt.list,n.list,x.full,x.null,s) ) # 4s

cat("p-values for true hits:",res$p.values[trueBetaAtTested!=0],"\n")
cat("true hit betas: ",trueBetaAtTested[trueBetaAtTested!=0],"\n")

hist( res$p.values[ trueBetaAtTested==0.0 ], main="p values for null sites")

i=which.max(trueBetaAtTested)
plot(x.full[[i]][,1],alt.list[[i]]/n.list[[i]],ylim=c(0,1))
