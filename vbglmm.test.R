# it: 99 re_var: 0.00119508 flips_prior_logodds 3.21112 lb: -1637.28
require(vbglmm)
set.seed(1)
n.loci=20
alt=list()
n=list()
xFull=list()
xNull=list()
fullToFlip=list()
nullToFlip=list()
for (i in 1:n.loci){
    fullToFlip[[i]]=c(1.0,1.0)
    nullToFlip[[i]]=c(1.0)
    n.samples=rpois(1,10)+2
    xFull[[i]]=cbind(rnorm(n.samples),numeric(n.samples)+1.0)
    xNull[[i]]=cbind(numeric(n.samples)+1.0)
    n[[i]]=rpois(n.samples,500)
    alt[[i]]=rbinom(n.samples,n[[i]],runif(n.samples)*.2*runif(1)+0.3)
    for (j in 1:n.samples)
        alt[[i]][j]=min(alt[[i]],n[[i]][j]-alt[[i]])
}
#numIts=c(1000)
#
#ni=1000
#repRep=50
res=list()
repRep=1
                                        #for (ni in numIts)
#repRepRange=c(0.1,0.2,0.5,1,2,5,10,100,1000,500,2000,20,50)
#for (repRepIndex in 1:length(repRepRange)){
#    repRep=repRepRange[repRepIndex]
res=run.all(alt,n,xFull,xNull,max.its=100000,tol=0.01,debug=F,flips="none",learn.rev=T,rev=1.0,traceEvery=100,rev.model="local.regression",null.first=T,storeAllCoeffs=F,repRep=repRep)
#}
#save(file="varyRepRepPriorAnneal.RData",res,repRepRange)

#mls=unlist( lapply( res, function(r) r$res.null$mlPerIteration[ 1000 ] ))

#pdf("~/Dropbox/epistatic_ase/optRepRep.pdf") 
#plot(repRepRange,mls,log="x")
#dev.off()

#pvalues=lapply(res,function(x) x$p.values)

#full
#it: 999 rep_rep 11223.2 rep_slope: 0.339587 rep_intercept: 5.23137 lb: -75421.3
#null
#it: 999 rep_rep 5561.16 rep_slope: 0.350873 rep_intercept: 5.18445 lb: -75446.9

