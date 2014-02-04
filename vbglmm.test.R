require(vbglmm)
n.loci=100
alt=list()
n=list()
xFull=list()
xNull=list()
for (i in 1:n.loci){
    n.samples=rpois(1,10)+3
    xFull[[i]]=cbind(rnorm(n.samples),numeric(n.samples)+1.0)
    xNull[[i]]=cbind(numeric(n.samples)+1.0)
    n[[i]]=rpois(n.samples,100)
    alt[[i]]=rbinom(n.samples,n[[i]],numeric(n.samples)+0.4)
}

a=run.all(alt,n,xFull,xNull,max.its=100,tol=1.0,debug=T,flips="structured",rev.model="global",learn.rev=T,coeff.reg=1.0)

