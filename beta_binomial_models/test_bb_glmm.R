source("betabinomial_glmm.R")
require(doMC)
registerDoMC(7)

testbbglmm=function(seed=1,random_effect_sd=.3, fixed_effect=.5,nsamp=50,av_counts=30) {
  set.seed(seed)
  ns=rpois(nsamp * 2,lambda=av_counts)
  z=rep(1:nsamp,2)
  u=rnorm(nsamp)
  x=c(rep(0,nsamp),rep(1,nsamp))
  ys=rbinom(nsamp * 2,ns,logistic(.3+.3*rnorm(2*nsamp)+random_effect_sd*u[z]+fixed_effect*x))
  xNull=cbind(numeric(nsamp*2)+1)
  xFull=cbind(xNull,x)
  c(betaBinomialGLMM(ys,ns,xFull,xNull,z)$lrtp, betaBinomialGLM(ys,ns,xFull,xNull)$lrtp[1])
}

under_null=foreach(i=1:100) %dopar% { testbbglmm(seed=i, fixed_effect=0) }
pv=do.call(rbind,under_null)
qplot(pv[,2],pv[,1]) + geom_abline(intercept=0,slope=1) + xlab("GLM -log10(p)") + ylab("GLMM -log10(p)") + theme_bw(base_size=16)

under_full=foreach(i=1:100) %dopar% { testbbglmm(seed=i, fixed_effect=.5) }
pvals=as.data.frame(do.call(rbind,under_full))
colnames(pvals)=c("GLMM","GLM")
qplot(-log10(pvals$GLM),-log10(pvals$GLMM))+geom_abline(intercept=0,slope=1)+theme_bw(base_size = 16)+ xlab("GLM -log10(p)") + ylab("GLMM -log10(p)")

