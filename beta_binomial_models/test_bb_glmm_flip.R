require(rstan)
require(doMC)

N=20L
K=3L
Ti=4L

ns=array(10, dim=c(N,Ti,K))
ys=array( rbinom(N*Ti*K, 10, .5), dim=c(N,Ti,K))
mode(ns)="integer"
mode(ys)="integer"
require(abind)
xNull=aperm( abind( foreach(i=seq_len(Ti)) %do% matrix(1,N,1), along=3 ), c(3,1,2) )
dat=list(N=N,P=1L,T=Ti,K=K,ys=ys,ns=ns,x=xNull,concShape=1.01,concRate=0.01)

sm=stan_model(file="~/Dropbox/eagle/eagle/beta_binomial_models/bb_glmm_flips_rep_prior_gene.stan", save_dso=F, auto_write=F)
o=optimizing(sm, dat=dat)
