#include <Eigen/Dense>
#include "vbglmm.h"
#include <vector> 
#include <math.h> 
#include <assert.h> 
#include <cmath> 
#include <limits>
#include <time.h> 

#include <R_ext/Applic.h>
#include <Rcpp.h>

#define REV_GLOBAL 0
#define REV_REGRESSION 1
#define REV_LOCAL 2
#define REV_LOCAL_REGRESSION 3

#define OVERRELAX 1

#define LB_NAN_ERROR 2

using namespace std ;
using namespace Rcpp ;
using namespace Eigen ;

double log2pi = log(2.0*M_PI); 
 
double sqr(double x) { return x*x; } 

double logistic(double x){
  return 1.0/(1.0+exp(-x)); 
}

double logit(double x){
  return (x==0.0) ? 
    -20.0 
    : ( (x==1.0) ? 
	20.0 
	: log( x / (1.0-x)) ) ;
}

double log1plusexp(double x){
  return (x > 15.0) ? x : log(1.0+exp(x)); 
}

NumericVector twoByTwoSolve(const NumericMatrix A, NumericVector x){
  double det = A(0,0)*A(1,1)-A(0,1)*A(1,0);
  NumericVector y(2);
  y[0]=(A(1,1)*x[0]-A(0,1)*x[1])/det; 
  y[1]=(-A(1,0)*x[0]+A(0,0)*x[1])/det;
  return y; 
}
  
double fit_gamma(double s){
  double init_shape=(3.0-s+sqrt((s-3.0)*(s-3.0)+24.0*s))/(12.0*s);
  double shape=init_shape; 
  double old_shape=-1.0;
  int counter=0; 
  while (abs(old_shape-shape)>0.000001){
    old_shape=shape;
    shape=shape- (log(shape)-::Rf_digamma(shape)-s)/(1.0/shape - ::Rf_trigamma(shape)); 
    counter++;
    if (counter > 1000){
      Rcerr << "Warning: fit_gamma(" << s << ") not converging, returning " << init_shape << endl; 
      return init_shape;
    }
  }
  return shape; 
}
  
double fit_gamma(NumericVector &shapes, NumericVector &rates, double &rate){
  double meanP=0.0, meanLogP=0.0; 
  int N=shapes.size();
  for (int i=0;i<N;i++){
    meanP+=shapes[i]/rates[i]; 
    meanLogP+=::Rf_digamma(shapes[i])-log(rates[i]); 
  }
  meanP /= (double)N;

  meanLogP /= (double)N; 
  double shape=fit_gamma( log(meanP) - meanLogP ); 
  rate = shape / meanP; 
  return shape; 
}

double lb_gamma(NumericVector &shapes, NumericVector &rates, double shape, double rate){
  double meanP=0.0, meanLogP=0.0; 
  int N=shapes.size();
  for (int i=0;i<N;i++){
    meanP+=shapes[i]/rates[i]; 
    meanLogP+=::Rf_digamma(shapes[i])-log(rates[i]); 
  }
  return (shape-1.0)*meanLogP-rate*meanP+(double)N*(shape*log(rate)-::Rf_lgammafn(shape)); 
}

NumericMatrix outer(NumericVector &x){
  int n=x.size(); 
  NumericMatrix res(n,n); 
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      res(i,j)=x[i]*x[j];
  return res; 
}

class Vbglmm {
  bool learn_betas; 
  int it;
  bool debug;
  int num_loci;
  int max_its; 
 
  double converge_tol; 
  bool learn_rep, learnRepRep; 
  int trace_every; 
  int rev_model; 
  bool return_aux_variables;
  bool store_all_coeffs;
  double coeff_regulariser;
  NumericVector normalised_depth; 
  NumericVector rep_shapes;
  NumericVector rep_rates; 

  // stuff for localreg
  NumericVector logrep_mp, logrep_p; 

  vector<NumericVector> a; // part of the variational approximation for the logistic factor

  // natural parameter of q(g), where g is the logistic regression auxiliary variable
  vector<NumericVector> mean_prec; 
  vector<NumericVector> prec; 

  double rep_rep; 

  double global_rep_shape, global_rep_rate;   
  double random_effect_precision; 
  double rep_intercept, rep_slope; 

  vector<NumericVector> beta; // coefficients
  NumericVector standard_errors; 
  NumericVector lower_bounds; // lower bound on the log likehood per locus
  
  vector<LDLT<MatrixXd> > choleskys;
  
   // C++ structure
  vector<NumericVector> alt; // counts of alternative SNP
  vector<NumericVector> n; // total reads mappable to either allele at each locus for each individual
  //vector<vector<NumericVector> > x; // covariate
  //  typedef ListOf< ListOf< NumericVector> > xType ; 
  typedef vector< List > xType ; 
  xType x; 


  
  double local_bound(int locus_index, int sample_index, double pred, double Erep, double Elogrep){
    if (prec[locus_index][sample_index] <= 0.0) return -std::numeric_limits<double>::infinity();
    double v=1.0/prec[locus_index][sample_index]; 
    double m=mean_prec[locus_index][sample_index]*v; 
    double Elog1plus_exp_g=.5*a[locus_index][sample_index]*a[locus_index][sample_index]*v+log1plusexp(m+(1.0-2.0*a[locus_index][sample_index])*v*.5); 
    // < log likelihood >
    double lb=(double)alt[locus_index][sample_index]*m-(double)n[locus_index][sample_index]*Elog1plus_exp_g; 
    double err = pred-m; 
    double Eerr2 = err*err+v;
    
    // < log normal >  [note: cancelling log2pi terms]
    lb+=-.5*Erep*Eerr2+.5*Elogrep; 
    // - <log q>
    lb+=.5*(1.0+log(v)); 
    if (isnan(lb)) ::Rf_error("Lower bound is nan"); 
    return lb; 
  }

  double per_locus_bound(int locus_index)
  {
    double lb=0.0; 
    double Erep, Elogrep; 
    switch (rev_model){ 
    case REV_GLOBAL:
      Erep = random_effect_precision; 
      Elogrep = log(Erep);  
      break; 
    case REV_REGRESSION:
      Elogrep = rep_slope * normalised_depth[locus_index] + rep_intercept; 
      Erep = exp( Elogrep ); 
      break; 
    case REV_LOCAL:
      Erep = rep_shapes[locus_index]/rep_rates[locus_index]; 
      Elogrep = ::Rf_digamma(rep_shapes[locus_index])-log(rep_rates[locus_index]); 
      // < log G(Erep; a,b) >
      lb+=(global_rep_shape-1.0)*Elogrep - global_rep_rate * Erep + global_rep_shape * log(global_rep_rate) - ::Rf_lgammafn(global_rep_shape); 
      lb-= (rep_shapes[locus_index]-1.0)*::Rf_digamma(rep_shapes[locus_index])+log(rep_rates[locus_index])-rep_shapes[locus_index]-::Rf_lgammafn(rep_shapes[locus_index]); 
      break;
    case REV_LOCAL_REGRESSION:
      if (logrep_p[locus_index] <= 0.0) 
	return -std::numeric_limits<double>::infinity();
      double v=1.0/logrep_p[locus_index] ;
      double m=logrep_mp[locus_index]*v; 
      Erep = exp(m+.5*v); 
      Elogrep = m; 
      if (Erep==0.0) return -std::numeric_limits<double>::infinity(); 
      double rpred=rep_slope*normalised_depth[locus_index]+rep_intercept; 
      double err=rpred-m; 
      double Eerr2 = err*err+v; 
      // < log N(log local_rep; mx+c,v) > [cancelling log2pi]
      lb+=-.5*rep_rep*Eerr2+.5*log(rep_rep); 
      // - <log q>
      lb+=.5*(1.0+log(v)); 
      break;
     }
    if (Erep==0.0) ::Rf_error("Expected random effect precision is 0"); 
    //if (coeff_regulariser != 0.0)
    //  lb -= 
    for (int sample_index=0;sample_index<alt[locus_index].size();sample_index++){
      NumericVector x_ns=x[locus_index][sample_index]; 
      double pred=sum(beta[locus_index]*x_ns); 
      lb+=local_bound(locus_index, sample_index, pred, Erep, Elogrep); 
    }
    return lb; 
  }

  double lower_bound(){
    for (int locus_index=0;locus_index<num_loci;locus_index++)
      lower_bounds[locus_index]=per_locus_bound(locus_index); 
    return accumulate(lower_bounds.begin(),lower_bounds.end(),0.0); 
  }

  NumericVector lower_bound_grad(double rep_slope, double rep_intercept, NumericVector &normalised_depth, NumericVector &expected_err, NumericVector &total_terms, NumericMatrix &hess) {
    NumericVector grad(2); 
    for (int locus_index=0;locus_index<expected_err.size();locus_index++){
      double c=.5 * expected_err[locus_index] * exp( rep_slope * normalised_depth[locus_index] + rep_intercept );
      hess(0,0) += normalised_depth[locus_index] * normalised_depth[locus_index] * c ; 
      hess(0,1) += normalised_depth[locus_index] * c ; 
      hess(1,1) += c;
      grad(0) += .5 * total_terms[locus_index] * normalised_depth[locus_index] - c * normalised_depth[locus_index] ;
      grad(1) += .5 * total_terms[locus_index] - c;
    }
    hess(1,0)=hess(0,1); 
    return grad; 
  }    

  
  double update_single_g(int locus_index, int sample_index, double local_rep, double log_local_rep){
    NumericVector x_ns=x[locus_index][sample_index]; 
    double regression_mean=sum(beta[locus_index]*x_ns); 
    double v=1.0/prec[locus_index][sample_index]; 
    double m=mean_prec[locus_index][sample_index]/prec[locus_index][sample_index]; 
    double a_ns=a[locus_index][sample_index]; 
    double sig=logistic(m+(1.0-2.0*a_ns)*v*.5);
    double old_prec=prec[locus_index][sample_index]; 
    double old_mean_prec=mean_prec[locus_index][sample_index];
    double pf=n[locus_index][sample_index]*sig*(1.0-sig); 
    double p=pf+local_rep;
    double mp=m*pf+alt[locus_index][sample_index]-n[locus_index][sample_index]*sig+local_rep*regression_mean; 
    double new_local_bound; 
    double old_local_bound=local_bound(locus_index, sample_index, regression_mean, local_rep, log_local_rep); 
    // update q(g) [NB: this is still correct for random q(rep)]
    double step=OVERRELAX ? 1.5 : 1.0; 

    prec[locus_index][sample_index]=step*p+(1.0-step)*old_prec;
    mean_prec[locus_index][sample_index]=step*mp+(1.0-step)*old_mean_prec;

    new_local_bound=local_bound(locus_index, sample_index, regression_mean, local_rep, log_local_rep); 
    while (new_local_bound<old_local_bound){
      step *= (OVERRELAX ? .666666666 : .5); 
      prec[locus_index][sample_index]=step*p+(1.0-step)*old_prec;
      mean_prec[locus_index][sample_index]=step*mp+(1.0-step)*old_mean_prec;
      new_local_bound=local_bound(locus_index, sample_index, regression_mean, local_rep, log_local_rep); 
      if (step<0.000001){
	//if (debug) Rcout << "Step is very small, reg mean: " << regression_mean << endl;
	prec[locus_index][sample_index]=old_prec; 
	mean_prec[locus_index][sample_index]=old_mean_prec; 
	new_local_bound=local_bound(locus_index, sample_index, regression_mean, local_rep, log_local_rep); 
	break; 
      }
    }
    if (debug) if ((new_local_bound + 1e-8)<old_local_bound) Rcout << "GDI bound got worse, old: " << old_local_bound << " new: " << new_local_bound << endl; 
    m=mean_prec[locus_index][sample_index]/prec[locus_index][sample_index]; 
    v=1.0/prec[locus_index][sample_index]; 
    // update a
    //double check=.5*a_ns*a_ns*v+log1plusexp(m+(1.0-2.0*a_ns)*v*.5); 
    a[locus_index][sample_index]=logistic(m+(1.0-2.0*a_ns)*v*.5);
    if ( isnan( a[locus_index][sample_index] ) ) ::Rf_error("Auxiliary variable a is nan"); 
    double old_lb = new_local_bound; 
    new_local_bound=local_bound(locus_index, sample_index, regression_mean, local_rep, log_local_rep);  

    if (debug) if ((new_local_bound + 1e-8) < old_lb) Rcout << "Warning: updating a decreased lower bound, old: " << old_lb << " new: " << new_local_bound << endl;
    return new_local_bound; 
    //double check2=.5*a_ns*a_ns*v+log1plusexp(m+(1.0-2.0*a_ns)*v*.5); 
    //if (check2>check) Rcout << "update of a was bad: before " << check << " after " << check2 << endl;
  
  }

  void get_locus_error_terms(int locus_index, NumericVector &expected_err, NumericVector &total_terms){
    expected_err[locus_index]=0.0; 
    //if (it < 100) // DEBUGDEBUG
    int num_samples=alt[locus_index].size(); 
    for (int sample_index=0;sample_index<num_samples;sample_index++){
      double m=mean_prec[locus_index][sample_index]/prec[locus_index][sample_index]; 
      double v=1.0/prec[locus_index][sample_index]; 
      NumericVector x_ns=x[locus_index][sample_index]; 
      double pred=sum(beta[locus_index]*x_ns);  
      double err=pred-m;
      if (!isfinite(err)) ::Rf_error("Residual error is infinite"); 
      expected_err[locus_index]+=err*err+v;  
    }
    if (!isfinite(expected_err[locus_index])) ::Rf_error("Total esidual error is infinite"); 
    total_terms[locus_index]=num_samples; 
  }

  void update_locus(int locus_index, NumericVector &expected_err, NumericVector &total_terms){
    double local_rep, log_local_rep;
    double old_lb ; 
    switch (rev_model){
    case REV_GLOBAL:
      local_rep = random_effect_precision; 
      log_local_rep = log(local_rep); 
      break; 
    case REV_REGRESSION:
      log_local_rep = rep_slope * normalised_depth[locus_index] + rep_intercept; 
      local_rep = exp( log_local_rep ); 
      break; 
    case REV_LOCAL:
      local_rep = rep_shapes[locus_index]/rep_rates[locus_index]; 
      log_local_rep = ::Rf_digamma(rep_shapes[locus_index])-log(rep_rates[locus_index]); 
      break; 
    case REV_LOCAL_REGRESSION:
      double v=1.0/logrep_p[locus_index] ;
      double m=logrep_mp[locus_index]*v; 
      local_rep = exp(m+.5*v); 
      log_local_rep = m; 
      break;
    }
    int num_cov=beta[locus_index].size(); 
    NumericVector xg(num_cov); 
    
    int num_samples=alt[locus_index].size(); 

    if (debug) old_lb=per_locus_bound(locus_index); 
    for (int sample_index=0;sample_index<num_samples;sample_index++){
      update_single_g(locus_index, sample_index, local_rep, log_local_rep); 
      NumericVector x_ns=x[locus_index][sample_index]; 
      double m=mean_prec[locus_index][sample_index]/prec[locus_index][sample_index]; 
      xg += x_ns * m; 
    }
    if (debug){
      double new_lb=per_locus_bound(locus_index); 
      if ((new_lb+1.0e-5)<old_lb) 
	Rcout << "Warning: lb got worse after optimizing q(g) and a, old:" << old_lb << " new: " << new_lb << endl;   	
      old_lb=new_lb; 
    }
    if (learn_betas & (num_cov>0)) { // && it < 100){ // DEBUGDEBUG
      NumericVector old_beta=beta[locus_index];
      // update beta and means
      // TODO: might want to check if determinant is too close to zero
      VectorXd xgv=VectorXd::Zero( num_cov ); 
      MatrixXd xx=MatrixXd::Zero(num_cov,num_cov); 
      
      VectorXd res; 
      for (int ii=0;ii<num_cov;ii++)
	xgv[ii]=xg[ii]; 
      res=choleskys[locus_index].solve(xgv);
            
      for (int ii=0;ii<beta[locus_index].size();ii++)
	beta[locus_index][ii]=res[ii]; 
      if (!isfinite(beta[locus_index][0])) ::Rf_error("Regression coefficient is infinite"); 
      /*if (debug && (num_cov==2)){
	NumericVector check = twoByTwoSolve(xx,xg); 
	if (abs(check[0]-beta[locus_index][0])>0.0001) throw 1; 
	}*/
      //Rcout << "beta: " << beta[locus_index][0] << "," << beta[locus_index][1] << endl; 
    }

    if (debug){
      double new_lb=per_locus_bound(locus_index); 
      if ((new_lb+1.0e-5)<old_lb) 
	Rcout << "Warning: lb got worse after optimizing beta old:" << old_lb << " new: " << new_lb << endl;   	
      old_lb=new_lb; 
    }

    get_locus_error_terms(locus_index, expected_err, total_terms); 

    //    Rcout << expected_err[locus_index] << endl;
    //if (it < 100) // DEBUGDEBUG
    switch (rev_model){
    case REV_LOCAL:
      rep_shapes[locus_index]=global_rep_shape+.5*(double)num_samples; 
      rep_rates[locus_index]=global_rep_rate+.5*expected_err[locus_index]; 
      if (!isfinite(rep_rates[locus_index])) ::Rf_error("b parameter is infinite");
      break; 
    case REV_LOCAL_REGRESSION:
      double v=1.0/logrep_p[locus_index] ;
      double m=logrep_mp[locus_index]*v; 
      double old_p=logrep_p[locus_index]; 
      double old_mp=logrep_mp[locus_index]; 
      double pred_mean=rep_slope*normalised_depth[locus_index]+rep_intercept; 
      double pf=.5*expected_err[locus_index]*exp(m+.5*v); 
      double mpf=(m-1.0)*pf+.5*num_samples; 
      double p=pf+rep_rep; 
      double mp=mpf+pred_mean*rep_rep; 
      double step=OVERRELAX ? 1.5 : 1.0; 
      double old_lb=per_locus_bound(locus_index); 
      logrep_mp[locus_index]=step*mp+(1.0-step)*old_mp; 
      logrep_p[locus_index]=step*p+(1.0-step)*old_p; 
      
      while (per_locus_bound(locus_index) < old_lb){
	step *= (OVERRELAX ? .666666666 : .5); 
	logrep_mp[locus_index]=step*mp+(1.0-step)*old_mp; 
	logrep_p[locus_index]=step*p+(1.0-step)*old_p; 
	if (step < 0.00001){
	  if (debug) Rcout << "Warning: step size is very small updating q(log rep)" << endl; 
	  step=0.0;
	  logrep_mp[locus_index]=old_mp; 
	  logrep_p[locus_index]=old_p;
	  break;
	}
      }
    }
    if (debug){
      double new_lb=per_locus_bound(locus_index); 
      if ((new_lb+1.0e-5)<old_lb) 
	Rcout << "Warning: lb got worse after optimising rev model, old:" << old_lb << " new: " << new_lb << endl;   	
      old_lb=new_lb; 
    }


  }

  double one_iteration(){
    NumericVector expected_err(num_loci);
    NumericVector total_terms(num_loci); 
    double old_lb, lb; 

    for (int locus_index=0;locus_index<num_loci;locus_index++){
	update_locus(locus_index, expected_err, total_terms); 
    }

    // double al=min( (double)it/100.0, 1.0); 
    // global_rep_shape= al * 4.0 + (1.0-al)*1.0; 
    // global_rep_rate=al * 0.003 + (1.0-al)*1.0;

    // update random_effect_precision
    if (learn_rep){
      switch (rev_model){
      case REV_GLOBAL:
	if (debug) old_lb=lower_bound(); 
	random_effect_precision=accumulate(total_terms.begin(),total_terms.end(),0.0)/accumulate(expected_err.begin(),expected_err.end(),0.0); 
	lb = lower_bound();
	break; 
      case REV_REGRESSION:
	{
	  NumericMatrix hess(2,2);
	  NumericVector grad = lower_bound_grad(rep_slope,rep_intercept,normalised_depth,expected_err,total_terms,hess); 
	  NumericVector dir = twoByTwoSolve(hess,grad); 
	  double step=1.0 ;
	  NumericVector current = NumericVector::create( rep_slope, rep_intercept ); 
	  old_lb=lower_bound();
	  NumericVector proposed = current + step * dir; 
	  rep_slope=proposed[0]; 
	  rep_intercept=proposed[1]; 
	  lb = lower_bound(); 
	  while (lb  < old_lb ){
	    step *= .5;
	    proposed = current + step * dir;
	    rep_slope=proposed[0]; 
	    rep_intercept=proposed[1]; 
	    lb = lower_bound();
	    if (step < 0.00001){
	      Rcout << "Warning: step size very small in learning REP model"<< endl; 
	      break; 
	    }
	  }
	}
	break;
      case REV_LOCAL:
	double check_old_lb, check_lb; 
	if (debug) {
	  old_lb=lower_bound(); 
	  check_old_lb=lb_gamma(rep_shapes,rep_rates,global_rep_shape,global_rep_rate); 
	}
	global_rep_shape=fit_gamma(rep_shapes,rep_rates,global_rep_rate); 
	lb = lower_bound();
	if (debug){
	  check_lb=lb_gamma(rep_shapes,rep_rates,global_rep_shape,global_rep_rate); 
	  if (check_lb < check_old_lb){
	    Rcout << "Warning: check old lb: " << check_old_lb << " new: " << check_lb << endl; 
	  }
	  if (debug && ((lb+1.0e-3)<old_lb)){
	    Rcout << "Warning: lb got worse after optimizing random effect var, old:" << old_lb << " new:" << lb << endl; 
	    Rcout << "....... check old lb: " << check_old_lb << " new: " << check_lb << endl; 
	  }

	}
	break; 
      case REV_LOCAL_REGRESSION:
	if (debug) 
	  old_lb=lower_bound(); 
	
	NumericMatrix xx(2,2);
	NumericVector xm(2); 
	double Eerr2=0; 
	for (int locus_index=0;locus_index<num_loci;locus_index++){
	  double v=1.0/logrep_p[locus_index];
	  double m=logrep_mp[locus_index]*v; 
	  xx(0,0)+=(normalised_depth[locus_index]*normalised_depth[locus_index]); 
	  xx(0,1)+=normalised_depth[locus_index];
	  xm[0]+=normalised_depth[locus_index]*m; 
	  xm[1]+=m; 
	}
	xx(1,0)=xx(0,1); 
	xx(1,1)=(double)num_loci; 
	NumericVector temp=twoByTwoSolve(xx,xm); 
	if (1) { // (it<100){ // DEBUGDEBUG
	  rep_slope=temp[0];
	  rep_intercept=temp[1];
	}
	if (isnan(rep_slope)) {
	  Rcout << xx(0,0) << " " << xx(1,0) << " " << xx(1,1) << " " << xm[0] << " " << xm[1] << endl ;
	  ::Rf_error("Slope is nan");
	}
	
	if (debug){
	  lb=lower_bound(); 
	  if ((lb+1.0e-3)<old_lb)
	    Rcout << "Warning: lb got worse after optimizing rep slope/intercept, old:" << old_lb << " new:" << lb << endl; 
	  old_lb=lb; 
	}
	for (int locus_index=0;locus_index<num_loci;locus_index++){
	  double v=1.0/logrep_p[locus_index];
	  double m=logrep_mp[locus_index]*v; 
	  double pred=rep_slope*normalised_depth[locus_index]+rep_intercept; 
	  double err=pred-m; 
	  Eerr2+=err*err+v; 
	}

	if (learnRepRep)
	  rep_rep=((double)num_loci+1.0)/(Eerr2+1.0);
	lb=lower_bound(); 
	
	if (debug && ((lb+1.0e-3)<old_lb)) Rcout << "Warning: lb got worse after optimizing reprep, old:" << old_lb << " new:" << lb << endl; 
      }
      
      if (!isfinite(random_effect_precision)) { Rcerr << "err " << expected_err << " n " << total_terms << endl; ::Rf_error("Random effect precision is nan"); }
      
      if (debug && ((lb+1.0e-3)<old_lb)) Rcout << "Warning: lb got worse after optimizing random effect var, old:" << old_lb << " new:" << lb << endl; 
    } else {
      lb=lower_bound(); 
    }

    return lb; 
  }

public:

  SEXP run(){
    double previous_it_lb=lower_bound(), lb;
    vector<double> mlPerIteration, repInterceptPerIteration, repSlopePerIt, repRepPerIt; 
    List allCoefs; 
    std::deque<NumericVector> thetas; 
    bool converged=false; 
    time_t start = time(NULL); 
    for (it=0;it<max_its;it++){
      lb = one_iteration(); 
      converged= abs(lb-previous_it_lb) < converge_tol;
      
      
      if ((trace_every<10000) && ((it % trace_every) == 0 || converged)){ 
	Rcout << "it: " << it << " lb: " << lb ; 
	switch (rev_model){
	case REV_GLOBAL:
	  Rcout << " re_var: " << 1.0/random_effect_precision;
	  break;
	case REV_LOCAL:
	  {
	    double sd = (global_rep_shape > 2.0) ? (global_rep_rate / ( (global_rep_shape - 1.0) * sqrt( global_rep_shape - 2.0 ) )) : 0.0; 
	    Rcout << " rep prior=G( " << global_rep_shape << "," << global_rep_rate << ") var~=" << global_rep_rate / global_rep_shape << "+/-" << sd; 
	  }
	  break;
	case REV_LOCAL_REGRESSION:
	  Rcout << " rep_rep " << rep_rep ; 
	case REV_REGRESSION:
	  Rcout << " rep_slope: " << rep_slope << " rep_intercept: " << rep_intercept;
	  break;
	}
	Rcout << " lb change:" << lb-previous_it_lb << " (tol=" << converge_tol << ")" ;
	Rcout << " (" << (double)difftime(time(NULL),start) << "seconds)" << endl; 
      }
      R_CheckUserInterrupt(); 

      if (converged){
	if (trace_every < 10000)
	  Rcout << "Converged!" << endl; 
	break; 
      }
      
      previous_it_lb=lb;
      mlPerIteration.push_back(lb); 
      repInterceptPerIteration.push_back(rep_intercept); 
      repRepPerIt.push_back(rep_rep) ;
      repSlopePerIt.push_back(rep_slope); 
      if (store_all_coeffs){
	allCoefs.push_back(clone(wrap(beta)));
      } 
      
      R_CheckUserInterrupt(); 
      
    }
    
    
    List result_list=List::create(_("coeffs") = beta,
				  _("mlPerIteration") = mlPerIteration, 
				  _("log.likelihoods") = lower_bounds);
    if (learn_rep){
      result_list["repInterceptPerIteration"]=repInterceptPerIteration;
      result_list["repSlopePerIt"]=repSlopePerIt; 
      result_list["repRepPerIt"]=repRepPerIt;
    }
    if (store_all_coeffs)
      result_list["coeffsPerIteration"]=allCoefs; 
    switch (rev_model){
    case REV_GLOBAL:
      result_list["random.effect.variance"] = 1.0/random_effect_precision;
      break;
    case REV_LOCAL:
      result_list["rep.shapes"] = rep_shapes;
      result_list["rep.rates"] = rep_rates; 
      result_list["rep.global.shape"] = global_rep_shape;
      result_list["rep.global.rate"] = global_rep_rate;
      break;
    case REV_LOCAL_REGRESSION:
      result_list["rep.rep"]=rep_rep;
      result_list["logrepP"]=logrep_p; 
      result_list["logrepMP"]=logrep_mp;
    case REV_REGRESSION:
      result_list["rep.slope"] = rep_slope;
      result_list["rep.intercept"] = rep_intercept;
      break;
    }
    
    if (return_aux_variables){
      result_list["aux.variables.mp"]=mean_prec;
      result_list["aux.variables.p"]=prec;
    }
    return result_list; 
  }

  void init(SEXP alt_sexp, SEXP n_sexp, SEXP x_sexp, SEXP settings_sexp){

    // convert R to Rcpp
    alt=as< vector< NumericVector > >(alt_sexp); 
    n=as< vector< NumericVector > >(n_sexp); 

    List alt_rlist(alt_sexp); 
    List n_rlist(n_sexp);
    //List x_rlist(x_sexp);
    List settings_list(settings_sexp);
    learn_betas=as<bool>(settings_list["learnBetas"]);
    debug=as<bool>(settings_list["debug"]);
    num_loci = alt.size();
    max_its = as<int>(settings_list["max.iterations"]); 
    converge_tol = as<double>(settings_list["convergence.tolerance"]); 
    learn_rep=as<bool>(settings_list["learn.rev"]); 
    learnRepRep=as<bool>(settings_list["learnRepRep"]);
    trace_every=as<int>(settings_list["traceEvery"]); 
    rev_model=as<int>(settings_list["rev.model"]);
    coeff_regulariser=as<double>(settings_list["coeff.regulariser"]); 
    return_aux_variables=as<bool>(settings_list["return.aux.variables"]); 
    store_all_coeffs=as<bool>(settings_list["storeAllCoeffs"]); 

    if (trace_every < 10000)
      Rcout << "VBGLMM: num_loci: " << num_loci << " max_its: " << max_its << " tolerance: " << converge_tol << endl; 

    NumericVector temp(num_loci); 
    rep_shapes = clone(temp);
    rep_rates = clone(temp); 
    standard_errors=clone(temp); 
    lower_bounds=clone(temp); 
    logrep_p=clone(temp); 
    logrep_mp=clone(temp); 
    
    switch (rev_model){
    case REV_GLOBAL:
      random_effect_precision=1.0/as<double>(settings_list["random.effect.variance"]); 
      break; 
    case REV_LOCAL:
      global_rep_rate=as<double>(settings_list["rep.global.rate"]); 
      global_rep_shape=as<double>(settings_list["rep.global.shape"]); 
      break;
    case REV_LOCAL_REGRESSION: // fall throgugh is delibrate
      rep_rep=as<double>(settings_list["rep.rep"]); 
      //final_rep_rep=rep_rep; 
      //rep_rep=1.0; // DEBUGDEBUG
    case REV_REGRESSION:
      normalised_depth=as<NumericVector>(settings_list["normalised.depth"]); 
      rep_intercept=as<double>(settings_list["rep.intercept"]) ; 
      rep_slope=as<double>(settings_list["rep.slope"]); 
      break;
    }
    NumericVector precs(num_loci); 
    x=as< xType >( x_sexp ); 

    for (int locus_index=0;locus_index<num_loci;locus_index++){ 

      if (trace_every<10000 && ((locus_index % 1000)==0))
	Rcout << "Loading... locus " << locus_index << "/" << num_loci << endl; 
      
      //      NumericMatrix xi(
      int num_samples=alt[locus_index].size(); 

      NumericVector x_=x[locus_index][0];
      int num_cov=x_.size();
      NumericVector b(num_cov,0.0);
      beta.push_back(b); 
      
      switch (rev_model){
      case REV_GLOBAL:
	precs[locus_index]=random_effect_precision;
	break; 
      case REV_LOCAL:
	rep_shapes[locus_index]=global_rep_shape; 
	rep_rates[locus_index]=global_rep_rate; 
	precs[locus_index]=global_rep_shape/global_rep_rate; 
	break; 
      case REV_REGRESSION:
	precs[locus_index]=exp(rep_slope*normalised_depth[locus_index]+rep_intercept); 
	break; 
      case REV_LOCAL_REGRESSION:
	double pred_mean=rep_slope*normalised_depth[locus_index]+rep_intercept; 
	logrep_mp[locus_index]=pred_mean*rep_rep; 
	logrep_p[locus_index]=rep_rep;
	precs[locus_index]=exp(pred_mean+1.0/(2.0*rep_rep)); 
	break;
	
      }

      // precompute XX and corresponding Cholesky decompositions
      NumericMatrix xx(num_cov,num_cov);
      NumericMatrix xMatrix(num_samples,num_cov); 
      for (int sample_index=0;sample_index<num_samples;sample_index++){
	NumericVector x_ns=x[locus_index][sample_index];
	xx += outer(x_ns); 
	
      }

      MatrixXd xxX(num_cov,num_cov);
      for (int ii=0;ii<num_cov;ii++)
	for (int jj=0;jj<num_cov;jj++){
	  xx(ii,jj) += (coeff_regulariser>0.0 && ii==jj) ? coeff_regulariser : 0.0; 
	  xxX(ii,jj)=xx(ii,jj); 
	}
      if (num_cov>0){
	LDLT<MatrixXd> llt(xxX);
	double det=MatrixXd(llt.matrixL()).diagonal().prod(); 
	det = det * det; 
	if (det < 1e-8) Rcerr << "Warning: X is ill conditioned for locus " << locus_index << endl; 
	choleskys.push_back(llt); 
      } else {
	LDLT<MatrixXd> lltnull; 
	choleskys.push_back(lltnull);
      }

      NumericVector ai(num_samples,0.5); 
      a.push_back(ai);
      NumericVector mp(num_samples,0.0) ;
      mean_prec.push_back(mp); 
      NumericVector p(num_samples,precs[locus_index]); 
      prec.push_back(p); 
	
    }

    
    if (trace_every < 10000) Rcout << "Loaded" << endl;
    //    throw 1; 
  }
};

  RcppExport SEXP runvb(SEXP alt_sexp, SEXP n_sexp, SEXP x_sexp, SEXP settings_sexp) { 
    Vbglmm vbglmm;
    vbglmm.init(alt_sexp, n_sexp, x_sexp, settings_sexp); 
    return vbglmm.run();
  }
