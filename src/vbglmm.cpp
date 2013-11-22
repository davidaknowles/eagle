#include "vbglmm.h"
#include <vector> 
#include <math.h> 
#include <assert.h> 
#include <cmath> 

#define REV_GLOBAL 0
#define REV_REGRESSION 1
#define REV_LOCAL 2
#define REV_LOCAL_REGRESSION 3

using namespace std ;
using namespace Rcpp ;

double log2pi = log(2.0*M_PI); 
 
double logistic(double x){
  return 1.0/(1.0+exp(-x)); 
}

double log1plusexp(double x){
  return log(1.0+exp(x)); 
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
      cerr << "Warning: fit_gamma(" << s << ") not converging, returning " << init_shape << endl; 
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

class Vbglmm {

  bool debug;
  int num_loci;
  int max_its; 
  double converge_tol; 
  bool learn_rep; 
  bool learn_coeffs; 
  bool allow_flips; 
  bool trace; 
  int rev_model; 
  NumericVector normalised_depth; 
 
  NumericVector rep_shapes;
  NumericVector rep_rates; 

  // stuff for localreg
  NumericVector delta_to_logrep_mp, delta_to_logrep_p; 
  NumericVector logrep_mp, logrep_p; 

  double rep_rep; 

  double global_rep_shape, global_rep_rate;   
  double random_effect_precision; 
  double rep_intercept, rep_slope; 

  NumericVector beta; // coefficients
  NumericVector standard_errors; 
  NumericVector means; 
  NumericVector lower_bounds; // lower bound on the log likehood per locus
  
  vector<NumericVector> a; // part of the variational approximation for the logistic factor
  // natural parameter of q(g), where g is the logistic regression auxiliary variable
  vector<NumericVector> g_mean_prec, mean_prec_f; 
  vector<NumericVector> g_prec, prec_f; 

  vector<NumericVector> flips; 

  // C++ structure
  vector<NumericVector> alt; // counts of alternative SNP
  vector<NumericVector> n; // total reads mappable to either allele at each locus for each individual
  vector<NumericVector> x; // covariate

  double local_bound(double pred, double g_prec, double g_mean_prec, double a, double Erep, double Elogrep, int alt, int n){
    double lb=0.0; 
    double v=1.0/g_prec; 
    double m=g_mean_prec*v; 
    double Elog1plus_exp_g=.5*a*a*v+log1plusexp(m+(1.0-2.0*a)*v*.5); 
    // < log likelihood >
    lb+=(double)alt*m-(double)n*Elog1plus_exp_g; 
    double err = pred-m; 
    double Eerr2 = err*err+v; 
    // < log normal >  [note: cancelling log2pi terms]
    lb+=-.5*Erep*Eerr2+.5*Elogrep; 
    // - <log q>
    lb+=.5*(v+log(v)); 
    return lb; 
  }
  
  double per_locus_bound_local_rep(double beta, double means, NumericVector &a, NumericVector &g_mean_prec, NumericVector &g_prec, NumericVector &alt, NumericVector &n, NumericVector &x, double rep_shape, double rep_rate, double global_rep_shape, double global_rep_rate, NumericVector flips){
    double lb=0.0; 
    double Erep=rep_shape/rep_rate; 
    double Elogrep=::Rf_digamma(rep_shape)-log(rep_rate); 
    // < log G(local_rep; a,b) >
    lb+=(global_rep_shape-1.0)*Elogrep + global_rep_rate * Erep + global_rep_shape * log(global_rep_rate) - ::Rf_lgammafn(global_rep_shape); 
    // - < log q >
    lb-= (rep_shape-1.0)*::Rf_digamma(rep_shape)+log(rep_rate)-rep_shape-::Rf_lgammafn(rep_shape); 
    for (int sample_index=0;sample_index<n.size();sample_index++){
      double pred=(beta*x[sample_index]+means)*flips[sample_index]; 
      lb+=local_bound(pred, g_prec[sample_index], g_mean_prec[sample_index], a[sample_index], Erep, Elogrep, alt[sample_index], n[sample_index]); 
    }
    return lb; 
  }

  double per_locus_bound_local_reg(double beta, double means, NumericVector &a, NumericVector &g_mean_prec, NumericVector &g_prec, NumericVector &alt, NumericVector &n, NumericVector &x, double log_rep_mp, double log_rep_p, double normalised_depth_here, NumericVector flips){
    double lb=0.0; 
    double rpred=rep_slope*normalised_depth_here+rep_intercept; 
    double v=1.0/log_rep_p; 
    double m=log_rep_mp*v; 
    double err=rpred-m; 
    double Eerr2 = err*err+v; 
    // < log N(log local_rep; mx+c,v) > [calling log2pi]
    lb+=-.5*rep_rep*Eerr2+.5*log(rep_rep); 
    // - <log q>
    lb+=.5*(v+log(v)); 
    double rep_shape=fit_gamma(.5*v); 
    double rep_rate=rep_shape/exp(m+.5*v);           
    double Erep=exp(m+.5*v); 
    double Elogrep=m; 
    for (int sample_index=0;sample_index<n.size();sample_index++){
      double pred=(beta*x[sample_index]+means)*flips[sample_index]; 
      lb+=local_bound(pred, g_prec[sample_index], g_mean_prec[sample_index], a[sample_index], Erep, Elogrep, alt[sample_index], n[sample_index]); 
    }
    return lb; 
  }


  double per_locus_bound(double beta, double means, NumericVector &a, NumericVector &g_mean_prec, NumericVector &g_prec, NumericVector &alt, NumericVector &n, NumericVector &x, double random_effect_precision, NumericVector flips){
    double lb=0.0; 
    for (int sample_index=0;sample_index<n.size();sample_index++){
      double pred=(beta*x[sample_index]+means)*flips[sample_index]; 
      lb+=local_bound(pred, g_prec[sample_index], g_mean_prec[sample_index], a[sample_index], random_effect_precision, log(random_effect_precision), alt[sample_index], n[sample_index]); 
    }
    return lb; 
  }

  double per_locus_bound(int locus_index){
    switch( rev_model ){
    case REV_GLOBAL:
      return per_locus_bound(beta[locus_index], means[locus_index], a[locus_index], g_mean_prec[locus_index], g_prec[locus_index], alt[locus_index], n[locus_index], x[locus_index], random_effect_precision, flips[locus_index]); 
    case REV_REGRESSION:
      {
	double log_local_rep = rep_slope * normalised_depth[locus_index] + rep_intercept; 
	double local_rep = exp( log_local_rep ); 
	return per_locus_bound(beta[locus_index], means[locus_index], a[locus_index], g_mean_prec[locus_index], g_prec[locus_index], alt[locus_index], n[locus_index], x[locus_index], local_rep, flips[locus_index]); 
      }
      break; 
    case REV_LOCAL:
      return per_locus_bound_local_rep(beta[locus_index],means[locus_index], a[locus_index], g_mean_prec[locus_index], g_prec[locus_index], alt[locus_index], n[locus_index], x[locus_index], rep_shapes[locus_index], rep_rates[locus_index], global_rep_shape, global_rep_rate, flips[locus_index]); 
      break; 
    case REV_LOCAL_REGRESSION:
      return per_locus_bound_local_reg(beta[locus_index],means[locus_index], a[locus_index], g_mean_prec[locus_index], g_prec[locus_index], alt[locus_index], n[locus_index], x[locus_index], logrep_mp[locus_index], logrep_p[locus_index], normalised_depth[locus_index], flips[locus_index]); 
    }
  }

  double lower_bound_local_rep(NumericVector &beta, NumericVector &means, vector<NumericVector> &a, vector<NumericVector> &g_mean_prec, vector<NumericVector> &g_prec, vector<NumericVector> &alt, vector<NumericVector> &n, vector<NumericVector> &x, NumericVector &rep_shapes, NumericVector &rep_rates, double global_rep_shape, double global_rep_rate, vector<NumericVector > &flips) {
    int num_loci = n.size(); 
    for (int locus_index=0;locus_index<num_loci;locus_index++){
      lower_bounds[locus_index]=per_locus_bound_local_rep(beta[locus_index],means[locus_index], a[locus_index], g_mean_prec[locus_index], g_prec[locus_index], alt[locus_index], n[locus_index], x[locus_index], rep_shapes[locus_index], rep_rates[locus_index], global_rep_shape, global_rep_rate, flips[locus_index]); 
    }
    return accumulate(lower_bounds.begin(),lower_bounds.end(),0.0); 
  }

double lower_bound_local_reg(NumericVector &beta, NumericVector &means, vector<NumericVector> &a, vector<NumericVector> &g_mean_prec, vector<NumericVector> &g_prec, vector<NumericVector> &alt, vector<NumericVector> &n, vector<NumericVector> &x, vector<NumericVector > &flips) {
    int num_loci = n.size(); 
    for (int locus_index=0;locus_index<num_loci;locus_index++){
      lower_bounds[locus_index]=per_locus_bound_local_reg(beta[locus_index],means[locus_index], a[locus_index], g_mean_prec[locus_index], g_prec[locus_index], alt[locus_index], n[locus_index], x[locus_index], logrep_mp[locus_index], logrep_p[locus_index], normalised_depth[locus_index], flips[locus_index]); 
    }
    return accumulate(lower_bounds.begin(),lower_bounds.end(),0.0); 
}
    
  double lower_bound(NumericVector &beta, NumericVector &means, vector<NumericVector> &a, vector<NumericVector> &g_mean_prec, vector<NumericVector> &g_prec, vector<NumericVector> &alt, vector<NumericVector> &n, vector<NumericVector> &x, double random_effect_precision, vector<NumericVector > &flips) {
  int num_loci = n.size(); 
  for (int locus_index=0;locus_index<num_loci;locus_index++){
    lower_bounds[locus_index]=per_locus_bound(beta[locus_index],means[locus_index], a[locus_index], g_mean_prec[locus_index], g_prec[locus_index], alt[locus_index], n[locus_index], x[locus_index], random_effect_precision, flips[locus_index]); 
  }
  return accumulate(lower_bounds.begin(),lower_bounds.end(),0.0); 
  }    
  
  double lower_bound(NumericVector &beta, NumericVector &means, vector<NumericVector> &a, vector<NumericVector> &g_mean_prec, vector<NumericVector> &g_prec, vector<NumericVector> &alt, vector<NumericVector> &n, vector<NumericVector> &x, double rep_slope, double rep_intercept, NumericVector &normalised_depth, vector<NumericVector > &flips) {
    int num_loci = n.size(); 
    for (int locus_index=0;locus_index<num_loci;locus_index++){
      double local_rep= exp( rep_slope * normalised_depth[locus_index] + rep_intercept );
      lower_bounds[locus_index]=per_locus_bound(beta[locus_index],means[locus_index], a[locus_index], g_mean_prec[locus_index], g_prec[locus_index], alt[locus_index], n[locus_index], x[locus_index], local_rep, flips[locus_index]); 
    }
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
  
 
  double update_single_g(int locus_index, int sample_index, double local_rep, double log_local_rep, double &x_g, double &g_1, double &xx, double &x1){
    
    double x_ns=x[locus_index][sample_index];
    double regression_mean=(beta[locus_index]*x_ns+means[locus_index])*flips[locus_index][sample_index]; 
    double v=1.0/g_prec[locus_index][sample_index]; 
    double m=g_mean_prec[locus_index][sample_index]/g_prec[locus_index][sample_index]; 
    double a_ns=a[locus_index][sample_index]; 
    double sig=logistic(m+(1.0-2.0*a_ns)*v*.5);
    double old_prec_f=prec_f[locus_index][sample_index]; 
    double old_mean_prec_f=mean_prec_f[locus_index][sample_index]; 	
    double pf=n[locus_index][sample_index]*sig*(1.0-sig); 
    double mpf=m*pf+alt[locus_index][sample_index]-n[locus_index][sample_index]*sig; 
    double old_lb; 
    if (debug) old_lb=local_bound(regression_mean, g_prec[locus_index][sample_index], g_mean_prec[locus_index][sample_index], a[locus_index][sample_index], local_rep, log_local_rep, alt[locus_index][sample_index], n[locus_index][sample_index]);
    // update q(g) based on changes to rep and regression model (i.e. beta/mean)
    g_prec[locus_index][sample_index]=local_rep+old_prec_f;
    g_mean_prec[locus_index][sample_index]=regression_mean*local_rep+old_mean_prec_f; 
    double old_local_bound=local_bound(regression_mean, g_prec[locus_index][sample_index], g_mean_prec[locus_index][sample_index], a[locus_index][sample_index], local_rep, log_local_rep, alt[locus_index][sample_index], n[locus_index][sample_index]); 
    if (debug && (old_lb > (old_local_bound+0.001))) cout << "Warning: updating message from normal to g lowered lb, old: " << old_lb << " new: " << old_local_bound << endl; 
    // update q(g) [NB: this is still correct for random q(rep)]
    g_prec[locus_index][sample_index]=local_rep+pf;
    g_mean_prec[locus_index][sample_index]=regression_mean*local_rep+mpf; 
	
    double new_local_bound=local_bound(regression_mean, g_prec[locus_index][sample_index], g_mean_prec[locus_index][sample_index], a[locus_index][sample_index], local_rep, log_local_rep, alt[locus_index][sample_index], n[locus_index][sample_index]); 
    double step=1.0; 
    while (new_local_bound<old_local_bound){
      step *= .5; 
      g_prec[locus_index][sample_index]=local_rep+step*pf+(1.0-step)*old_prec_f;
      g_mean_prec[locus_index][sample_index]=regression_mean*local_rep+step*mpf+(1.0-step)*old_mean_prec_f;
      new_local_bound=local_bound(regression_mean, g_prec[locus_index][sample_index], g_mean_prec[locus_index][sample_index], a[locus_index][sample_index], local_rep, log_local_rep, alt[locus_index][sample_index], n[locus_index][sample_index]); 
      if (step<0.000001){
	if (debug) cout << "Step is very small" << endl;
	g_prec[locus_index][sample_index]=local_rep+old_prec_f;
	g_mean_prec[locus_index][sample_index]=regression_mean*local_rep+old_mean_prec_f;
	new_local_bound=local_bound(regression_mean, g_prec[locus_index][sample_index], g_mean_prec[locus_index][sample_index], a[locus_index][sample_index], local_rep, log_local_rep,  alt[locus_index][sample_index], n[locus_index][sample_index]); 
	break; 
      }
    } 
    prec_f[locus_index][sample_index]=g_prec[locus_index][sample_index]-local_rep; 
    mean_prec_f[locus_index][sample_index]=g_mean_prec[locus_index][sample_index]-local_rep*regression_mean; 
	  
    //if (new_local_bound<old_local_bound) cout << "GDI bound got worse, old: " << old_local_bound << " new: " << new_local_bound << endl; 
    m=g_mean_prec[locus_index][sample_index]/g_prec[locus_index][sample_index]; 
    v=1.0/g_prec[locus_index][sample_index]; 
    // update a
    //double check=.5*a_ns*a_ns*v+log1plusexp(m+(1.0-2.0*a_ns)*v*.5); 
    a[locus_index][sample_index]=logistic(m+(1.0-2.0*a_ns)*v*.5);
    a_ns=a[locus_index][sample_index]; 
    //double check2=.5*a_ns*a_ns*v+log1plusexp(m+(1.0-2.0*a_ns)*v*.5); 
    //if (check2>check) cout << "update of a was bad: before " << check << " after " << check2 << endl;
    if (allow_flips){
      // TODO could just look at abs(pred-m) i think
      double pure_reg_mean=beta[locus_index]*x_ns+means[locus_index];
      double flipOn = local_bound(-pure_reg_mean,g_prec[locus_index][sample_index], g_mean_prec[locus_index][sample_index], a[locus_index][sample_index], local_rep, log_local_rep, alt[locus_index][sample_index], n[locus_index][sample_index]); 
      double flipOff = local_bound(pure_reg_mean,g_prec[locus_index][sample_index], g_mean_prec[locus_index][sample_index], a[locus_index][sample_index], local_rep, log_local_rep, alt[locus_index][sample_index], n[locus_index][sample_index]); 
      flips[locus_index][sample_index]=(flipOn > flipOff) ? -1.0 : 1.0; 
    }
    x_g+=x_ns*m*flips[locus_index][sample_index];
    g_1+=m*flips[locus_index][sample_index]; 
    xx+=(x_ns*x_ns); // no flips here since -1*-1=1
    x1+=x_ns*flips[locus_index][sample_index];
  }

  double update_locus(int locus_index, NumericVector &expected_err, NumericVector &total_terms){
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
    case REV_LOCAL_REGRESSION:
      double v=1.0/logrep_p[locus_index] ;
      double m=logrep_mp[locus_index]*v; 
      local_rep = exp(m+.5*v); 
      log_local_rep = m; 
      break;
    }
    double x_g=0.0, g_1=0.0; 
    double xx=0.0,x1=0.0;
    int num_samples=alt[locus_index].size(); 
	
    if (debug) old_lb=per_locus_bound(locus_index); 
    for (int sample_index=0;sample_index<num_samples;sample_index++){
      update_single_g(locus_index, sample_index, local_rep, log_local_rep, x_g, g_1, xx, x1); 
    }
    if (debug){
      double new_lb=per_locus_bound(locus_index); 
      if ((new_lb+1.0e-3)<old_lb) 
	cout << "Warning: lb got worse after optimizing q(g) and a, old:" << old_lb << " new: " << new_lb << endl;   	
      old_lb=new_lb; 
    }
    double old_beta=beta[locus_index], old_mean=means[locus_index]; 
    // update beta and means
    double one1=(double)num_samples;
    double det=xx*one1-x1*x1; 
    if ((!learn_coeffs) || (det<1.0e-8)) { // x is constant
      beta[locus_index]=0.0; 
      means[locus_index]= (one1 > 0.0) ? g_1/one1 : 0.0; 
    }
    else {
      beta[locus_index]=(one1*x_g-x1*g_1)/det;
      means[locus_index]=-(x1*x_g-xx*g_1)/det; 
    }
    if (debug && (det > 1.0e-8) && (one1 > 0.0)){
      double new_lb=per_locus_bound(locus_index); 
      if ((new_lb+1.0e-4)<old_lb) {
	cout << "Warning: lb got worse after optimizing beta and mean, old:" << old_lb << " new: " << new_lb << endl; 
	cout << "x_g " << x_g << " g_1 " << g_1 << " xx " << xx << " x1 " << x1 << " ns " << num_samples << " det " << det << endl;
	cout << "x ";
	for (int sample_index=0;sample_index<num_samples;sample_index++) cout << x[locus_index][sample_index] << ",";
	cout << endl;
	cout << "locus: " << locus_index << " beta " << beta[locus_index] << " ("<<old_beta << ") means " << means[locus_index] << " (" << old_mean << ")" << endl; 
      }
      old_lb=new_lb; 
    }
      
    if ((!isfinite(beta[locus_index]) || (!isfinite(means[locus_index])))) { cerr << "x_g " << x_g << " g_1 " << g_1 << " xx " << xx << " x1 " << x1 << " ns " << num_samples << " det " << det << endl;  throw 1;  } 
    expected_err[locus_index]=0.0; 
    for (int sample_index=0;sample_index<num_samples;sample_index++){
      double m=g_mean_prec[locus_index][sample_index]/g_prec[locus_index][sample_index]; 
      double v=1.0/g_prec[locus_index][sample_index]; 
      double err=(beta[locus_index]*x[locus_index][sample_index]+means[locus_index])*flips[locus_index][sample_index]-m;
      if (!isfinite(err)) { cerr << "m " << m << " v " << v << " beta " << beta[locus_index] << " mean " << means[locus_index] << endl; throw 1; } 
      expected_err[locus_index]+=err*err+v; 
    }
    total_terms[locus_index]=num_samples; 
    //    cout << expected_err[locus_index] << endl;
    switch (rev_model){
    case REV_LOCAL:
      rep_shapes[locus_index]=global_rep_shape+.5*num_samples; 
      rep_rates[locus_index]=global_rep_rate+.5*expected_err[locus_index]; 
      break; 
    case REV_LOCAL_REGRESSION:
      double v=1.0/logrep_p[locus_index] ;
      double m=logrep_mp[locus_index]*v; 
      double b=.5*expected_err[locus_index]; 
      double aMinus1=.5*num_samples; 
      double pf=b*exp(m+.5*v); 
      double mpf=m*pf+aMinus1-pf; 
      double old_pf=delta_to_logrep_p[locus_index]; 
      double old_mpf=delta_to_logrep_mp[locus_index]; 
      double pred_mean=rep_slope*normalised_depth[locus_index]+rep_intercept; 
      double step=1.0; 
      double old_lb=per_locus_bound(locus_index); 
      logrep_mp[locus_index]=pred_mean*rep_rep+mpf; 
      logrep_p[locus_index]=rep_rep+pf; 
      while (per_locus_bound(locus_index) < old_lb){
	step*=0.5; 
	logrep_mp[locus_index]=pred_mean*rep_rep+step*mpf+(1.0-step)*old_mpf; 
	logrep_p[locus_index]=rep_rep+step*pf+(1.0-step)*old_pf; 
	if (step < 0.00001){
	  if (debug) cout << "Warning: step size is very small updating q(log rep)" << endl; 
	  step=0.0;
	  logrep_mp[locus_index]=pred_mean*rep_rep+old_mpf; 
	  logrep_p[locus_index]=rep_rep+old_pf;
	  break;
	}
      }
      delta_to_logrep_mp[locus_index]=step*mpf+(1.0-step)*old_mpf; 
      delta_to_logrep_p[locus_index]=step*pf+(1.0-step)*old_pf; 
    }

  }

  double one_iteration(){
    NumericVector expected_err(num_loci);
    NumericVector total_terms(num_loci); 
    double old_lb, lb; 
    for (int locus_index=0;locus_index<num_loci;locus_index++){
      update_locus(locus_index, expected_err, total_terms); 
    }

    // update random_effect_precision
    
    if (learn_rep){
      switch (rev_model){
      case REV_GLOBAL:
	if (debug) old_lb=lower_bound(beta,means,a,g_mean_prec,g_prec,alt,n,x,random_effect_precision,flips);
	random_effect_precision=accumulate(total_terms.begin(),total_terms.end(),0.0)/accumulate(expected_err.begin(),expected_err.end(),0.0); 
	lb = lower_bound(beta,means,a,g_mean_prec,g_prec,alt,n,x,random_effect_precision,flips);
	break; 
      case REV_REGRESSION:
	{
	  NumericMatrix hess(2,2);
	  NumericVector grad = lower_bound_grad(rep_slope,rep_intercept,normalised_depth,expected_err,total_terms,hess); 
	  NumericVector dir = twoByTwoSolve(hess,grad); 
	  double step=1.0 ;
	  NumericVector current = NumericVector::create( rep_slope, rep_intercept ); 
	  old_lb=lower_bound(beta,means,a,g_mean_prec,g_prec,alt,n,x,current[0],current[1],normalised_depth,flips);
	  NumericVector proposed = current + step * dir; 
	  lb = lower_bound(beta,means,a,g_mean_prec,g_prec,alt,n,x,proposed[0],proposed[1],normalised_depth,flips); 
	  while (lb  < old_lb ){
	    step *= .5;
	    proposed = current + step * dir;
	    lb = lower_bound(beta,means,a,g_mean_prec,g_prec,alt,n,x,proposed[0],proposed[1],normalised_depth,flips); 
	    if (step < 0.00001){
	      cout << "Warning: step size very small in learning REP model"<< endl; 
	      break; 
	    }
	  }
	  
	  rep_slope=proposed[0]; 
	  rep_intercept=proposed[1]; 
	}
	break;
      case REV_LOCAL:
	double check_old_lb, check_lb; 
	if (debug) {
	  old_lb=lower_bound_local_rep(beta,means,a,g_mean_prec,g_prec,alt,n,x,rep_shapes,rep_rates,global_rep_shape,global_rep_rate,flips); 
	  check_old_lb=lb_gamma(rep_shapes,rep_rates,global_rep_shape,global_rep_rate); 
	}
	global_rep_shape=fit_gamma(rep_shapes,rep_rates,global_rep_rate); 
	lb = lower_bound_local_rep(beta,means,a,g_mean_prec,g_prec,alt,n,x,rep_shapes,rep_rates,global_rep_shape,global_rep_rate,flips); 
	if (debug){
	  check_lb=lb_gamma(rep_shapes,rep_rates,global_rep_shape,global_rep_rate); 
	  if (check_lb < check_old_lb){
	    cout << "Warning: check old lb: " << check_old_lb << " new: " << check_lb << endl; 
	  }
	  if (debug && ((lb+1.0e-3)<old_lb)){
	    cout << "Warning: lb got worse after optimizing random effect var, old:" << old_lb << " new:" << lb << endl; 
	    cout << "....... check old lb: " << check_old_lb << " new: " << check_lb << endl; 
	  }

	}
	break; 
      case REV_LOCAL_REGRESSION:
	if (debug) 
	  old_lb=lower_bound_local_reg(beta,means,a,g_mean_prec,g_prec,alt,n,x,flips); 
	
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
	rep_slope=temp[0];
	rep_intercept=temp[1];
	if (isnan(rep_slope)) {
	  cout << xx(0,0) << " " << xx(1,0) << " " << xx(1,1) << " " << xm[0] << " " << xm[1] << endl ;
	  throw 1; 
	}
	for (int locus_index=0;locus_index<num_loci;locus_index++){
	  double v=1.0/logrep_p[locus_index];
	  double m=logrep_mp[locus_index]*v; 
	  double pred=rep_slope*normalised_depth[locus_index]+rep_intercept; 
	  double err=pred-m; 
	  Eerr2+=err*err+v; 
	}
	rep_rep=(double)num_loci/Eerr2; 
	for (int locus_index=0;locus_index<num_loci;locus_index++){
	  double pred_mean=rep_slope*normalised_depth[locus_index]+rep_intercept; 
	  logrep_mp[locus_index]=pred_mean*rep_rep+delta_to_logrep_mp[locus_index]; 
	  logrep_p[locus_index]=rep_rep+delta_to_logrep_p[locus_index]; 
	}
	lb=lower_bound_local_reg(beta,means,a,g_mean_prec,g_prec,alt,n,x,flips); 

      }
      
      if (!isfinite(random_effect_precision)) { cerr << "err " << expected_err << " n " << total_terms << endl; throw 1; }
      
      if (debug && ((lb+1.0e-3)<old_lb)) cout << "Warning: lb got worse after optimizing random effect var, old:" << old_lb << " new:" << lb << endl; 
    } else {
      switch (rev_model){
      case REV_GLOBAL:
	lb = lower_bound(beta,means,a,g_mean_prec,g_prec,alt,n,x,random_effect_precision,flips);
	break; 
      case REV_REGRESSION:
	lb = lower_bound(beta,means,a,g_mean_prec,g_prec,alt,n,x,rep_slope,rep_intercept,normalised_depth,flips); 
	break;
      case REV_LOCAL:
	lb = lower_bound_local_rep(beta,means,a,g_mean_prec,g_prec,alt,n,x,rep_shapes,rep_rates,global_rep_shape,global_rep_rate,flips); 
	break;
      case REV_LOCAL_REGRESSION:
	lb=lower_bound_local_reg(beta,means,a,g_mean_prec,g_prec,alt,n,x,flips); 
	break;
      }
    }
    return lb; 
  }

public:

  SEXP run(){
    double previous_it_lb=-1.0e30;   
    for (int it=0;it<max_its;it++){
      double lb = one_iteration(); 
      double beta_l2 = inner_product(beta.begin(),beta.end(),beta.begin(),0.0); 
      double means_l2 = inner_product(means.begin(),means.end(),means.begin(),0.0); 
      if (trace){ 
	cout << "it: " << it; 
	switch (rev_model){
	case REV_GLOBAL:
	  cout << " re_var: " << 1.0/random_effect_precision;
	  break;
	case REV_LOCAL:
	  {
	    double sd = (global_rep_shape > 2.0) ? (global_rep_rate / ( (global_rep_shape - 1.0) * sqrt( global_rep_shape - 2.0 ) )) : 0.0; 
	    cout << " rep prior=G( " << global_rep_shape << "," << global_rep_rate << ")~=" << global_rep_rate / global_rep_shape << "+/-" << sd; 
	  }
	  break;
	case REV_LOCAL_REGRESSION:
	  cout << " rep_rep " << rep_rep ; 
	case REV_REGRESSION:
	  cout << " rep_slope: " << rep_slope << " rep_intercept: " << rep_intercept;
	  break;
	}	
	cout << " beta sd " << sqrt(beta_l2/(double)num_loci) << " mean sd " << sqrt(means_l2/(double)num_loci) << " lb: " << lb << endl; 
      }
      R_CheckUserInterrupt(); 
      if (abs(lb-previous_it_lb) < converge_tol){
	cout << "Converged!" << endl; 
	break; 
      }
      previous_it_lb=lb;
    }

    NumericVector se2(num_loci); 
    // calculate standard errors and Wald p-values
    for (int locus_index=0;locus_index<num_loci;locus_index++){
      double xx=0.0,x1=0.0,one1=0.0;
      int num_samples=alt[locus_index].size(); 
    
      for (int sample_index=0;sample_index<num_samples;sample_index++){
	double x_ns=x[locus_index][sample_index]; 
	// next line uses that a=sigma(m+...)
	double c=((double)n[locus_index][sample_index])*a[locus_index][sample_index]*(1.0-a[locus_index][sample_index]); 
	xx+=c*x_ns*x_ns; 
	x1+=c*x_ns; 
	one1+=c; 
      }
      double det=xx*one1-x1*x1; 
      standard_errors[locus_index]=one1 / det;
    }
  
    // alternative calculation
    for (int locus_index=0;locus_index<num_loci;locus_index++){
      double xx=0.0,x1=0.0,one1=0.0;
      int num_samples=alt[locus_index].size(); 
    
      for (int sample_index=0;sample_index<num_samples;sample_index++){
	double x_ns=x[locus_index][sample_index]; 
	double m=g_mean_prec[locus_index][sample_index]/g_prec[locus_index][sample_index]; 
	double v=1.0/g_prec[locus_index][sample_index]; 
	double a_ns=a[locus_index][sample_index]; 
	double sig=logistic(m+(1.0-2.0*a_ns)*v*.5);
	double c=n[locus_index][sample_index]*sig*(1.0-sig);  
	double f = random_effect_precision/(random_effect_precision+c); 
	c*= f*f; 
	xx+=c*x_ns*x_ns; 
	x1+=c*x_ns; 
	one1+=c; 
      }
      double det=xx*one1-x1*x1; 
      se2[locus_index]=one1 / det;
    }

    NumericVector se_ep(num_loci), coeffs_ep(num_loci);

    for (int locus_index=0;locus_index<num_loci;locus_index++){
      double xx=0.0,x1=0.0,one1=0.0,mp1=0.0,mp2=0.0;
      int num_samples=alt[locus_index].size(); 
    
      for (int sample_index=0;sample_index<num_samples;sample_index++){
	double x_ns=x[locus_index][sample_index]; 
	double regression_mean=beta[locus_index]*x_ns+means[locus_index];
	double prec_f=g_prec[locus_index][sample_index]-random_effect_precision; 
	double mean_prec_f=g_mean_prec[locus_index][sample_index]-regression_mean*random_effect_precision; 
	double R = random_effect_precision/(random_effect_precision+prec_f); 
	double c=R*prec_f;
	mp1+=R*mean_prec_f*x_ns;
	mp2+=R*mean_prec_f; 
	xx+=c*x_ns*x_ns; 
	x1+=c*x_ns; 
	one1+=c; 
      }
      double det=xx*one1-x1*x1; 
      se_ep[locus_index]=one1 / det;
      coeffs_ep[locus_index]= (mp1*one1-mp2*x1)/det;
    }

    return List::create(_("coeffs") = beta,
			_("means") = means,
			_("standard.errors") = standard_errors, 
			_("alt.se") = se2, 
			_("se.ep") = se_ep,
			_("coeffs.ep") = coeffs_ep,
			_("flips") = flips,
			_("rep.slope") = rep_slope,
			_("rep.intercept") = rep_intercept,
			_("rep.rep")=rep_rep, 
			_("rep.shapes") = rep_shapes,
			_("rep.rates") = rep_rates, 
			_("random.effect.variance") = 1.0/random_effect_precision,
			_("rep.global.shape") = global_rep_shape,
			_("rep.global.rate") = global_rep_rate,
			_("log.likelihoods") = lower_bounds);
  }

  Vbglmm(SEXP alt_sexp, SEXP n_sexp, SEXP x_sexp, SEXP settings_sexp){

    // convert R to Rcpp
    List alt_rlist(alt_sexp); 
    List n_rlist(n_sexp);
    List x_rlist(x_sexp);
    
    List settings_list(settings_sexp);
    
    debug=as<bool>(settings_list["debug"]);   
    num_loci = alt_rlist.size(); 
    max_its = as<int>(settings_list["max.iterations"]); 
    converge_tol = as<double>(settings_list["convergence.tolerance"]); 
    learn_rep=as<bool>(settings_list["learn.rev"]); 
    learn_coeffs=as<bool>(settings_list["learn.coeffs"]);
    allow_flips=as<bool>(settings_list["allow.flips"]); 
    trace=as<bool>(settings_list["trace"]); 
    rev_model=as<int>(settings_list["rev.model"]); 
 
    cout << "VBGLMM: num_loci: " << num_loci << " max_its: " << max_its << " tolerance: " << converge_tol << endl; 

    NumericVector temp(num_loci); 
    rep_shapes = clone(temp);
    rep_rates = clone(temp); 
    beta=clone(temp); 
    standard_errors=clone(temp); 
    means=clone(temp); 
    lower_bounds=clone(temp); 
    delta_to_logrep_p=clone(temp); 
    delta_to_logrep_mp=clone(temp); 
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
    case REV_LOCAL_REGRESSION: // fall through is delibrate
      rep_rep=as<double>(settings_list["rep.rep"]); 
    case REV_REGRESSION:
      normalised_depth=as<NumericVector>(settings_list["normalised.depth"]); 
      rep_intercept=as<double>(settings_list["rep.intercept"]) ; 
      rep_slope=as<double>(settings_list["rep.slope"]); 
      break;
    }

    for (int locus_index=0;locus_index<num_loci;locus_index++){
      NumericVector alti((SEXP)alt_rlist[locus_index]); 
      alt.push_back(alti); 
      NumericVector ni((SEXP)n_rlist[locus_index]); 
      n.push_back(ni); 
      NumericVector xi((SEXP)x_rlist[locus_index]); 
      x.push_back(xi); 
      int num_samples=alti.size(); 
      NumericVector ai(num_samples,0.5); 
      a.push_back(ai);
      NumericVector mp(num_samples,0.0) ;
      g_mean_prec.push_back(mp); 
      NumericVector p(num_samples,random_effect_precision); 
      g_prec.push_back(p); 
      NumericVector mpf(num_samples,0.0) ;
      mean_prec_f.push_back(mpf); 
      NumericVector pf(num_samples,0.0); 
      prec_f.push_back(pf); 
      NumericVector f(num_samples,1.0); 
      flips.push_back(f); 
      switch (rev_model){
      case REV_LOCAL:
	rep_shapes[locus_index]=global_rep_shape; 
	rep_rates[locus_index]=global_rep_rate; 
	break; 
      case REV_LOCAL_REGRESSION:
	double pred_mean=rep_slope*normalised_depth[locus_index]+rep_intercept; 
	logrep_mp[locus_index]=pred_mean*rep_rep; 
	logrep_p[locus_index]=rep_rep;
	break;
      }
      for (int sample_index=0;sample_index<num_samples;sample_index++){ // TODO: probably don't need to do this? 
	g_mean_prec[locus_index][sample_index]=0.0; 
	g_prec[locus_index][sample_index]=1.0;
	a[locus_index][sample_index]=.5; 
	flips[locus_index][sample_index]=1.0; 
      }
    }
  }
};

  RcppExport SEXP runvb(SEXP alt_sexp, SEXP n_sexp, SEXP x_sexp, SEXP settings_sexp) { 
    Vbglmm vbglmm(alt_sexp, n_sexp, x_sexp, settings_sexp); 
    return vbglmm.run();
  }
