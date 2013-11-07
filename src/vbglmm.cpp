#include "vbglmm.h"
#include <vector> 
#include <math.h> 
#include <assert.h> 
#include <cmath> 

using namespace std ;
using namespace Rcpp ;
 
double logistic(double x){
  return 1.0/(1.0+exp(-x)); 
}

double log1plusexp(double x){
  return log(1.0+exp(x)); 
}

double local_bound(double pred, double g_prec, double g_mean_prec, double a, double random_effect_precision, int alt, int n){
  double lb=0.0; 
  double v=1.0/g_prec; 
  double m=g_mean_prec*v; 
  double Elog1plus_exp_g=.5*a*a*v+log1plusexp(m+(1.0-2.0*a)*v*.5); 
  // < log likelihood >
  lb+=(double)alt*m-(double)n*Elog1plus_exp_g; 
  // < log normal > 
  double err = pred-m; 
  double Eerr2 = err*err+v; 
  lb+=-.5*random_effect_precision*Eerr2-.5*log(2.0*M_PI/random_effect_precision); 
  // - <log q>
  lb+=.5*(v+log(2.0*M_PI*v)); 
  return lb; 
}

double per_locus_bound(double beta, double means, NumericVector &a, NumericVector &g_mean_prec, NumericVector &g_prec, NumericVector &alt, NumericVector &n, NumericVector &x, double random_effect_precision, NumericVector flips){
  double lb=0.0; 
  for (int sample_index=0;sample_index<n.size();sample_index++){
    double pred=(beta*x[sample_index]+means)*flips[sample_index]; 
    lb+=local_bound(pred, g_prec[sample_index], g_mean_prec[sample_index], a[sample_index], random_effect_precision, alt[sample_index], n[sample_index]); 
  }
  return lb; 
}

double lower_bound(NumericVector &beta, NumericVector &means, vector<NumericVector> &a, vector<NumericVector> &g_mean_prec, vector<NumericVector> &g_prec, vector<NumericVector> &alt, vector<NumericVector> &n, vector<NumericVector> &x, double random_effect_precision, vector<NumericVector > &flips) {
  int num_loci = n.size(); 
  double lb=0.0; 
  for (int locus_index=0;locus_index<num_loci;locus_index++){
    lb+=per_locus_bound(beta[locus_index],means[locus_index], a[locus_index], g_mean_prec[locus_index], g_prec[locus_index], alt[locus_index], n[locus_index], x[locus_index], random_effect_precision, flips[locus_index]); 
  }
  return lb; 
}    

double lower_bound(NumericVector &beta, NumericVector &means, vector<NumericVector> &a, vector<NumericVector> &g_mean_prec, vector<NumericVector> &g_prec, vector<NumericVector> &alt, vector<NumericVector> &n, vector<NumericVector> &x, double rep_slope, double rep_intercept, NumericVector &normalised_depth, vector<NumericVector > &flips) {
  int num_loci = n.size(); 
  double lb=0.0; 
  for (int locus_index=0;locus_index<num_loci;locus_index++){
    double local_rep= exp( rep_slope * normalised_depth[locus_index] + rep_intercept );
    lb+=per_locus_bound(beta[locus_index],means[locus_index], a[locus_index], g_mean_prec[locus_index], g_prec[locus_index], alt[locus_index], n[locus_index], x[locus_index], local_rep, flips[locus_index]); 
  }
  return lb; 
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


NumericVector twoByTwoSolve(NumericMatrix A, NumericVector x){
  double det = A(0,0)*A(1,1)-A(0,1)*A(1,0);
  NumericVector y(2);
  y[0]=(A(1,1)*x[0]-A(0,1)*x[1])/det; 
  y[1]=(-A(1,0)*x[0]+A(0,0)*x[1])/det;
  return y; 
}

RcppExport SEXP runvb(SEXP alt_sexp, SEXP n_sexp, SEXP x_sexp, SEXP settings_sexp) { 
  // convert R to Rcpp
  List alt_rlist(alt_sexp); 
  List n_rlist(n_sexp);
  List x_rlist(x_sexp);

  List settings_list(settings_sexp);

  bool debug=as<bool>(settings_list["debug"]);   
  int num_loci = alt_rlist.size(); 
  int max_its = as<int>(settings_list["max.iterations"]); 
  double converge_tol = as<double>(settings_list["convergence.tolerance"]); 
  double random_effect_precision=1.0/as<double>(settings_list["random.effect.variance"]); 
  bool learn_rep=as<bool>(settings_list["learn.rev"]); 
  bool learn_coeffs=as<bool>(settings_list["learn.coeffs"]);
  bool allow_flips=as<bool>(settings_list["allow.flips"]); 
  bool trace=as<bool>(settings_list["trace"]); 
  bool dependent_rev=as<bool>(settings_list["dependent.rev"]); 
  NumericVector normalised_depth; 
  
  cout << "VBGLMM: num_loci: " << num_loci << " max_its: " << max_its << " tolerance: " << converge_tol << endl; 

  double rep_intercept, rep_slope; 
  if (dependent_rev){
    normalised_depth=as<NumericVector>(settings_list["normalised.depth"]); 
    rep_intercept=as<double>(settings_list["rep.intercept"]) ; 
    rep_slope=as<double>(settings_list["rep.slope"]); 
  }

  NumericVector beta(num_loci); // coefficients
  NumericVector standard_errors(num_loci); 
  NumericVector means(num_loci); 
  NumericVector lower_bounds(num_loci); // lower bound on the log likehood per locus
  
  vector<NumericVector> a; // part of the variational approximation for the logistic factor
  // natural parameter of q(g), where g is the logistic regression auxiliary variable
  vector<NumericVector> g_mean_prec, mean_prec_f; 
  vector<NumericVector> g_prec, prec_f; 

  vector<NumericVector> flips; 

  // C++ structure
  vector<NumericVector> alt; // counts of alternative SNP
  vector<NumericVector> n; // total reads mappable to either allele at each locus for each individual
  vector<NumericVector> x; // covariate
  
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
    for (int sample_index=0;sample_index<num_samples;sample_index++){
      g_mean_prec[locus_index][sample_index]=0.0; 
      g_prec[locus_index][sample_index]=1.0;
      a[locus_index][sample_index]=.5; 
      flips[locus_index][sample_index]=1.0; 
    }
  }
  double previous_it_lb=-1.0e30;   
  for (int it=0;it<max_its;it++){

    NumericVector expected_err(num_loci);
    NumericVector total_terms(num_loci); 
    double old_lb; 
    for (int locus_index=0;locus_index<num_loci;locus_index++){
      double x_g=0.0, g_1=0.0; 
      double xx=0.0,x1=0.0;
      int num_samples=alt[locus_index].size(); 

      if (debug) old_lb=per_locus_bound(beta[locus_index],means[locus_index], a[locus_index], g_mean_prec[locus_index], g_prec[locus_index], alt[locus_index], n[locus_index], x[locus_index], random_effect_precision, flips[locus_index]); 
      for (int sample_index=0;sample_index<num_samples;sample_index++){
	double x_ns=x[locus_index][sample_index];
	double local_rep= dependent_rev ? exp( rep_slope * normalised_depth[locus_index] + rep_intercept ) : random_effect_precision; 
	double regression_mean=(beta[locus_index]*x_ns+means[locus_index])*flips[locus_index][sample_index]; 
	double old_local_bound=local_bound(regression_mean, g_prec[locus_index][sample_index], g_mean_prec[locus_index][sample_index], a[locus_index][sample_index], local_rep, alt[locus_index][sample_index], n[locus_index][sample_index]); 
	double v=1.0/g_prec[locus_index][sample_index]; 
	double m=g_mean_prec[locus_index][sample_index]/g_prec[locus_index][sample_index]; 
	double a_ns=a[locus_index][sample_index]; 
	double sig=logistic(m+(1.0-2.0*a_ns)*v*.5);
	double old_prec_f=prec_f[locus_index][sample_index]; 
	double old_mean_prec_f=mean_prec_f[locus_index][sample_index]; 
	
	double pf=n[locus_index][sample_index]*sig*(1.0-sig); 
	double mpf=m*pf+alt[locus_index][sample_index]-n[locus_index][sample_index]*sig; 

	// update q(g)
	g_prec[locus_index][sample_index]=local_rep+pf;
	g_mean_prec[locus_index][sample_index]=regression_mean*local_rep+mpf; 
	
	double new_local_bound=local_bound(regression_mean, g_prec[locus_index][sample_index], g_mean_prec[locus_index][sample_index], a[locus_index][sample_index], local_rep, alt[locus_index][sample_index], n[locus_index][sample_index]); 
	double step=1.0; 
	while (new_local_bound<old_local_bound){
	  step *= .5; 
	  g_prec[locus_index][sample_index]=local_rep+step*pf+(1.0-step)*old_prec_f;
	  g_mean_prec[locus_index][sample_index]=regression_mean*local_rep+step*mpf+(1.0-step)*old_mean_prec_f;
	  new_local_bound=local_bound(regression_mean, g_prec[locus_index][sample_index], g_mean_prec[locus_index][sample_index], a[locus_index][sample_index], local_rep, alt[locus_index][sample_index], n[locus_index][sample_index]); 
	  if (step<0.000001){
	    //cout << "Step is very small" << endl;
	    g_prec[locus_index][sample_index]=local_rep+old_prec_f;
	    g_mean_prec[locus_index][sample_index]=regression_mean*local_rep+old_mean_prec_f;
	    new_local_bound=local_bound(regression_mean, g_prec[locus_index][sample_index], g_mean_prec[locus_index][sample_index], a[locus_index][sample_index], local_rep, alt[locus_index][sample_index], n[locus_index][sample_index]); 
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
	  double flipOn = local_bound(-pure_reg_mean,g_prec[locus_index][sample_index], g_mean_prec[locus_index][sample_index], a[locus_index][sample_index], local_rep, alt[locus_index][sample_index], n[locus_index][sample_index]); 
	  double flipOff = local_bound(pure_reg_mean,g_prec[locus_index][sample_index], g_mean_prec[locus_index][sample_index], a[locus_index][sample_index], local_rep, alt[locus_index][sample_index], n[locus_index][sample_index]); 
	  flips[locus_index][sample_index]=(flipOn > flipOff) ? -1.0 : 1.0; 
	}
	x_g+=x_ns*m*flips[locus_index][sample_index];
	g_1+=m*flips[locus_index][sample_index]; 
	xx+=(x_ns*x_ns); // no flips here since -1*-1=1
	x1+=x_ns*flips[locus_index][sample_index];
      }
      if (debug){
	double new_lb=per_locus_bound(beta[locus_index],means[locus_index], a[locus_index], g_mean_prec[locus_index], g_prec[locus_index], alt[locus_index], n[locus_index], x[locus_index], random_effect_precision, flips[locus_index]); 
	if ((new_lb+1.0e-4)<old_lb) 
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
	double new_lb=per_locus_bound(beta[locus_index],means[locus_index], a[locus_index], g_mean_prec[locus_index], g_prec[locus_index], alt[locus_index], n[locus_index], x[locus_index], random_effect_precision, flips[locus_index]); 
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
    }

    // update random_effect_precision

    double old_rep = random_effect_precision; 
    if (debug) old_lb=lower_bound(beta,means,a,g_mean_prec,g_prec,alt,n,x,random_effect_precision,flips);
    if (learn_rep){
      if (dependent_rev){
	NumericMatrix hess(2,2);
	NumericVector grad = lower_bound_grad(rep_slope,rep_intercept,normalised_depth,expected_err,total_terms,hess); 
	// double eps=.0001;
	// double ngrad = ( lower_bound(beta,means,a,g_mean_prec,g_prec,alt,n,x,rep_slope+eps,rep_intercept,normalised_depth,flips) - lower_bound(beta,means,a,g_mean_prec,g_prec,alt,n,x,rep_slope-eps,rep_intercept,normalised_depth,flips) ) / (2.0 * eps) ; 
	// cout << "ngrad " << ngrad << " my grad: " << grad[0] << endl; 
	// ngrad = ( lower_bound(beta,means,a,g_mean_prec,g_prec,alt,n,x,rep_slope,rep_intercept+eps,normalised_depth,flips) - lower_bound(beta,means,a,g_mean_prec,g_prec,alt,n,x,rep_slope,rep_intercept-eps,normalised_depth,flips) ) / (2.0 * eps) ; 
	// cout << "ngrad " << ngrad << " my grad: " << grad[1] << endl; 
	NumericVector dir = twoByTwoSolve(hess,grad); 
	double step=1.0 ;
	NumericVector current = NumericVector::create( rep_slope, rep_intercept ); 
	double current_lb=lower_bound(beta,means,a,g_mean_prec,g_prec,alt,n,x,current[0],current[1],normalised_depth,flips);
	NumericVector proposed = current + step * dir; 
	while ( lower_bound(beta,means,a,g_mean_prec,g_prec,alt,n,x,proposed[0],proposed[1],normalised_depth,flips) < current_lb ){
	  step *= .5;
	  proposed = current + step * dir;
	  if (step < 0.00001){
	    cout << "Warning: step size very small in learning REP model"<< endl; 
	    break; 
	  }
	}
	//	cout << step << endl;
	rep_slope=proposed[0]; 
	rep_intercept=proposed[1]; 
      } else {
	random_effect_precision=accumulate(total_terms.begin(),total_terms.end(),0.0)/accumulate(expected_err.begin(),expected_err.end(),0.0); 
      }
    }
    if (!isfinite(random_effect_precision)) { cerr << "err " << expected_err << " n " << total_terms << endl; throw 1; }
    for (int locus_index=0;locus_index<num_loci;locus_index++){
      double local_rep= dependent_rev ? exp( rep_slope * normalised_depth[locus_index] + rep_intercept ) : random_effect_precision;       
      lower_bounds[locus_index]=per_locus_bound(beta[locus_index],means[locus_index], a[locus_index], g_mean_prec[locus_index], g_prec[locus_index], alt[locus_index], n[locus_index], x[locus_index], local_rep, flips[locus_index]); 
    }
    double lb=accumulate(lower_bounds.begin(),lower_bounds.end(),0.0); 
    if (debug && (lb<old_lb)) cout << "Warning: lb got worse after optimizing random effect var" << endl; 
    double beta_l2 = inner_product(beta.begin(),beta.end(),beta.begin(),0.0); 
    double means_l2 = inner_product(means.begin(),means.end(),means.begin(),0.0); 
    if (trace){ 
      cout << "it: " << it; 
      if (dependent_rev)
	cout << " rep_slope: " << rep_slope << " rep_intercept: " << rep_intercept;
      else
	cout << " re_var: " << 1.0/random_effect_precision;
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
		      _("random.effect.variance") = 1.0/random_effect_precision,
		      _("log.likelihoods") = lower_bounds);
}

double per_locus_bound_null(double mean, NumericVector &a, NumericVector &g_mean_prec, NumericVector &g_prec, NumericVector &alt, NumericVector &n, double random_effect_precision){
  double lb=0.0; 
  for (int sample_index=0;sample_index<n.size();sample_index++){
      lb+=local_bound(mean, g_prec[sample_index], g_mean_prec[sample_index], a[sample_index], random_effect_precision, alt[sample_index], n[sample_index]); 
  }
  return lb; 
}

double lower_bound_null(NumericVector &means, vector<NumericVector> &a, vector<NumericVector> &g_mean_prec, vector<NumericVector> &g_prec, vector<NumericVector> &alt, vector<NumericVector> &n, double random_effect_precision) {
  int num_loci = n.size(); 
  double lb=0.0; 
  for (int locus_index=0;locus_index<num_loci;locus_index++){
    lb+=per_locus_bound_null(means[locus_index], a[locus_index], g_mean_prec[locus_index], g_prec[locus_index], alt[locus_index], n[locus_index], random_effect_precision); 
  }
  return lb; 
}    


RcppExport SEXP fitnull(SEXP alt_sexp, SEXP n_sexp, SEXP max_its_sexp, SEXP tol_sexp, SEXP rev_sexp, SEXP debug_sexp) {
  
  bool debug=as<bool>(debug_sexp); 
  double tol = as<double>(tol_sexp) ; 
  List alt_rlist(alt_sexp);
  List n_rlist(n_sexp);
  
  int num_loci = alt_rlist.size(); 
  int max_its = as<int>(max_its_sexp); 

  cout << "VBGLMM fit null: num_loci: " << num_loci << " max_its: " << max_its << endl; 

  NumericVector means(num_loci); 
  NumericVector lower_bounds(num_loci); 
  
  vector<NumericVector> a; 
  vector<NumericVector> g_mean_prec; 
  vector<NumericVector> g_prec; 

  vector<NumericVector> alt; 
  vector<NumericVector> n;

  double random_effect_precision=1.0/as<double>(rev_sexp);

  for (int locus_index=0;locus_index<num_loci;locus_index++){
    NumericVector alti((SEXP)alt_rlist[locus_index]); 
    alt.push_back(alti); 
    NumericVector ni((SEXP)n_rlist[locus_index]); 
    n.push_back(ni); 
    int num_samples=alti.size(); 
    NumericVector ai(num_samples,0.5); 
    a.push_back(ai);
    NumericVector mp(num_samples,0.0) ;
    g_mean_prec.push_back(mp); 
    NumericVector p(num_samples,random_effect_precision); 
    g_prec.push_back(p); 

    /*    for (int sample_index=0;sample_index<num_samples;sample_index++){
      g_mean_prec[locus_index][sample_index]=0.0; 
      g_prec[locus_index][sample_index]=1.0;
      a[locus_index][sample_index]=.5; 
      }*/
  }
  double previous_it_lb=-1.0e30; 
  for (int it=0;it<max_its;it++){
    double old_lb; 
    for (int locus_index=0;locus_index<num_loci;locus_index++){
      double g_1=0.0; 
      int num_samples=alt[locus_index].size(); 

      if (debug) old_lb=per_locus_bound_null(means[locus_index], a[locus_index], g_mean_prec[locus_index], g_prec[locus_index], alt[locus_index], n[locus_index], random_effect_precision); 
      for (int sample_index=0;sample_index<num_samples;sample_index++){
	double regression_mean=means[locus_index];
	double old_local_bound=local_bound(regression_mean, g_prec[locus_index][sample_index], g_mean_prec[locus_index][sample_index], a[locus_index][sample_index], random_effect_precision, alt[locus_index][sample_index], n[locus_index][sample_index]); 
	double v=1.0/g_prec[locus_index][sample_index]; 
	double m=g_mean_prec[locus_index][sample_index]/g_prec[locus_index][sample_index]; 
	double a_ns=a[locus_index][sample_index]; 
	double sig=logistic(m+(1.0-2.0*a_ns)*v*.5);
	double prec_f=n[locus_index][sample_index]*sig*(1.0-sig); 
	
	double mean_prec_f=m*prec_f+alt[locus_index][sample_index]-n[locus_index][sample_index]*sig; 

	// update q(g)
	double old_prec_f=g_prec[locus_index][sample_index]-random_effect_precision; 
	double old_mean_prec_f=g_mean_prec[locus_index][sample_index]-regression_mean*random_effect_precision; 
	g_prec[locus_index][sample_index]=random_effect_precision+prec_f;
	g_mean_prec[locus_index][sample_index]=regression_mean*random_effect_precision+mean_prec_f; 
	
	double new_local_bound=local_bound(regression_mean, g_prec[locus_index][sample_index], g_mean_prec[locus_index][sample_index], a[locus_index][sample_index], random_effect_precision, alt[locus_index][sample_index], n[locus_index][sample_index]); 
	double step=1.0; 
	while (new_local_bound<old_local_bound){
	  step *= .5; 
	  g_prec[locus_index][sample_index]=random_effect_precision+step*prec_f+(1.0-step)*old_prec_f;
	  g_mean_prec[locus_index][sample_index]=regression_mean*random_effect_precision+step*mean_prec_f+(1.0-step)*old_mean_prec_f;
	  new_local_bound=local_bound(regression_mean, g_prec[locus_index][sample_index], g_mean_prec[locus_index][sample_index], a[locus_index][sample_index], random_effect_precision, alt[locus_index][sample_index], n[locus_index][sample_index]); 
	  if (step<0.000001){
	    //cout << "Step is very small" << endl; 
	    g_prec[locus_index][sample_index]=random_effect_precision+old_prec_f;
	    g_mean_prec[locus_index][sample_index]=regression_mean*random_effect_precision+old_mean_prec_f;
	    new_local_bound=local_bound(regression_mean, g_prec[locus_index][sample_index], g_mean_prec[locus_index][sample_index], a[locus_index][sample_index], random_effect_precision, alt[locus_index][sample_index], n[locus_index][sample_index]); 
	    break; 
	  }
	} 
	  
	//if (new_local_bound<old_local_bound) cout << "GDI bound got worse, old: " << old_local_bound << " new: " << new_local_bound << endl; 
	m=g_mean_prec[locus_index][sample_index]/g_prec[locus_index][sample_index]; 
	v=1.0/g_prec[locus_index][sample_index]; 
	// update a
	//double check=.5*a_ns*a_ns*v+log1plusexp(m+(1.0-2.0*a_ns)*v*.5); 
	a[locus_index][sample_index]=logistic(m+(1.0-2.0*a_ns)*v*.5);
	a_ns=a[locus_index][sample_index]; 
	//double check2=.5*a_ns*a_ns*v+log1plusexp(m+(1.0-2.0*a_ns)*v*.5); 
	//if (check2>check) cout << "update of a was bad: before " << check << " after " << check2 << endl;
	g_1+=m; 
      }
      if (debug){
	double new_lb=per_locus_bound_null(means[locus_index], a[locus_index], g_mean_prec[locus_index], g_prec[locus_index], alt[locus_index], n[locus_index],  random_effect_precision); 
	if ((new_lb+1.0e-4)<old_lb) 
	  cout << "Warning: lb got worse after optimizing q(g) and a, old:" << old_lb << " new: " << new_lb << endl;   	
	old_lb=new_lb; 
      }
      // update beta and means
      double one1=(double)num_samples;
      means[locus_index]= (one1 > 0.0) ? g_1/one1 : 0.0; 

    }
    // update random_effect_precision
    double means_l2 = inner_product(means.begin(),means.end(),means.begin(),0.0); 
    for (int locus_index=0;locus_index<num_loci;locus_index++)
      lower_bounds[locus_index]=per_locus_bound_null(means[locus_index], a[locus_index], g_mean_prec[locus_index], g_prec[locus_index], alt[locus_index], n[locus_index], random_effect_precision); 
    double lb=accumulate(lower_bounds.begin(),lower_bounds.end(),0.0); 
    cout << "it: " << it << " re_var: " << 1.0/random_effect_precision << " mean l2 " << means_l2 << " lb: " << lb << endl; 
    R_CheckUserInterrupt(); 
    if (abs(lb-previous_it_lb) < tol){
      cout << "Converged!" << endl; 
      break; 
    }
    previous_it_lb=lb; 
  }
        
  return List::create(_("log.likelihoods") = lower_bounds, _("means")=means);
}

