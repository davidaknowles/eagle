#include <Eigen/Dense>
#include "vbglmm.h"
#include <vector> 
#include <math.h> 
#include <assert.h> 
#include <cmath> 

#define REV_GLOBAL 0
#define REV_REGRESSION 1
#define REV_LOCAL 2
#define REV_LOCAL_REGRESSION 3

#define FLIPS_NONE 0
#define FLIPS_HARD 1 
#define FLIPS_SOFT 2
#define FLIPS_STRUCTURED 3

using namespace std ;
using namespace Rcpp ;
using namespace Eigen ;

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

NumericMatrix outer(NumericVector &x){
  int n=x.size(); 
  NumericMatrix res(n,n); 
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      res(i,j)=x[i]*x[j];
  return res; 
}

class G_stuff {
public : 
  vector<NumericVector> a; // part of the variational approximation for the logistic factor

  // natural parameter of q(g), where g is the logistic regression auxiliary variable
  vector<NumericVector> mean_prec; 
  vector<NumericVector> prec; 
  
  void init(vector<NumericVector> &alt, NumericVector &precisions){
    int num_loci=alt.size();
    for (int locus_index=0;locus_index<num_loci;locus_index++){
      int num_samples=alt[locus_index].size(); 
      NumericVector ai(num_samples,0.5); 
      a.push_back(ai);
      NumericVector mp(num_samples,0.0) ;
      mean_prec.push_back(mp); 
      NumericVector p(num_samples,precisions[locus_index]); 
      prec.push_back(p); 
    }
  }

};

class Vbglmm {
  bool learn_betas; 
  int it, its_to_hold_flips_fixed; 
  bool debug;
  int num_loci;
  int max_its; 
  double converge_tol; 
  bool learn_rep; 
  bool trace; 
  int rev_model; 
  int flips_setting; 
  bool learn_flips_prior; 
  bool return_aux_variables;
  bool store_all_coeffs; 
  double coeff_regulariser;
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

  vector<NumericVector> beta; // coefficients
  NumericVector standard_errors; 
  NumericVector lower_bounds; // lower bound on the log likehood per locus
  
  G_stuff g_stuff; 
  G_stuff g_flip; 

  vector<NumericVector> flips; 
  vector<NumericVector> flips_log_odds; 
  double flips_log_odds_prior, flips_log1minusP; 
  vector<NumericVector> to_flip; 

  // C++ structure
  vector<NumericVector> alt; // counts of alternative SNP
  vector<NumericVector> n; // total reads mappable to either allele at each locus for each individual
  //vector<vector<NumericVector> > x; // covariate
  vector<NumericMatrix> x; 

  NumericVector selective_flip(const NumericVector &x_ns, const NumericVector &index, double flip) {
    NumericVector result(x_ns.size()); 
    for (int i=0;i<x_ns.size();i++)
      result[i] = (index[i]==1.0) ? x_ns[i]*flip : x_ns[i]; 
    return result; 
  }

  double local_bound(int locus_index, int sample_index, double pred_times_flip, double Erep, double Elogrep, G_stuff &g, double pred_to_flip){
    double v=1.0/g.prec[locus_index][sample_index]; 
    double m=g.mean_prec[locus_index][sample_index]*v; 
    double Elog1plus_exp_g=.5*g.a[locus_index][sample_index]*g.a[locus_index][sample_index]*v+log1plusexp(m+(1.0-2.0*g.a[locus_index][sample_index])*v*.5); 
    // < log likelihood >
    double lb=(double)alt[locus_index][sample_index]*m-(double)n[locus_index][sample_index]*Elog1plus_exp_g; 
    double err = pred_times_flip-m; 
    double Eerr2 = err*err+v;
    if (flips_setting == FLIPS_SOFT){
      double flip_prob=.5*(1.0-flips[locus_index][sample_index]);
      Eerr2+=pred_to_flip*pred_to_flip*4.0*flip_prob*(1.0-flip_prob); 
    }
    // < log normal >  [note: cancelling log2pi terms]
    lb+=-.5*Erep*Eerr2+.5*Elogrep; 
    // - <log q>
    lb+=.5*(v+log(v)); 
    if (isnan(lb)) throw 1; 
    return lb; 
  }

  double local_bound(int locus_index, int sample_index, double pred_no_flip, double Erep, double Elogrep, double pred_to_flip){
    double lb=0.0;
    if (flips_setting == FLIPS_SOFT  || flips_setting == FLIPS_STRUCTURED){
      double flip_prob=.5*(1.0-flips[locus_index][sample_index]);
      // entropy TODO: more numerically stable to use log odds? 
      lb -= (flip_prob==0.0 || flip_prob==1.0) ? 0.0 : (flip_prob*log(flip_prob)+(1.0-flip_prob)*log(1.0-flip_prob) ); 
      //if (isnan(lb)) { cout << flip_prob << " " << flips_log_odds << endl; throw 1; }
      lb += flips_log1minusP + flip_prob * flips_log_odds_prior;
    } 
    if (flips_setting == FLIPS_HARD)
      lb += ((flips[locus_index][sample_index]==-1.0) ? flips_log_odds_prior : 0.0) + flips_log1minusP ; 
    if (flips_setting == FLIPS_STRUCTURED){
      double flip_prob=.5*(1.0-flips[locus_index][sample_index]);
      double lb_flipped = local_bound(locus_index, sample_index, pred_no_flip - pred_to_flip, Erep, Elogrep, g_flip, pred_to_flip); 
      double lb_not = local_bound(locus_index, sample_index, pred_no_flip + pred_to_flip, Erep, Elogrep, g_stuff, pred_to_flip); 
      lb += flip_prob * lb_flipped + (1.0 - flip_prob) * lb_not; 
      if (isnan(lb)) throw 1; 
      return lb; 
    } else {
      return lb + local_bound(locus_index, sample_index, pred_to_flip * flips[locus_index][sample_index] + pred_no_flip, Erep, Elogrep, g_stuff, pred_to_flip); 
    }
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
      lb+=(global_rep_shape-1.0)*Elogrep + global_rep_rate * Erep + global_rep_shape * log(global_rep_rate) - ::Rf_lgammafn(global_rep_shape); 
      // - < log q >
      lb-= (rep_shapes[locus_index]-1.0)*::Rf_digamma(rep_shapes[locus_index])+log(rep_rates[locus_index])-rep_shapes[locus_index]-::Rf_lgammafn(rep_shapes[locus_index]); 
      break;
    case REV_LOCAL_REGRESSION:
      double v=1.0/logrep_p[locus_index] ;
      double m=logrep_mp[locus_index]*v; 
      Erep = exp(m+.5*v); 
      Elogrep = m; 
      double rpred=rep_slope*normalised_depth[locus_index]+rep_intercept; 
      double err=rpred-m; 
      double Eerr2 = err*err+v; 
      // < log N(log local_rep; mx+c,v) > [cancelling log2pi]
      lb+=-.5*rep_rep*Eerr2+.5*log(rep_rep); 
      // - <log q>
      lb+=.5*(v+log(v)); 
      break;
    }

    for (int sample_index=0;sample_index<alt[locus_index].size();sample_index++){
      NumericVector x_ns=x[locus_index]( sample_index, _ );
      double pred_to_flip=sum(beta[locus_index]*x_ns*to_flip[locus_index]); 
      double pred_no_flip=sum(beta[locus_index]*x_ns*(-to_flip[locus_index]+1.0));
      lb+=local_bound(locus_index, sample_index, pred_no_flip, Erep, Elogrep,pred_to_flip); 
    }
    return lb; 
  }

  double lower_bound(){
    int num_loi = n.size(); 
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

  
  double update_single_g(int locus_index, int sample_index, double local_rep, double log_local_rep, G_stuff &g){
    NumericVector x_ns=x[locus_index]( sample_index, _ );
    double pred_to_flip=sum(beta[locus_index]*x_ns*to_flip[locus_index]); 
    double regression_mean=pred_to_flip*flips[locus_index][sample_index]+sum(beta[locus_index]*x_ns*(-to_flip[locus_index]+1.0)); 
    double v=1.0/g.prec[locus_index][sample_index]; 
    double m=g.mean_prec[locus_index][sample_index]/g.prec[locus_index][sample_index]; 
    double a_ns=g.a[locus_index][sample_index]; 
    double sig=logistic(m+(1.0-2.0*a_ns)*v*.5);
    double old_prec=g.prec[locus_index][sample_index]; 
    double old_mean_prec=g.mean_prec[locus_index][sample_index];
    double pf=n[locus_index][sample_index]*sig*(1.0-sig); 
    double p=pf+local_rep;
    double mp=m*pf+alt[locus_index][sample_index]-n[locus_index][sample_index]*sig+local_rep*regression_mean; 
    double new_local_bound; 
    double old_local_bound=local_bound(locus_index, sample_index, regression_mean, local_rep, log_local_rep, g, pred_to_flip); 
    // update q(g) [NB: this is still correct for random q(rep)]
    g.prec[locus_index][sample_index]=p;
    g.mean_prec[locus_index][sample_index]=mp;

    
    new_local_bound=local_bound(locus_index, sample_index, regression_mean, local_rep, log_local_rep, g, pred_to_flip); 
    double step=1.0; 
    while (new_local_bound<old_local_bound){
      step *= .5; 
      g.prec[locus_index][sample_index]=step*p+(1.0-step)*old_prec;
      g.mean_prec[locus_index][sample_index]=step*mp+(1.0-step)*old_mean_prec;
      new_local_bound=local_bound(locus_index, sample_index, regression_mean, local_rep, log_local_rep, g, pred_to_flip); 
      if (step<0.000001){
	if (debug) cout << "Step is very small, reg mean: " << regression_mean << endl;
	g.prec[locus_index][sample_index]=old_prec; 
	g.mean_prec[locus_index][sample_index]=old_mean_prec; 
	new_local_bound=local_bound(locus_index, sample_index, regression_mean, local_rep, log_local_rep, g, pred_to_flip); 
	break; 
      }
    }
    //if (new_local_bound<old_local_bound) cout << "GDI bound got worse, old: " << old_local_bound << " new: " << new_local_bound << endl; 
    m=g.mean_prec[locus_index][sample_index]/g.prec[locus_index][sample_index]; 
    v=1.0/g.prec[locus_index][sample_index]; 
    // update a
    //double check=.5*a_ns*a_ns*v+log1plusexp(m+(1.0-2.0*a_ns)*v*.5); 
    g.a[locus_index][sample_index]=logistic(m+(1.0-2.0*a_ns)*v*.5);

    double old_lb = new_local_bound; 
    new_local_bound=local_bound(locus_index, sample_index, regression_mean, local_rep, log_local_rep, g, pred_to_flip); 

    // if (debug) {
    // if ((new_local_bound + 0.001) < old_lb)
    //	cout << "Warning: updating a decreased lower bound, old: " << old_lb << " new: " << new_local_bound << endl;
    //} 
    return new_local_bound; 
    //double check2=.5*a_ns*a_ns*v+log1plusexp(m+(1.0-2.0*a_ns)*v*.5); 
    //if (check2>check) cout << "update of a was bad: before " << check << " after " << check2 << endl;
  
  }

  double converge_local(int locus_index, int sample_index, double local_rep, double log_local_rep){
    NumericVector x_ns=x[locus_index]( sample_index, _ );
    g_stuff.prec[locus_index][sample_index]=local_rep; 
    g_stuff.mean_prec[locus_index][sample_index]=0.0 ;
    g_stuff.a[locus_index][sample_index]=0.5;
    double pred_to_flip=sum(beta[locus_index]*x_ns*to_flip[locus_index]); 
    double regression_mean1=pred_to_flip*flips[locus_index][sample_index]+sum(beta[locus_index]*x_ns*(-to_flip[locus_index]+1.0)); 
    int max_inner=10; 
    double threshold=0.1; 
    double old_lb=local_bound(locus_index, sample_index, regression_mean1, local_rep, log_local_rep, pred_to_flip); 
    double new_lb;
    for (int ii=0;ii<max_inner;ii++){
      update_single_g(locus_index, sample_index, local_rep, log_local_rep, g_stuff); 
      new_lb=local_bound(locus_index, sample_index, regression_mean1, local_rep, log_local_rep, pred_to_flip); 
      if (abs(new_lb-old_lb)<threshold)
	break; 
      old_lb=new_lb; 
    }
    return new_lb; 
  }

  double update_flip(int locus_index, int sample_index, double local_rep, double log_local_rep){
    NumericVector x_ns=x[locus_index]( sample_index, _ );
    double pred_to_flip=sum(beta[locus_index]*x_ns*to_flip[locus_index]); 
    double pred_no_flip=sum(beta[locus_index]*x_ns*(-to_flip[locus_index]+1.0));
    double old_lb; 
    if (debug) {
      old_lb=local_bound(locus_index, sample_index, pred_no_flip, local_rep, log_local_rep, pred_to_flip); 
    }
      
    if (flips_setting==FLIPS_SOFT || flips_setting==FLIPS_HARD){
      double regression_mean=pred_to_flip*flips[locus_index][sample_index]+pred_no_flip; 
      double m=g_stuff.mean_prec[locus_index][sample_index]/g_stuff.prec[locus_index][sample_index]; 
      flips_log_odds[locus_index][sample_index]=flips_log_odds_prior - 2.0*local_rep*m*regression_mean;
      flips[locus_index][sample_index]=1.0-2.0*logistic(flips_log_odds[locus_index][sample_index]);  
      double mmax=.999; 
      if ((abs(flips[locus_index][sample_index])>mmax) || isnan(flips[locus_index][sample_index]))
	{
	  flips[locus_index][sample_index]=min(flips[locus_index][sample_index],mmax); 
	  flips[locus_index][sample_index]=max(flips[locus_index][sample_index],-mmax); 
	}
      if (flips_setting == FLIPS_HARD){
	// flips[locus_index][sample_index]= (flips[locus_index][sample_index] > 0.0) ? 1.0 : -1.0; 
	flips[locus_index][sample_index]=1.0; 
	double noflip_lb=converge_local(locus_index,sample_index,local_rep,log_local_rep); // TODO start from scratch? 
	flips[locus_index][sample_index]=-1.0; 
	double flip_lb=converge_local(locus_index,sample_index,local_rep,log_local_rep); 
	flips[locus_index][sample_index]=((flip_lb - noflip_lb) > 0.0) ? -1.0 : 1.0; // NOTE: prior included in bound now
	converge_local(locus_index,sample_index,local_rep,log_local_rep); // TODO save result rather than recomputing
      }
    } else if (flips_setting == FLIPS_STRUCTURED){
        double old_flips=flips[locus_index][sample_index]; 
	flips[locus_index][sample_index]=1.0; 
	double noflip_lb=update_single_g(locus_index,sample_index,local_rep,log_local_rep, g_stuff);
	flips[locus_index][sample_index]=-1.0; 
	double flip_lb=update_single_g(locus_index,sample_index,local_rep,log_local_rep, g_flip);       
	if (it < its_to_hold_flips_fixed){
	  flips[locus_index][sample_index]=old_flips; 
	} else {
	  //cout << "ll flip: " << flip_lb << " ll no flip " << noflip_lb <<endl; 
	  flips_log_odds[locus_index][sample_index]=flip_lb + flips_log_odds_prior - noflip_lb;
	  flips[locus_index][sample_index]=1.0-2.0*logistic(flips_log_odds[locus_index][sample_index]);  
	}
    }
    if (debug){
      double new_local_bound=local_bound(locus_index, sample_index, pred_no_flip, local_rep, log_local_rep, pred_to_flip); 
      if ((new_local_bound+0.001) < old_lb)
	cout << "Warning: updating flips decreased lower bound, old: " << old_lb << " new: " << new_local_bound << endl;
    }

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
      break; 
    case REV_LOCAL_REGRESSION:
      double v=1.0/logrep_p[locus_index] ;
      double m=logrep_mp[locus_index]*v; 
      local_rep = exp(m+.5*v); 
      log_local_rep = m; 
      break;
    }
    int num_cov=x[locus_index].ncol(); 
    NumericVector xg(num_cov); 
    NumericMatrix xx(num_cov,num_cov); ;
    int num_samples=alt[locus_index].size(); 

    if (debug) old_lb=per_locus_bound(locus_index); 
    for (int sample_index=0;sample_index<num_samples;sample_index++){
      if (flips_setting != FLIPS_STRUCTURED) // in that case updating q(g) and a is done with update_flip
	update_single_g(locus_index, sample_index, local_rep, log_local_rep, g_stuff); 
      if (flips_setting != FLIPS_NONE) 
	update_flip(locus_index, sample_index, local_rep, log_local_rep); 
      NumericVector x_ns=x[locus_index]( sample_index, _ );
      double out; 
      double m=g_stuff.mean_prec[locus_index][sample_index]/g_stuff.prec[locus_index][sample_index]; 
      if (flips_setting == FLIPS_STRUCTURED){
	double m_flip=g_flip.mean_prec[locus_index][sample_index]/g_flip.prec[locus_index][sample_index]; 
	double flip_prob=.5*(1.0-flips[locus_index][sample_index]);
	NumericVector x_flipped = selective_flip(x_ns, to_flip[locus_index], -1.0); 
	xg += (x_flipped * flip_prob * m_flip) + (x_ns * (1.0-flip_prob) * m); 
	xx += (flip_prob * outer(x_flipped)) + ((1.0 - flip_prob) * outer(x_ns)); 
      } else {
	xg += x_ns * flips[locus_index][sample_index] * m; 
	xx += outer(x_ns);// TODO: cache XX
      }
    }
    if (debug){
      double new_lb=per_locus_bound(locus_index); 
      if ((new_lb+1.0e-3)<old_lb) 
	cout << "Warning: lb got worse after optimizing q(g) and a, old:" << old_lb << " new: " << new_lb << endl;   	
      old_lb=new_lb; 
    }
    if (learn_betas){
      NumericVector old_beta=beta[locus_index];
      // update beta and means
      // TODO: might want to check if determinant is too close to zero
      MatrixXd xxX(num_cov,num_cov);
      for (int ii=0;ii<num_cov;ii++)
	for (int jj=0;jj<num_cov;jj++){
	  xx(ii,jj) += (coeff_regulariser>0.0 && ii==jj) ? coeff_regulariser : 0.0; 
	  xxX(ii,jj)=xx(ii,jj); 
	}
      LLT<MatrixXd> llt(xxX); // TODO cache these solves if slow
      
      VectorXd xgv(num_cov);
      for (int ii=0;ii<beta[locus_index].size();ii++)
	xgv[ii]=xg[ii];
      VectorXd res=llt.solve(xgv);
      // beta[locus_index]=wrap(res);
      
      for (int ii=0;ii<beta[locus_index].size();ii++)
	beta[locus_index][ii]=res[ii]; 
      if (!isfinite(beta[locus_index][0])) throw 1; 
      if (debug && (num_cov==2)){
	NumericVector check = twoByTwoSolve(xx,xg); 
	if (abs(check[0]-beta[locus_index][0])>0.0001) throw 1; 
      }
      //cout << "beta: " << beta[locus_index][0] << "," << beta[locus_index][1] << endl; 
    }
    expected_err[locus_index]=0.0; 
    for (int sample_index=0;sample_index<num_samples;sample_index++){
      double m=g_stuff.mean_prec[locus_index][sample_index]/g_stuff.prec[locus_index][sample_index]; 
      double v=1.0/g_stuff.prec[locus_index][sample_index]; 
      NumericVector x_ns=x[locus_index]( sample_index, _ );
      double pred_to_flip=sum(beta[locus_index]*x_ns*to_flip[locus_index]); 
      double pred_no_flip=sum(beta[locus_index]*x_ns*(-to_flip[locus_index]+1.0)); 
      double pred=pred_to_flip*flips[locus_index][sample_index]+pred_no_flip; 
      if (flips_setting == FLIPS_STRUCTURED){
	double err_noflip = (pred_no_flip + pred_to_flip) - m; 
	err_noflip = err_noflip * err_noflip + v; 
	m=g_flip.mean_prec[locus_index][sample_index]/g_flip.prec[locus_index][sample_index]; 
	v=1.0/g_flip.prec[locus_index][sample_index]; 
	double flip_prob=.5*(1.0-flips[locus_index][sample_index]);
      	double err_flip = (pred_no_flip - pred_to_flip) - m; 
	err_flip = err_flip * err_flip + v; 
	expected_err[locus_index]+=flip_prob * err_flip + (1.0-flip_prob) * err_noflip; 
      } else {
	double err=pred*flips[locus_index][sample_index]-m;
	if (!isfinite(err)) throw 1;
	expected_err[locus_index]+=err*err+v; 
	if (flips_setting==FLIPS_SOFT){
	  double p=.5*(1.0-flips[locus_index][sample_index]); 
	  expected_err[locus_index]+=pred*pred*4.0*p*(1.0-p); 
	}
      }
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
	      cout << "Warning: step size very small in learning REP model"<< endl; 
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
	// cout << "Eerr2: " << Eerr2 << endl; 
	for (int locus_index=0;locus_index<num_loci;locus_index++){
	  double pred_mean=rep_slope*normalised_depth[locus_index]+rep_intercept; 
	  logrep_mp[locus_index]=pred_mean*rep_rep+delta_to_logrep_mp[locus_index]; 
	  logrep_p[locus_index]=rep_rep+delta_to_logrep_p[locus_index]; 
	}
	lb=lower_bound(); 

      }
      
      if (!isfinite(random_effect_precision)) { cerr << "err " << expected_err << " n " << total_terms << endl; throw 1; }
      
      if (debug && ((lb+1.0e-3)<old_lb)) cout << "Warning: lb got worse after optimizing random effect var, old:" << old_lb << " new:" << lb << endl; 
    } else {
      lb=lower_bound(); 
    }

    if (learn_flips_prior){ // note: this should work for hard assignment too
      double sumQ=0.0,total=0.0,minf=0.0, maxf=0.0;
      for (int locus_index=0;locus_index<num_loci;locus_index++){
	sumQ += sum(flips[locus_index]); 
	total += (double)flips[locus_index].size(); 
      }
      sumQ = .5*(total-sumQ); 
      flips_log_odds_prior = log(sumQ+1.0) - log(1.0+total - sumQ);
      if (!isfinite(flips_log_odds_prior)) throw 1; 
      flips_log1minusP = log(1.0 - logistic(flips_log_odds_prior)); 
    }

    return lb; 
  }

public:

  SEXP run(){
    double previous_it_lb=-1.0e30;
    vector<double> mlPerIteration, repInterceptPerIteration, repSlopePerIt, repRepPerIt; 
    List allCoefs; 
    for (it=0;it<max_its;it++){
      double lb = one_iteration(); 
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
	if (flips_setting != FLIPS_NONE)
	  cout << " flips_prior_logodds " << flips_log_odds_prior ; 
	cout << " lb: " << lb << endl; 
      }
      R_CheckUserInterrupt(); 
      if (abs(lb-previous_it_lb) < converge_tol){
	if (trace) cout << "Converged!" << endl; 
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
    }

    List result_list=List::create(_("coeffs") = beta,
				  _("mlPerIteration") = mlPerIteration, 
				  _("repInterceptPerIteration")=repInterceptPerIteration,
				  _("repSlopePerIt")=repSlopePerIt, 
				  _("repRepPerIt")=repRepPerIt, 
				  _("log.likelihoods") = lower_bounds);
    if (flips_setting != FLIPS_NONE){
      result_list["flips.logodds.prior"] = flips_log_odds_prior;
      result_list["flips"]=flips; 
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
    case REV_REGRESSION:
      result_list["rep.slope"] = rep_slope;
      result_list["rep.intercept"] = rep_intercept;
      break;
    }
    
    if (return_aux_variables){
      result_list["aux.variables.mp"]=g_stuff.mean_prec;
      result_list["aux.variables.p"]=g_stuff.prec;
    }
    return result_list; 
  }

  Vbglmm(SEXP alt_sexp, SEXP n_sexp, SEXP x_sexp, SEXP settings_sexp){

    // convert R to Rcpp
    List alt_rlist(alt_sexp); 
    List n_rlist(n_sexp);
    List x_rlist(x_sexp);
    List settings_list(settings_sexp);
    learn_betas=as<bool>(settings_list["learnBetas"]);
    debug=as<bool>(settings_list["debug"]);
    num_loci = alt_rlist.size();
    max_its = as<int>(settings_list["max.iterations"]); 
    converge_tol = as<double>(settings_list["convergence.tolerance"]); 
    learn_rep=as<bool>(settings_list["learn.rev"]); 
    flips_setting=as<int>(settings_list["flips.setting"]);
    learn_flips_prior=as<bool>(settings_list["learn.flips.prior"]); 
    trace=as<bool>(settings_list["trace"]); 
    rev_model=as<int>(settings_list["rev.model"]);
    coeff_regulariser=as<double>(settings_list["coeff.regulariser"]); 
    return_aux_variables=as<bool>(settings_list["return.aux.variables"]); 
    store_all_coeffs=as<bool>(settings_list["storeAllCoeffs"]); 
    bool init_flips=as<bool>(settings_list["init.flips"]); 
    List to_flip_list((SEXP)settings_list["toFlip"]); 
    List flips_rlist; 
    if (init_flips){
      flips_rlist=as<List>(settings_list["flips"]); 
      its_to_hold_flips_fixed=10;
    } else {
      its_to_hold_flips_fixed=0; 
    }
    if (flips_setting != FLIPS_NONE){
      flips_log_odds_prior = as<double>(settings_list["flips.logodds.prior"]); 
      flips_log1minusP = log(1.0-logistic(flips_log_odds_prior)); 
    }
 
    if (trace) cout << "VBGLMM: num_loci: " << num_loci << " max_its: " << max_its << " tolerance: " << converge_tol << endl; 

    NumericVector temp(num_loci); 
    rep_shapes = clone(temp);
    rep_rates = clone(temp); 
    standard_errors=clone(temp); 
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
    NumericVector precs(num_loci); 
    x.resize(num_loci); 
    to_flip.resize(num_loci); 
    for (int locus_index=0;locus_index<num_loci;locus_index++){ 
     if (trace && ((locus_index % 1000)==0))
	cout << "Loading... locus " << locus_index << "/" << num_loci << endl; 
      NumericVector to_flipi((SEXP)to_flip_list[locus_index]); 
      to_flip[locus_index]=as<NumericVector>(to_flipi); 
      NumericVector alti((SEXP)alt_rlist[locus_index]); 
      alt.push_back(alti); 
      NumericVector ni((SEXP)n_rlist[locus_index]); 
      n.push_back(ni); 
      //      NumericMatrix xi(
      int num_samples=alti.size(); 
      x[locus_index]=as<NumericMatrix>((SEXP)x_rlist[locus_index]); 
      int num_cov=x[locus_index].ncol();
      NumericVector b(num_cov,0.0);
      beta.push_back(b); 
      if (init_flips){
	NumericVector flipsi((SEXP)flips_rlist[locus_index]); 
	flips.push_back(flipsi); 
	NumericVector probs=-.5*(flipsi-1.0); 
	NumericVector fl=log(probs)-log(-(probs-1.0)); 
	flips_log_odds.push_back(fl); 
      } else {
	switch (flips_setting){
	case FLIPS_STRUCTURED:
	case FLIPS_SOFT:
	{
	  double init=0.8; 
	  NumericVector f(num_samples,init);
	  double prob=.5*(1.0-init);
	  flips.push_back(f); 
	  NumericVector fl(num_samples,log(prob)-log(1.0-prob)); 
	  flips_log_odds.push_back(fl); 
	}
	break; 
	case FLIPS_HARD:
	case FLIPS_NONE:
	{
	  NumericVector f(num_samples,1.0); 
	  flips.push_back(f); 
	  NumericVector fl(num_samples,0.0); // not used 
	  flips_log_odds.push_back(fl); 
	}
	break;
	}
      }
         
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
    }

    
    g_stuff.init(alt,precs); 
    if (flips_setting == FLIPS_STRUCTURED)
      g_flip.init(alt,precs);
    if (trace) cout << "Loaded" << endl; 
  }
};

  RcppExport SEXP runvb(SEXP alt_sexp, SEXP n_sexp, SEXP x_sexp, SEXP settings_sexp) { 
    Vbglmm vbglmm(alt_sexp, n_sexp, x_sexp, settings_sexp); 
    return vbglmm.run();
  }
