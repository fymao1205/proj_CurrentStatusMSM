// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppNumerical.h>
//#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Rmath.h>
#include <iterator>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace arma;
using namespace Numer;



// M(t) = E(exp(t * X)) = int exp(t * x) * f(x) dx, f(x) is the p.d.f.
class Mintegrand: public Func
{
private:
  const double a;
  const double b;
  const double c;
  const double d;
public:
  Mintegrand(double a_, double b_, double c_, double d_) : a(a_), b(b_), c(c_), d(d_) {}
  
  double operator()(const double& x) const
  {
    return R::pweibull(x, c, d, 0, 0) * R::dweibull(x, a, b, 0);
    //return a*b*std::pow(a*x, b-1) * std::exp(-std::pow(a*x, b)) * std::exp(-std::pow(c*x, d));
  }
};

// [[Rcpp::export]]
double dbwei_mgf(double u, double a, double b, double c, double d)
{
  Mintegrand f(a, b, c, d);
  double err_est;
  int err_code;
  //double res=integrate(f, u1, u2, err_est, err_code);
  return 1-integrate(f, 0.0, u, err_est, err_code);
}

// [[Rcpp::export]]
double dgengam(double u, double mu, double sigma, double Q){
  
  if(Q!=0){
    
    double y=std::log(u), abs_q=std::abs(Q); 
    double w=(y-mu)/sigma, qi=1/(Q*Q);
    double qw=Q*w;
    double res0=-std::log(sigma*u)+std::log(abs_q)*(1-2*qi)+qi*(qw-std::exp(qw))-R::lgammafn(qi);
    return exp(res0); 
    
  }else{
    return R::dlnorm(u, mu, sigma, 0);
  }
  
}

// [[Rcpp::export]]
double hwei_double(double t, String family, NumericVector params0, int logInd){
  
  double res=0.0;
  //int len_cov = X.size();
  if(family=="weibull"){
    
    double y = params0[0]*params0[1]*pow(params0[0]*t, params0[1]-1);
    //res = h12_0*exp(beta1*X1+mu_X2);
    if(logInd){
      res = log(y);
    }else{
      res=y;
    }
  }
  
  return(res);
}

// [[Rcpp::export]]
double survwei_double(double t, String family, NumericVector params0, int logInd){
  
  double res=0.0;
  //int len_cov = X.size();
  if(family=="weibull"){
    //double cov_reg = exp(beta1*X1+mu_X2);
    double H0 = pow(params0[0]*t, params0[1]);
    if(logInd){
      res = -H0;
    }else{
      res = exp(-H0);
    }
    
  }
  
  return(res);
  
}

//[[Rcpp::export()]]
double hpwc_double(double x, NumericVector cuts, NumericVector levels, int logInd)
{
  //cuts and levels are changing with every age of fsex: cutsL=c(L,cuts[>L]); levels=c(0,levels[>L])
  
  NumericVector cut = sort_unique(cuts);
  int p = levels.size();
  double y;
  cut.push_front(0);
  cut.push_back(R_PosInf);
  
  if((cut[0] < x) & (x <= cut[1])){
    y = levels[0];
  }
  if (p > 1.5) {
    for (int i=1; i<p; i++) {
      //y[(cut[i] <= x) & (x <= cut[i + 1])] = levels[i];
      if((cut[i] < x) & (x <= cut[i+1])){
        y = levels[i];
      }
    }
  }
  if (logInd)
    y = log(y);
  return(y);
}

// [[Rcpp::export]]
double Hpwc_double(double x, NumericVector cuts, NumericVector levels, int logInd)
{
  NumericVector cut = sort_unique(cuts);
  int p = levels.size();
  //NumericVector y(x.size());
  double y;
  cut.push_front(0);
  cut.push_back(R_PosInf);
  NumericVector z;
  bool which = (cut[0] < x) & (x <= cut[0 + 1]);
  if (which) {
    //y[which] = x[which];
    y = x*levels[0];
  }
  double wt = levels[0] * cut[1];
  if (p > 1.5) {
    for (int i = 1; i<p; i++) {
      which = (cut[i] < x) & (x <= cut[i + 1]);
      if (which) {
        //NumericVector xwhich= x[which];
        double tmpx = wt + levels[i] * (x - cut[i]);
        y = tmpx;
      }
      wt = wt + levels[i] * (cut[i + 1] - cut[i]);
    }
  }
  
  if (logInd)
    y = log(y);
  return(y);
}



//[[Rcpp::export()]]
double ppwc_double(double q, NumericVector cuts, NumericVector levels, int lower, int logInd)
{
  double y;
  if (cuts[0]==0) {
    y = R::pexp(q, 1/levels[0], 0.0, 0.0);
  }else{
    
    y = Hpwc_double(q, cuts, levels, 0.0);
    if (logInd) {
      if (lower) {
        y = log(-(exp(-y)-1));
      }
      else {
        y = -y;
      }
    }
    else {
      if (lower) {
        y = -(exp(-y)-1);
      }
      else {
        y = exp(-y);
      }
    }
  }
  return(y);
}


//[[Rcpp::export()]]
Rcpp::NumericMatrix gaussleg(int n, double x1, double x2){
  
  double EPS = 3e-14;
  double m = 0.5*(n+1);
  double xm = 0.5 * (x2 + x1);
  double xl = 0.5 * (x2 - x1);
  
  Rcpp::NumericVector x=no_init(n);
  Rcpp::NumericVector w=no_init(n);
  double z;
  //double tmp;
  for (int i=0; i<m; i++) {
    
    double tmp = M_PI * ( 1.0*(i+1) - 0.25)/(n*1.0 + 0.5);
    z = cos(tmp);
    
    double tol = 9999;
    //double p1, p2, p3;
    double pp;
    //double z1;
    while (tol > EPS) {
      double p1 = 1.0;
      double p2 = 0.0;
      for (int j=0; j<n; j++) {
        double p3 = p2;
        p2 = p1;
        p1 = ((2 * j + 1) * z * p2 - j  * p3)/(j+1);
      }
      
      pp = n * (z * p1 - p2)/((z-1)*(1+z));
      double z1 = z;
      z = z1 - p1/pp;
      tol = abs(z - z1);
      //z -= p1/pp;
      //tol= abs(p1/pp);
    }
    
    int s=n -1 - i;
    x[i] = xm - xl * z;
    x[s] = xm + xl * z;
    w[i] = (2 * xl)/((1-z)*(1+z) * pp * pp);
    w[s] = w[i];
  }
  
  Rcpp::NumericMatrix res=no_init_matrix(n, 2);
  res(_,0) = x;
  res(_,1) =w;
  
  //double tmp1 = cos(tmp);
  //Rcout << "tmp is " << tmp << std::endl;
  //Rcout << "tmp1 is " << tmp1 << std::endl;
  return(res);
}

// [[Rcpp::export()]]
double survF01_gL_pwc01(double u, NumericVector cutsL, NumericVector levs01L){
  double res = ppwc_double(u, cutsL, levs01L, 0, 0);
  return res;
}

// [[Rcpp::export()]]
double int_ele_gau(double u1, double u2, double R, NumericVector cutsL, NumericVector levs01L, NumericVector lam12, NumericVector lam10){
  
  //double m=R-floor(L); 
  double res=0.0;
  
  Rcpp::NumericMatrix gau_a = gaussleg(20, u1, u2);
  
  for(int i=0; i<20; i++){
    res+=survF01_gL_pwc01(gau_a(i,0), cutsL, levs01L)*dbwei_mgf(R-gau_a(i,0), lam10[1], lam10[0], lam12[1], lam12[0])*hpwc_double(gau_a(i,0), cutsL, levs01L, 0)*gau_a(i,1);
  }
  
  //Rcout << "res is " << res << std::endl;
  
  return res; 
}

// [[Rcpp::export()]]
Rcpp::NumericVector int_vec(double L, double R, NumericVector cuts, NumericVector levs01, NumericVector lam12, NumericVector lam10){
  
  double m=R-floor(L); 
  NumericVector res=no_init(floor(m));
  
  LogicalVector ind=(cuts>L);
  LogicalVector inde=ind;
  inde.push_back(std::isfinite(L));
  
  NumericVector cutsL, levsL;
  if(is_true(any(ind))){
    cutsL=cuts[ind];
    levsL=levs01[inde];
    cutsL.push_front(L);
    levsL.push_front(0);
  }else{
    cutsL={L};;
    levsL={0,levs01[levs01.size()-1]};//{0, 0.075};//{0,levs01[levs01.size()-1]};
  }
  
  //Rcout << "cutsL is " << L << std::endl;
  //Rcout << "levsL is " << levsL << std::endl;
  
  double tmp=0.0;
  
  Rcpp::NumericMatrix gau_a = gaussleg(20, L, floor(L)+1);
  
  for(int i=0; i<20; i++){
    tmp+=survF01_gL_pwc01(gau_a(i,0), cutsL, levsL)*dbwei_mgf(R-gau_a(i,0), lam10[1], lam10[0], lam12[1], lam12[0])*hpwc_double(gau_a(i,0), cutsL, levsL, 0)*gau_a(i,1);
    //Rcout << "tmp1 is " << survF01_gL_pwc01(gau_a(i,0), cutsL, levsL) << std::endl;
    //Rcout << "tmp2 is " << dbwei_mgf(R-gau_a(i,0), lam10[1], lam10[0], lam12[1], lam12[0]) << std::endl;
    //Rcout << "gau_a(i,0) is " << gau_a(i,0) << std::endl;
    //Rcout << "tmp3 is " << hpwc_double(30.0034, 30, levsL, 0) << std::endl;
    //Rcout << "tmp4 is " << hpwc_double(gau_a(i,0), cutsL, levsL, 0) << std::endl;
  }
  
  res[0]=tmp;
  //res[0]=int_ele_gau(L, floor(L)+1, R, cutsL, levsL, lam12, lam10);
  //Rcout << "res[0] is " << tmp << std::endl;
  
  if(m>1){
    
    for(int k=1; k<m; k++){
      double r=floor(L)+k;
      
      double denom=survF01_gL_pwc01(r, cutsL, levsL);//ppwc_double(r, cutsL, levsL, 0, 0); //
      if(denom>0){
        res[k]=int_ele_gau(r, r+1, R, cutsL, levsL, lam12, lam10)/denom; 
        
        //Rcout << "int ele gau is " << int_ele_gau(r, r+1, R, cutsL, levsL, lam12, lam10) << std::endl;
        //Rcout << "test1 is " << 99 << std::endl;
        
      }else{
        res[k]=0.0;
        
        //Rcout << "test2 is " << 999 << std::endl;
      }
      
      //Rcout << "denom is " << denom << std::endl;
      //Rcout << "res[k] is " << res[k] << std::endl;
    }
    
  }
  
  return res; 
  
}


// [[Rcpp::export()]]
double p1_g_L(double L, int R, NumericVector cuts, NumericVector levs01, NumericVector lam12, NumericVector lam10){
  
  
  //Rcpp::List res0(K);
  NumericVector res(1);
  res[0]=1;
  double res1=0.0;
  
  if(R>L){
    
    int K=R-floor(L); 
    
    for(int k=0; k<K; k++){
      NumericVector intk=int_vec(L,floor(L)+k+1, cuts, levs01, lam12, lam10);
      double temp=1-sum(res*intk); 
      res.push_back(temp); 
    }
    
    res1=1-res[K];
    
  }
  
  if(res1<0 | res1>1){
    Rcout << "The value is " << res1 << std::endl;
    Rcout << "R is " << R << std::endl;
    
    if(std::abs(res1-1)<0.1){
      res1=1;
    }
    
    if(std::abs(res1)<1e-03){
      res1=0.0;
    }
    
  }
  
  //Rcout << "The value is " << res << std::endl;
  
  return res1;
}


// [[Rcpp::export()]]
Rcpp::NumericVector p1_g_L_pwc01_vec(NumericVector L, IntegerVector R, NumericVector cuts, NumericVector levs01, NumericVector lam12, NumericVector lam10){
  
  int n=R.size();
  NumericVector res=no_init(n);
  for(int i=0; i<n; i++){
    
    res[i]=p1_g_L(L[i], R[i], cuts, levs01, lam12, lam10);
    
  }
  
  return res;
}


// [[Rcpp::export()]]
double p1(int R, NumericVector cuts, NumericVector levs01, NumericVector lam12, NumericVector lam10, double meanlog, double sdlog, double Q){
  
  NumericMatrix gau_a = gaussleg(20, 0, R);
  NumericVector u = gau_a(_,0);
  NumericVector w = gau_a(_,1);
  
  double res=0.0;
  for(int i=0; i<20; i++){
    
    //LogicalVector ind=(cuts>u[i]);
    //LogicalVector inde=ind;
    //inde.push_back(std::isfinite(u[i]));
    
    //NumericVector cutsL, levsL;
    //if(is_true(any(ind))){
    //  cutsL=cuts[ind];
    //  levsL=levs01[inde];
    //  cutsL.push_front(u[i]);
    //  levsL.push_front(0);
    //}else{
    //  cutsL[0]=u[i];
    //  levsL={0,levs01[levs01.size()-1]};
    //}
    
    //double dens1=R::dlnorm(u[i],meanlog,sdlog,0);
    //double dens2=R::dlnorm(1/u[i],meanlog,sdlog,0);
    //double tmp=dgengam(u[i], meanlog, sdlog, Q);
    double temp1=p1_g_L(u[i], R, cuts, levs01, lam12, lam10)*w[i]*dgengam(u[i], meanlog, sdlog, Q);//R::dlnorm(u[i],meanlog,sdlog,0);
    //double temp2=p1_g_L_exp(1/u[i], R, lam01, lam12, lam10)*w[i]*pow(u[i],-2)*R::dlnorm(1/u[i],meanlog,sdlog,0);
    //res+=temp1+temp2;
    res+=temp1;
    
    //if(res<0 | res>1){
    //  Rcout << "temp1 is " << temp1 << std::endl;
    //  Rcout << "temp2 is " << temp2 << std::endl;
    //}
  }
  
  if(res<0 | res>1){
    Rcout << "The value is " << res << std::endl;
    Rcout << "R is " << R << std::endl;
    
    if(std::abs(res-1)<0.1){
      res=1;
    }
  }
  
  return res;
}


// [[Rcpp::export()]]
NumericVector p1_vec(NumericVector R, NumericVector cuts, NumericVector levs01, NumericVector lam12, NumericVector lam10, double meanlog, double sdlog, double Q){
  
  int n=R.size();
  NumericVector res=no_init(n);
  for(int i=0; i<n; i++){
    
    res[i]=p1(R[i], cuts, levs01, lam12, lam10, meanlog, sdlog, Q);
    
  }
  
  return res;
}
