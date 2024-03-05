#include <RcppArmadillo.h>
using namespace std;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List cov_cal_long( arma::colvec beta,
                         arma::colvec M,
                         arma::colvec rho01,
                         arma::colvec rho02,
                         arma::mat XS,
                         arma::colvec dati,
                         arma::colvec datj) {
  
  
  int n0 = M.size();
  int nij = rho01.size();
  int nb = beta.size();
  
  
  arma::mat sigma1(nb,nb);
  sigma1.fill(0.0);
  arma::mat sigma2(nb,nb);
  sigma2.fill(0.0);
  arma::cube f_store(n0,n0,nb);
  f_store.fill(0.0);
  arma::vec poly(nb);
  poly.fill(0.0);
  
  
  for (int i=0; i<nij; ++i){
    poly = XS.col(i);
    double temp = exp(beta.t()* poly)[0];
    sigma1=sigma1- (poly * poly.t()) * (rho02[i] + rho01[i]) *temp / pow((1+temp),2.0);
    arma::vec f_storeij =  poly /(1+temp) * ( rho01[i] - rho02[i] * temp);
    f_store.subcube(dati[i],datj[i],0,dati[i],datj[i],nb-1)=f_storeij;
  }
  
  for (int i=0; i<nij; ++i){
    int i2 = dati[i];
    int j2 = datj[i];
    for (int jj=0;jj<n0; ++jj){
      bool keep=(i2!=j2)&&(i2!=jj)&&(j2!=jj);
      if(keep){
        arma::vec vleft = f_store.subcube(i2,j2,0,i2,j2,nb-1)+f_store.subcube(j2,i2,0,j2,i2,nb-1);
        arma::vec vright = f_store.subcube(i2,jj,0,i2,jj,nb-1)+f_store.subcube(jj,i2,0,jj,i2,nb-1);
        sigma2+=vleft * vright.t();
        arma::vec vleft1 = f_store.subcube(j2,i2,0,j2,i2,nb-1)+f_store.subcube(i2,j2,0,i2,j2,nb-1);
        arma::vec vright1 = f_store.subcube(j2,jj,0,j2,jj,nb-1)+f_store.subcube(jj,j2,0,jj,j2,nb-1);
        sigma2+=vleft1 * vright1.t();
      }
    }           
  }
  return Rcpp::List::create(Named("sigma1")=sigma1,Named("sigma2")=sigma2);
}
