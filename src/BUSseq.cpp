//Bayesian Zero-Inflated Negative Binomial Mixture Model
//Compared with v2, put the updating of w later
//

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Print.h>
#include <Rmath.h>
#include <math.h>
//#include <gsl/gsl_sf_gamma.h>
//#include <fstream>

//#include <Rdefines.h>
//#define PI 3.14159265358
//#define TRUN 0.64

using namespace std;
/* get the list element named str, or return NULL */
extern "C"{
  //function statement
  SEXP getListElement(SEXP list, const char* str);
  
  //sampler
  void rdirichlet(double* xi, int k, double* rn);
  int rcate(double* prob, int k);
  
  //operate
  double vec_max(double* vec,int n);
  
  //update 
  void _update_zx(int _B, int *_Nb, int _N, int _G,
                  double *_gamma, double *_phi, double *_log_mu,
                  int *Y, 
                  int *Z, int *X);//to be updated
  
  void _update_gamma(int _B, int *_Nb, int _N, int _G,//dimension
                     double *_sigma_z,//prior
                     int *Z, int *X,// double *_omega, //latent variable
                     double *_gamma);//to be updated
  
  void _update_alpha(int _B, int *_Nb, int _N, int _G,//dimension
                     double* _mu_a, double _sigma_a,//prior
                     int *_w, double *_beta, double *_nu, double *_delta, double *_phi, //parameter
                     int *X, //latent variable
                     double *_alpha);
  
  void _update_l(int _G, int _K,
                 double *_ptau0, double _tau1,//hyperparameter
                 double *_beta,//parameter
                 int *L);//to be updated
  
  void _update_ptau0(int _G, int _K,
                     double _a_tau, double _b_tau, double _a_p, double _b_p,
                     int *L,//hyperparameter
                     double *_beta,//parameter
                     double *_ptau0);
  
  void _update_beta(int _B, int *_Nb, int _N, int _G, int _K,//dimension
                    double _tau1,//prior
                    double *_ptau0, int *L,//hyperparameter
                    int *_w, double *_alpha, double *_nu, double *_delta, double *_phi, //parameter
                    int *X, //latent variable
                    double *_beta);
  
  void _update_nu(int _B, int *_Nb, int _N, int _G,
                  double* _mu_c, double _sigma_c,//prior
                  int *_w, double *_alpha, double *_beta, double *_delta, double *_phi,//parameter
                  int *X, //latent variable
                  double *_nu);//to be updated
  
  void _update_delta(int _B, int *_Nb, int _N, int _G,
                     double* _mu_d, double _sigma_d,//prior
                     int *_w, double *_alpha, double *_beta, double *_nu, double *_phi,//parameter
                     int *X, //latent variable
                     double *_delta);//to be updated
  
  void _update_logmu(int _B, int *_Nb, int _N, int _G,
                     int *_w, double *_alpha, double *_beta, double *_nu, double *_delta, //parameter
                     double *_log_mu);//to be updated
  
  void _update_phi(int _B, int *_Nb, int _N, int _G,
                   double _kappa, double _tau,//prior
                   double *_log_mu,//parameter
                   int *X, //latent variable
                   double *_phi);//to be updated
  
  void _update_w(int _B, int *_Nb, int _N, int _G, int _K,
                 double *_pi,//prior
                 double *_gamma, double *_alpha, double *_beta, double *_nu, double *_delta, double *_phi,//parameter
                 int *X,//latent variable
                 //int up,
                 int *_w, int *_count_w);//to be updated
  
  void _update_pi(int _B, int *_Nb, int _K,
                  double _xi,//prior
                  int *_w,//parameter
                  double *_pi);//to be updated   
  
  //function content
  SEXP getListElement(SEXP list, const char* str){
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    
    for (int i = 0; i < length(list); i++)
      if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
        elmt = VECTOR_ELT(list, i);
        break;
      }
      return elmt;
  }
  
  //Sample from Dirichlet distribution
  void rdirichlet(double* xi, int k, double* rn) {
    double sum = 0.0;
    for (int l = 0; l < k; l++) {
      rn[l] = rgamma(xi[l], 1.0);
      sum = sum + rn[l];
    }
    for (int l = 0; l < k; l++) {
      rn[l] = rn[l] / sum;
    }
  }
  
  //Sample from categorical distribution
  int rcate(double* prob, int k) {
    int rn[k];
    int res = 0;
    
    rmultinom(1, prob, k, rn);
    while (rn[res] == 0.0) {
      res++;
    }
    return res;
  }
  
  //Max value in a vector
  double vec_max(double* vec,int n){
    double res = vec[0];
    for (int i = 1; i < n; i++) {
      if (res < vec[i]) {
        res = vec[i];
      }
    }
    return res;
  }
  
  //update z_{big} and x_{big} 
  void _update_zx(int _B, int *_Nb, int _N, int _G,
                  double *_gamma, double *_phi, double *_log_mu,
                  int *Y, int *Z, int *X){
    
    int ind, ind_nu;
    int ind_n = 0; //index the row of (b,i)	
    
    
    for (int b = 0; b < _B; b ++){
      for(int i = 0; i < _Nb[b]; i ++){
        for(int j = 0; j < _G; j ++){
          ind = ind_n + j * _N;
          ind_nu = b + j * _B;
          if (Y[ind] == 0) {
            //Rprintf("y %d %d %d = 0\n", b , i , g);
            //update z_{big}
            if (X[ind] == 0) {
              double log_rat = -_gamma[b];
              Z[ind] = rbinom(1, 1 / (1 + exp(log_rat)));
            }
            else {
              Z[ind] = 0;
            }
            
            //Rprintf("z %d %d %d = %d\n", b, i, j, Dropout[index]);
            //update x_{big}
            if (Z[ind] == 0) {
              //MH sampling
              double pb = _phi[ind_nu] / (_phi[ind_nu] + exp(_log_mu[ind]));
              double temp_x = rnbinom(_phi[ind_nu], pb);
              double u = runif(0, 1);
              //make the exponential be compuational
              double temp_max = 0.0;
              if(temp_max<_gamma[b] + _gamma[b + _B] * temp_x){
                temp_max = _gamma[b] + _gamma[b + _B] * temp_x;
              }
              if(temp_max<_gamma[b] + _gamma[b + _B] * X[ind]){
                temp_max = _gamma[b] + _gamma[b + _B] * X[ind];
              }
              
              double ratio = (exp(-temp_max) + exp(_gamma[b] + _gamma[b + _B] * X[ind] - temp_max)) / (exp(-temp_max) + exp(_gamma[b] + _gamma[b + _B] * temp_x - temp_max));
              if (u < ratio) {
                X[ind] = int(temp_x);
              }
            }
            else {
              X[ind] = 0;
            }
            //Rprintf("x %d %d %d= %d\n", b , i , j, X[ind]);
          }
        }
        //Rprintf("Finish sampling the %d-th cell\n", ind_n);
        ind_n = ind_n + 1;
      }
    }
  }
  
  //update gamma_b
  void _update_gamma(int _B, int *_Nb, int _N, int _G,//dimension
                     double *_sigma_z,//prior
                     int *Z, int *X,// double *_omega, //latent variable
                     double *_gamma){
    
    int ind;
    int ind_n = 0; //index the row of (b,i)
    double gamma_iter[2];
    double logr;
    //update gamma_0
    for (int b = 0; b < _B; b++) {
      
      gamma_iter[0] = rnorm(_gamma[b], 0.1);
      logr = 0.0;
      
      //prior
      logr = logr - pow(gamma_iter[0], 2.0) / 2 / _sigma_z[0] + pow(_gamma[b], 2.0) / 2 / _sigma_z[0];
      
      for(int i = 0; i < _Nb[b]; i++ ){
        for(int j = 0; j < _G; j++){
          ind = ind_n + j * _N;
          //numerator
          double temp_iter = gamma_iter[0] + X[ind] * _gamma[b + _B];//prevent temp_iter is extremely large
          if(temp_iter > 0){
            logr = logr + gamma_iter[0] * Z[ind] - temp_iter  - log(1 + exp(-temp_iter) );
          }else{
            logr = logr + gamma_iter[0] * Z[ind] - log(1 + exp(temp_iter));
          }
          
          //denomerator
          double temp_pres = _gamma[b] + X[ind] * _gamma[b + _B];
          if(temp_pres>0){
            logr = logr - _gamma[b] * Z[ind] + temp_pres  + log(1 + exp(-temp_iter));
          }else{
            logr = logr - _gamma[b] * Z[ind]  + log(1 + exp(temp_iter) );
          }
        }
        ind_n = ind_n + 1;
      }
      if (logr > log(runif(0, 1))) {
        _gamma[b] = gamma_iter[0];
      }
    }
    
    //update gamma_1
    ind_n = 0;
    for (int b = 0; b < _B; b++) {
      
      //proposal
      gamma_iter[1] = rgamma(10 * _gamma[b + _B], 0.1);
      logr = 0.0;
      
      //prior
      logr = logr - pow(gamma_iter[1], 2.0) / 2 / _sigma_z[1] + pow(_gamma[b + _B], 2.0) / 2 / _sigma_z[1];
      
      //proposal
      //numerator
      logr = logr - lgamma(gamma_iter[1]) + (gamma_iter[1] - 1) * log(_gamma[b + _B]) - _gamma[b + _B];
      //denomerator
      logr = logr + lgamma(_gamma[b + _B]) + (_gamma[b + _B] - 1) * log(gamma_iter[1]) - gamma_iter[1];
      
      for (int i = 0; i < _Nb[b]; i++) {
        for (int j = 0; j < _G; j++) {
          ind = ind_n + j * _N;
          //numerator
          double temp_iter = _gamma[b] + X[ind] * gamma_iter[1];//prevent temp_iter is extremely large
          if(temp_iter > 0){
            logr = logr + X[ind] * gamma_iter[1] * Z[ind] - temp_iter  - log(1 + exp(-temp_iter));
          }else{
            logr = logr + X[ind] * gamma_iter[1] * Z[ind] - log(1 + exp(temp_iter));
          }
          
          //denomerator
          double temp_pres = _gamma[b] + X[ind] * _gamma[b + _B];
          if(temp_pres>0){
            logr = logr - X[ind] * _gamma[b + _B] * Z[ind] + temp_pres + log(1 + exp(-temp_pres));
          }else{
            logr = logr - X[ind] * _gamma[b + _B] * Z[ind] + log(1 + exp(temp_pres));
          }
        }
        ind_n = ind_n + 1;
      }
      if (logr > log(runif(0, 1))) {
        _gamma[b + _B] = gamma_iter[1];
      }
    }
  }
  
  //update alpha_g
  void _update_alpha(int _B, int *_Nb, int _N, int _G,//dimension
                     double* _mu_a, double _sigma_a,//prior
                     int *_w, double *_beta, double *_nu, double *_delta, double *_phi, //parameter
                     int *X, //latent variable
                     double *_alpha){
    
    int ind, ind_beta, ind_nu;
    int ind_n; //index the row of (b,i)
    
    for (int j = 0; j < _G; j++) {
      //proposal
      double alpha_iter = rnorm(_alpha[j], 0.1);
      double logr = 0.0;
      
      //prior
      logr = logr - pow(alpha_iter - _mu_a[j], 2.0) / 2 / pow(_sigma_a, 2.0) + pow(_alpha[j] - _mu_a[j], 2.0) / 2 / pow(_sigma_a, 2.0);
      
      ind_n = 0;
      for (int b = 0; b < _B; b++) {
        for (int i = 0; i < _Nb[b]; i++) {
          ind = ind_n + j * _N;
          ind_beta = j + _w[ind_n] * _G;//bgw_i
          ind_nu = b + j * _B;
          
          //numerator
          logr = logr + alpha_iter * X[ind] - (_phi[ind_nu] + X[ind]) * log(_phi[ind_nu] + exp(alpha_iter + _beta[ind_beta] + _nu[ind_nu] + _delta[ind_n]));
          //denomerator
          logr = logr - _alpha[j] * X[ind] + (_phi[ind_nu] + X[ind]) * log(_phi[ind_nu] + exp(_alpha[j] + _beta[ind_beta] + _nu[ind_nu] + _delta[ind_n]));
          
          ind_n = ind_n + 1;
        }
      }
      
      if (logr > log(runif(0, 1))) {
        _alpha[j] = alpha_iter;
      }
    }
  }
  
  //update L_gk
  void _update_l(int _G, int _K,
                 double *_ptau0, double _tau1,//hyperparameter
                 double *_beta,//parameter
                 int *L){
    int ind_beta;
    for (int j = 0; j < _G; j++) {
      for (int k = 1; k < _K; k++) {
        ind_beta = j + k * _G;
        double log_rat = 0.0; //the ratio between L_gk = 1 and L_gk = 0
        log_rat = log_rat + log(_ptau0[0]) - log(1 - _ptau0[0]);
        log_rat = log_rat - log(_tau1)/2.0 + log(_ptau0[1])/2.0;
        log_rat = log_rat - pow(_beta[ind_beta], 2) / 2.0 / _tau1 + pow(_beta[ind_beta], 2) / 2.0 / _ptau0[1];
        //Rprintf("rate = %f\n", log_rat);
        L[ind_beta] = rbinom(1, 1 / (1 + exp(-log_rat)));
        //Rprintf("L %d %d= %d\n", j, k, Indicator[index_beta]);
      }
    }
  }
  
  //update p and tau0
  void _update_ptau0(int _G, int _K,
                     double _a_tau, double _b_tau, double _a_p, double _b_p,
                     int *L,//hyperparameter
                     double *_beta,//parameter
                     double *_ptau0){
    //update p and tau_0
    int ind_beta;
    int sum_L = 0;
    double sum_beta = 0.0;
    for (int j = 0; j < _G; j++) {
      for (int k = 1; k < _K; k++) {
        ind_beta = j + k * _G;
        sum_L = sum_L + L[ind_beta];
        if (L[ind_beta] == 0) {
          sum_beta = sum_beta + pow(_beta[ind_beta], 2);
        }
      }
    }
    double ap_post = _a_p + sum_L;
    double bp_post = _b_p + _G * (_K - 1) - sum_L;
    _ptau0[0] = rbeta(ap_post, bp_post);
    double atau_post = _a_tau + (_G * (_K - 1) - sum_L) / 2.0;
    double btau_post = _b_tau + sum_beta / 2.0;
    _ptau0[1] = 1.0 / rgamma(atau_post, 1 / btau_post);
  }
  
  //update beta_gk
  void _update_beta(int _B, int *_Nb, int _N, int _G, int _K,//dimension
                    double _tau1,//prior
                    double* _ptau0, int *L,//hyperparameter
                    int *_w, double *_alpha, double *_nu, double *_delta, double *_phi, //parameter
                    int *X, //latent variable
                    double *_beta){
    
    int ind, ind_beta, ind_nu;
    int ind_n; //index the row of (b,i)
    double beta_iter[_K];
    double logr[_K];
    for (int j = 0; j < _G; j++) {
      for (int k = 1; k < _K; k++) {
        //proposal
        ind_beta = j + k * _G;
        beta_iter[k] = rnorm(_beta[ind_beta], 0.1);
        logr[k] = 0.0;
        
        //prior
        if (L[ind_beta] == 1) {
          logr[k] = logr[k] - pow(beta_iter[k], 2.0) / 2 / _tau1 + pow(_beta[ind_beta], 2.0) / 2 / _tau1;
        }
        else {
          logr[k] = logr[k] - pow(beta_iter[k], 2.0) / 2 / _ptau0[1] + pow(_beta[ind_beta], 2.0) / 2 / _ptau0[1];
        }
      }
      
      ind_n = 0;
      for (int b = 0 ; b < _B; b++){
        for (int i = 0; i < _Nb[b]; i++) {
          int k_temp = _w[ind_n];
          ind = j * _N + ind_n;
          ind_beta = j + k_temp * _G;//bgw_bi
          ind_nu = b + j * _B;
          
          //numerator
          logr[k_temp] = logr[k_temp] + beta_iter[k_temp] * X[ind] - (_phi[ind_nu] + X[ind]) * log(_phi[ind_nu] + exp(_alpha[j] + beta_iter[k_temp] + _nu[ind_nu] + _delta[ind_n]));
          //denomerator
          logr[k_temp] = logr[k_temp] - _beta[ind_beta] * X[ind] + (_phi[ind_nu] + X[ind]) * log(_phi[ind_nu] + exp(_alpha[j] + _beta[ind_beta] + _nu[ind_nu] + _delta[ind_n]));
          
          ind_n = ind_n + 1;
        }	
      }
      
      for (int k = 1; k < _K; k++) {
        //proposal
        ind_beta = j + k * _G;
        if (logr[k] > log(runif(0, 1))) {
          _beta[ind_beta] = beta_iter[k];
        }
      }
      
    }
  }
  
  //update nu_bg
  void _update_nu(int _B, int *_Nb, int _N, int _G,
                  double* _mu_c, double _sigma_c,//prior
                  int *_w, double *_alpha, double *_beta, double *_delta, double *_phi,//parameter
                  int *X, //latent variable
                  double *_nu){//to be updated
    
    int ind, ind_beta, ind_nu, ind_n;
    for(int j = 0; j < _G; j++){
      ind_n = _Nb[0];
      for(int b = 1; b < _B; b ++){
        ind_nu = b + j * _B;
        
        //proposal
        double nu_iter = rnorm(_nu[ind_nu], 0.1);
        double logr = 0.0;
        
        //prior
        logr = logr - pow(nu_iter - _mu_c[b], 2.0) / 2 / pow(_sigma_c, 2.0) + pow(_nu[ind_nu] - _mu_c[b], 2.0) / 2 / pow(_sigma_c, 2.0);
        
        for(int i = 0; i < _Nb[b]; i ++){
          
          ind_beta = j + _w[ind_n] * _G;
          ind = j * _N + ind_n;
          //numerator
          logr = logr + nu_iter * X[ind] - (_phi[ind_nu] + X[ind]) * log(_phi[ind_nu] + exp(_alpha[j] + _beta[ind_beta] + nu_iter + _delta[ind_n]));
          //denomerator
          logr = logr - _nu[ind_nu] * X[ind] + (_phi[ind_nu] + X[ind]) * log(_phi[ind_nu] + exp(_alpha[j] + _beta[ind_beta] + _nu[ind_nu] + _delta[ind_n]));
          
          ind_n = ind_n + 1;
        }
        
        if (logr > log(runif(0, 1))) {
          _nu[ind_nu] = nu_iter;
        }			
      }
    }
  }
  
  //update delta_bi
  void _update_delta(int _B, int *_Nb, int _N, int _G,
                     double* _mu_d, double _sigma_d,//prior
                     int *_w, double *_alpha, double *_beta, double *_nu, double *_phi,//parameter
                     int *X, //latent variable
                     double *_delta){//to be updated
    
    int ind, ind_beta, ind_nu, ind_n;
    ind_n = 0;
    for(int b = 0; b < _B; b++){
      for(int i = 0; i < _Nb[b]; i ++){
        
        if(i > 0){//let ind_n be the cell index
          double delta_iter = rnorm(_delta[ind_n], 0.1);
          double logr = 0.0;
          
          //prior
          logr = logr - pow(delta_iter - _mu_d[ind_n], 2.0) / 2 / pow(_sigma_d, 2.0) + pow(_delta[ind_n] - _mu_d[ind_n], 2.0) / 2 / pow(_sigma_d, 2.0);
          
          for(int j = 0; j < _G;j ++){
            ind_nu = b + j * _B;
            ind_beta = j + _w[ind_n] * _G;
            ind = j * _N + ind_n;
            
            //numerator
            logr = logr + delta_iter * X[ind] - (_phi[ind_nu] + X[ind]) * log(_phi[ind_nu] + exp(_alpha[j] + _beta[ind_beta] + _nu[ind_nu] + delta_iter));
            //denomerator
            logr = logr - _delta[ind_n] * X[ind] + (_phi[ind_nu] + X[ind]) * log(_phi[ind_nu] + exp(_alpha[j] + _beta[ind_beta] + _nu[ind_nu] + _delta[ind_n]));
            
          }
          if (logr > log(runif(0, 1))) {
            _delta[ind_n] = delta_iter;
          }	
          
        }
        ind_n = ind_n + 1;
      }
    }
  }
  
  //update log_mu
  void _update_logmu(int _B, int *_Nb, int _N, int _G,
                     int *_w, double *_alpha, double *_beta, double *_nu, double *_delta,//parameter
                     double *_log_mu){//to be updated
    
    int ind, ind_n, ind_beta, ind_nu;
    ind_n = 0;
    //auxiliry
    for(int b = 0; b < _B; b++) {
      for (int i = 0; i < _Nb[b]; i++) {
        for (int j = 0; j < _G; j++) {
          ind = j * _N + ind_n;
          ind_beta = j + _w[ind_n] * _G;
          ind_nu = b + j * _B;
          _log_mu[ind] = _alpha[j] + _beta[ind_beta] + _nu[ind_nu] + _delta[ind_n];
        }
        ind_n = ind_n + 1;
      }
    }
  }
  
  //update phi_bg
  void _update_phi(int _B, int *_Nb, int _N, int _G,
                   double _kappa, double _tau,//prior
                   double *_log_mu,//parameter
                   int *X, //latent variable
                   double *_phi){//to be updated
    
    int ind, ind_nu, ind_n;
    for (int j = 0; j < _G; j++)	{
      ind_n = 0;
      for (int b = 0; b < _B; b++) {
        ind_nu = b + j * _B;
        double phi_iter = rgamma(_phi[ind_nu], 1);
        double logr = 0.0;
        
        for (int i = 0; i < _Nb[b]; i++) {
          ind = j * _N + ind_n;
          
          //numerator
          logr = logr + lgamma(phi_iter + X[ind]) - lgamma(phi_iter);
          logr = logr + phi_iter * log(phi_iter) - (phi_iter + X[ind]) * log(phi_iter + exp(_log_mu[ind]));
          //denomerator
          logr = logr - lgamma(_phi[ind_nu] + X[ind]) + lgamma(_phi[ind_nu]);
          logr = logr - _phi[ind_nu] * log(_phi[ind_nu]) + (_phi[ind_nu] + X[ind]) * log(_phi[ind_nu] + exp(_log_mu[ind]));
          ind_n = ind_n + 1;
        }
        
        //prior
        logr = logr + (_kappa - 1) * log(phi_iter) - _tau * phi_iter;
        logr = logr - (_kappa - 1) * log(_phi[ind_nu]) + _tau * _phi[ind_nu];
        //proposal
        logr = logr + (phi_iter - 1.0) * log(_phi[ind_nu]) + phi_iter - lgamma(phi_iter);
        logr = logr - (_phi[ind_nu] - 1.0) * log(phi_iter) - _phi[ind_nu] + lgamma(_phi[ind_nu]);
        
        if (logr > log(runif(0, 1))) {
          _phi[ind_nu] = phi_iter;
        }
      }
    }
  }
  
  //update w_bi
  void _update_w(int _B, int *_Nb, int _N, int _G, int _K,
                 double *_pi,//prior
                 double *_gamma, double *_alpha, double *_beta, double *_nu, double *_delta, double *_phi,//parameter
                 int *X,//latent variable
                 //int up,
                 int *_w, int *_count_w){//to be updated
    
    double log_proposal, log_current;
    int w_proposal, w_current;
    double proposal_pi[_K];
    for (int k = 0; k < _K; k++) {
      proposal_pi[k]=1.0/_K;
    }
    
    int ind, ind_nu, ind_beta;
    int ind_n = 0;
    double temp_logmu;
    double prob_y0;
    
    //Rprintf("Start updating w_bi:\n");
    for (int b = 0; b < _B; b++) {
      
      for (int i = 0; i < _Nb[b]; i++) {//update w_bi
        
        //get a proposal
        w_proposal = rcate(proposal_pi, _K);
        w_current = _w[ind_n];
        log_proposal = log(_pi[w_proposal + b * _K]);
        log_current = log(_pi[w_current + b * _K]);
        //calculate the posterior ratio in log scale
        for (int j = 0; j < _G; j++) {
          ind = j * _N + ind_n;
          ind_nu = b + j * _B;
          ind_beta = j + w_proposal * _G;
          temp_logmu = _alpha[j] + _beta[ind_beta] + _nu[ind_nu] + _delta[ind_n];
          log_proposal = log_proposal + _beta[ind_beta] * X[ind];
          log_proposal = log_proposal - (_phi[ind_nu] + X[ind]) * log(_phi[ind_nu] + exp(temp_logmu));
          
          ind_beta = j + w_current * _G;
          temp_logmu = _alpha[j] + _beta[ind_beta] + _nu[ind_nu] + _delta[ind_n];
          log_current = log_current + _beta[ind_beta] * X[ind]; 
          log_current = log_current - (_phi[ind_nu] + X[ind]) * log(_phi[ind_nu] + exp(temp_logmu));

        }
        double log_ratio;
        log_ratio = log_proposal - log_current;
        if (log_ratio > log(runif(0, 1))) {
          _w[ind_n] = w_proposal;
        }
        
        _count_w[_w[ind_n]] = _count_w[_w[ind_n]] + 1;
        ind_n = ind_n + 1;
      }
    }
  }

  //update pi_bk
  void _update_pi(int _B, int *_Nb, int _K,
                  double _xi,//prior
                  int *_w,//parameter
                  double *_pi){
    
    int ind_n = 0;
    for(int b = 0; b < _B; b++){
      double post_xi[_K];
      for (int k = 0; k < _K; k++) {
        post_xi[k] = _xi;
      }
      for (int i = 0; i < _Nb[b]; i++) {
        post_xi[_w[ind_n]] = post_xi[_w[ind_n]] + 1.0;
        ind_n = ind_n + 1;
      }
      rdirichlet(post_xi, _K, _pi + b * _K);
    }
  }
  
  /* Main MCMC function*/
  
  SEXP BZINBBUS(SEXP args) {
    
    int nProtect = 0;
    //Rprintf("Input data.\n");
    
    // load list from R into Cpp 
    // load observed data
    int *read = INTEGER(getListElement(args, "Reads"));
    //[Y_11,Y_21,...,Y_N1,Y_12,Y_22,...,Y_NG]
    int B = INTEGER(getListElement(args, "BNum"))[0];
    int K = INTEGER(getListElement(args, "K"))[0];
    int* Nb = INTEGER(getListElement(args, "CNum"));
    int G = INTEGER(getListElement(args, "GNum"))[0];

    int N = 0;
    for (int b = 0; b < B; b++) {
      N = N + Nb[b];
      //Rprintf("The number of Cells in %d-th batch: %d\n", b + 1, Nb[b]);
    }
    //Rprintf("The total number of Cells: %d\n",N);
    //Rprintf("The number of Genes: %d\n", G);
    
    //load the indicator to control whether or not to update p and tau0
    int IND_UPDATE_PTAU0 = INTEGER(getListElement(args,"ind_update_ptau"))[0];
    
    // Intial value of paramter
    double* alpha = REAL(getListElement(args,"alpha.init"));
    //the form of alpha [alpha_1,alpha_2,...,alpha_G] G
    double* beta = REAL(getListElement(args, "beta.init"));
    //the form of beta [beta_11,beta_21,...,beta_G1,beta_12,beta_22,...,beta_GK] G * K
    double* nu = REAL(getListElement(args, "nu.init"));
    //the form of nu [nu_11,nu_21,...,beta_B1,beta_12,beta_22,...,beta_BG] B * G
    double* delta = REAL(getListElement(args, "delta.init"));
    //the form of delta [delta_11,delta_12,...,delta_1n_1,delta_21,...,delta_Bn_B] N
    double* gamma = REAL(getListElement(args, "gamma.init"));
    //the form of gamma [gamma_11,gamma_21,...,gamma_B1,gamma_12,gamma_22,...,gamma_B2] B * 2
    double* phi = REAL(getListElement(args, "phi.init"));
    //the form of phi [phi_11,phi_21,...,phi_B1,phi_12,phi_22,...,phi_BG] B * G
    double* pi = REAL(getListElement(args, "pi.init"));
    //the form of pi [pi_11,pi_21,...,pi_K1,pi_12,pi_22,...,pi_KB] K * B
    double p = REAL(getListElement(args, "p.init"))[0];
    double tau0 = REAL(getListElement(args, "tau0.init"))[0];
    double ptau0[2];
    ptau0[0] = p;
    ptau0[1] = tau0;
    
    //load initial value of w, z and x
    int* w = INTEGER(getListElement(args, "w.init"));
    int* Dropout = INTEGER(getListElement(args, "z.init"));
    int* Trueread = INTEGER(getListElement(args, "x.init"));
    int* Indicator = INTEGER(getListElement(args, "l.init"));
    //double Like = 0.0;
    
    //load prior
    double xi = REAL(getListElement(args, "xi"))[0];
    double* mu_a = REAL(getListElement(args, "mu.a"));
    double sigma_a = REAL(getListElement(args, "sigma.a"))[0];
    //double sigma_b = REAL(getListElement(args, "sigma.b"))[0];
    double* mu_c = REAL(getListElement(args, "mu.c"));
    double sigma_c = REAL(getListElement(args, "sigma.c"))[0];
    
    SEXP Mu_d = PROTECT(allocVector(REALSXP, N));
    nProtect++;
    double* mu_d = REAL(Mu_d);
    for (int q = 0; q < N; q++) {
      mu_d[q] = delta[q];
    }
    
    double sigma_d = REAL(getListElement(args, "sigma.d"))[0];
    double* sigma_z = REAL(getListElement(args, "sigma.z"));
    double kappa = REAL(getListElement(args, "kappa"))[0];
    double tau = REAL(getListElement(args, "tau"))[0];
    double a_tau = REAL(getListElement(args, "a.tau"))[0];
    double b_tau = REAL(getListElement(args, "b.tau"))[0];
    double a_p = REAL(getListElement(args, "a.p"))[0];
    double b_p = REAL(getListElement(args, "b.p"))[0];
    double tau1 = REAL(getListElement(args, "tau1"))[0];
    
    //Rprintf("Input has finished.\n");
    
    // MCMC parameter
    int iter_num;
    iter_num = INTEGER(getListElement(args, "n.out"))[0];
    
    SEXP par_post;
    //containing the posterior sampling of all parameters:
    // G*iter_num for alpha
    // G*K*iter_num for beta
    // G*iter_num for phi
    // K*iter_num for pi
    // (G + G*K +  G + K)*iter_num
    
    //record the posterior sampling of alpha
    SEXP alpha_post = PROTECT(allocVector(REALSXP, iter_num * G));
    nProtect++;
    
    //record the posterior sampling of beta
    SEXP beta_post = PROTECT(allocVector(REALSXP, iter_num * G * K));
    nProtect++;
    
    //record the posterior sampling of nu
    SEXP nu_post = PROTECT(allocVector(REALSXP, iter_num * B * G));
    nProtect++;
    
    //record the posterior sampling of delta
    SEXP delta_post = PROTECT(allocVector(REALSXP, iter_num * N));
    nProtect++;
    
    //record the posterior saampling of gamma
    SEXP gamma_post = PROTECT(allocVector(REALSXP, iter_num * 2 * B));
    nProtect++;
    
    //record the posterior sampling of phi
    SEXP phi_post = PROTECT(allocVector(REALSXP, iter_num * B * G));
    nProtect++;
    
    //record the posterior sampling of pi
    SEXP pi_post = PROTECT(allocVector(REALSXP, iter_num * K * B));
    nProtect++;
    
    //record the posterior sampling of w
    SEXP w_post = PROTECT(allocVector(INTSXP, iter_num * N));
    nProtect++;
    
    //record the posterior sampling of p
    SEXP p_post = PROTECT(allocVector(REALSXP, iter_num));
    nProtect++;
    
    //record the posterior sampling of tau0
    SEXP tau0_post = PROTECT(allocVector(REALSXP, iter_num));
    nProtect++;
    
    //record the posterior sampling of l
    SEXP ind_post = PROTECT(allocVector(INTSXP, iter_num * G * K));
    nProtect++;
    
    PROTECT(par_post = allocVector(VECSXP, 11));
    SET_VECTOR_ELT(par_post, 0 , alpha_post);
    SET_VECTOR_ELT(par_post, 1, beta_post);
    SET_VECTOR_ELT(par_post, 2, nu_post);
    SET_VECTOR_ELT(par_post, 3, delta_post);
    SET_VECTOR_ELT(par_post, 4, gamma_post);
    SET_VECTOR_ELT(par_post, 5, phi_post);
    SET_VECTOR_ELT(par_post, 6, pi_post);
    SET_VECTOR_ELT(par_post, 7, w_post);
    SET_VECTOR_ELT(par_post, 8, p_post);
    SET_VECTOR_ELT(par_post, 9, tau0_post);
    SET_VECTOR_ELT(par_post, 10, ind_post);
    nProtect++;
    
    
    //Rprintf("MCMC starts\n");
    GetRNGstate();
    
    SEXP logmu = PROTECT(allocVector(REALSXP, N * G));
    nProtect++;
    
    double* log_mu = REAL(logmu);
    _update_logmu(B, Nb, N, G,
                  w, alpha, beta, nu, delta,//parameter
                  log_mu);//to be updated
    
    int count_w[K];
    
    for (int iter = 0; iter < iter_num; iter++) {
      
      //update z_{big} and x_{big}
      _update_zx(B, Nb, N, G,//dimension
                 gamma, phi, log_mu,//parameter
                 read, //observed data
                 Dropout, Trueread);//latent variable to be updated
      
      //Rprintf("Finish %d-th sampling of z and x.\n", iter + 1);
      
      //update gamma_{b}
      _update_gamma(B, Nb, N, G,//dimension
                    sigma_z,//prior
                    Dropout, Trueread,// omega, //latent variable
                    gamma);//to be updated
      
      //Rprintf("Finish %d-th Gibbs sampling of gamma.\n", iter + 1);
      
      //update alpha_g by MH
      _update_alpha(B, Nb, N, G,//dimension
                    mu_a, sigma_a,//prior
                    w, beta, nu, delta, phi, //parameter
                    Trueread, //latent variable
                    alpha);//to be updated
      
      //Rprintf("Finish %d-th Gibbs sampling of alpha.\n",iter + 1);
      
      //update L_gk
      _update_l(G, K,
                ptau0, tau1,//hyperparameter
                beta,//parameter
                Indicator);
      
      if (IND_UPDATE_PTAU0 == 1){
        //update p and tau0
        _update_ptau0(G, K,
                      a_tau, b_tau, a_p, b_p,//hyperparameter
                      Indicator,
                      beta,//parameter
                      ptau0);
        
        //Rprintf("Finish %d-th Gibbs sampling of p and tau0.\n", iter + 1);
      }
      
      
      //update beta_bg
      _update_beta(B, Nb, N, G, K,//dimension 
                   tau1,//prior
                   ptau0, Indicator,//hyperparameter
                   w, alpha, nu, delta, phi, //parameter 	
                   Trueread,//latent variable
                   beta);
      
      //Rprintf("Finish %d-th Gibbs sampling of beta.\n",iter + 1);
      
      //update nu_bg
      _update_nu(B, Nb, N, G,
                 mu_c, sigma_c,//prior
                 w, alpha, beta, delta, phi,//parameter
                 Trueread, //latent variable
                 nu);//to be updated
      
      //Rprintf("Finish %d-th Gibbs sampling of nu.\n",iter + 1);
      
      //update nu_bg
      _update_delta(B, Nb, N, G,
                    mu_d, sigma_d,//prior
                    w, alpha, beta, nu, phi,//parameter
                    Trueread, //latent variable
                    delta);//to be updated
      
      //Rprintf("Finish %d-th Gibbs sampling of delta.\n",iter + 1);
      
      //update log_mu
      _update_logmu(B, Nb, G, K,
                    w, alpha, beta, nu, delta,//parameter
                    log_mu);//to be updated
      
      //Rprintf("Finish %d-th logmu.\n", iter + 1); 
      
      //update phi_bg
      _update_phi(B, Nb, N, G,
                  kappa, tau,//prior
                  log_mu,//parameter
                  Trueread, //latent variable
                  phi);//to be updated
      
      //Rprintf("Finish %d-th MH of phi.\n",iter + 1);
      
      //sample subtype effect for the empty subtype
      for (int k = 0; k < K; k++) {
        count_w[k] = 0;
      }

      //update w_bi
      for (int k = 0; k < K; k++) {
        count_w[k] = 0;
      }
      
      _update_w(B, Nb, N, G, K,
                pi,//prior
                gamma, alpha, beta, nu, delta, phi,//parameter
                Trueread,//latent variable
                //1000,
                w, count_w);
      
      //check empty subtype
      for(int k = 0; k < K; k++){
        //Rprintf("There are %d cells in the subtype %d\n", count_w[k],k + 1);
      }
      //Rprintf("Finish %d-th MH sampling of w.\n", iter + 1);
      
      //update log_mu
      _update_logmu(B, Nb, G, K,
                    w, alpha, beta, nu, delta,//parameter
                    log_mu);//to be updated
      
      //update pi_bk
      _update_pi(B, Nb, K,
                 xi,//prior
                 w,//parameter
                 pi);//to be updated
      
      //Rprintf("Finish %d-th Gibbs sampling of pi.\n",iter + 1); 
      
      //record the posterior sampling
      //int index_record;
      //record alpha G
      for(int j = 0; j < G;j++) {
        REAL(alpha_post)[iter * G + j] = *(alpha + j);
      }
      
      //record beta G * K
      for (int q = 0; q < G * K; q++) {
        REAL(beta_post)[iter * K * G + q] = *(beta + q);
      }
      
      //record nu B * G
      for (int q = 0; q < G * B; q++) {
        REAL(nu_post)[iter * B * G + q] = *(nu + q);
      }
      
      //record nu B * G
      for (int q = 0; q < N; q++) {
        REAL(delta_post)[iter * N + q] = *(delta + q);
      }
      
      //record gamma B * 2
      for (int q = 0; q < B * 2; q++) {
        REAL(gamma_post)[iter * 2 * B + q] = *(gamma + q);
      }
      
      //record phi B * G
      for (int q = 0; q < B * G; q++) {
        REAL(phi_post)[iter * B * G + q] = *(phi + q);
      }
      
      //record pi B * K
      for (int q = 0; q < B * K; q++) {
        REAL(pi_post)[iter * B * K + q] = *(pi + q);
      }
      
      //record w N
      for (int i = 0; i < N; i++) {
        INTEGER(w_post)[iter * N + i] = *(w + i);
      }	
      
      //record p 1
      REAL(p_post)[iter] = ptau0[0];
      
      //record tau0 1
      REAL(tau0_post)[iter] = ptau0[1];
      
      //record L_gk
      for (int q = 0; q < G * K; q++) {
        INTEGER(ind_post)[iter * G * K + q] = *(Indicator + q);
      }
      
      
    }
    
    PutRNGstate();
    
    UNPROTECT(nProtect);
    return par_post;
    }
  
  SEXP Cal_LogLike(SEXP args) {
    int nProtect = 0;
    //Rprintf("Input data.\n");
    
    // load list from R into Cpp 
    // load observed data
    int *Reads = INTEGER(getListElement(args, "Reads"));
    //[Y_11,Y_21,...,Y_N1,Y_12,Y_22,...,Y_NG]
    int B = INTEGER(getListElement(args, "BNum"))[0];
    int K = INTEGER(getListElement(args, "K"))[0];
    int* Nb = INTEGER(getListElement(args, "CNum"));
    int G = INTEGER(getListElement(args, "GNum"))[0];
    int N = 0;
    for (int b = 0; b < B; b++) {
      N = N + Nb[b];
      //Rprintf("The number of Cells in %d-th batch: %d\n", b, Nb[b]);
    }
    
    // Values of paramter
    double* mu = REAL(getListElement(args, "mu"));
    double* nu = REAL(getListElement(args, "nu"));
    double* delta = REAL(getListElement(args, "delta"));
    double* gamma = REAL(getListElement(args, "gamma"));
    double* phi = REAL(getListElement(args, "phi"));
    double* pi = REAL(getListElement(args, "pi"));
    
    SEXP Loglikelihood = PROTECT(allocVector(REALSXP, 1));
    nProtect++;
    double *like = REAL(Loglikelihood);
    
    //auxiliary variable
    SEXP logproby = PROTECT(allocVector(REALSXP, K));
    nProtect++;
    double *lpy = REAL(logproby);
    //lpy = log(pi_bk) + sum_g=1^G log(Pr(y_big|-))
    
    //int largest_x = INTEGER(getListElement(args, "largest.x"))[0];
    //SEXP logread0 = PROTECT(allocVector(REALSXP, largest_x));
    //nProtect ++;
    //double *lr0 = REAL(logread0);
    //lr0 = log(Pr(x_{big}=x| y_{big} = 0) x= 0,1,2, ... , 999
    
    int index_n = 0;
    int index, index_mu, index_nu;
    int read;
    double logmubikg, pbgk, lr0_temp, sum_lr0, lpy_max, sum_lpy;
    double loglike = 0.0;
    int x_max;
    for (int b = 0; b < B; b++) {
      for (int i = 0; i < Nb[b]; i++) {
        //if(index_n == 155){
        for (int k = 0; k < K; k++) {
          lpy[k] = log(pi[b * K + k]);
        }
        for (int k = 0; k < K; k++) {
          for (int j = 0; j < G; j++) {
            index = j * N + index_n;
            index_mu = k * G + j;
            index_nu = j * B + b;
            read = Reads[index];
            logmubikg = mu[index_mu] + nu[index_nu] + delta[index_n];
            pbgk = exp(logmubikg) / (exp(logmubikg) + phi[index_nu]);
            if (read>0) {
              lpy[k] = lpy[k] - log(1 + exp(-gamma[b] - gamma[b + B] * read));
              lpy[k] = lpy[k] + lgamma(phi[index_nu] + read) - lgamma(read + 1) - lgamma(phi[index_nu]) + read * log(pbgk) + phi[index_nu] * log(1 - pbgk);
            }
            else {
              
              x_max = (int)3 * exp(logmubikg);
              lr0_temp = phi[index_nu] * log(1 - pbgk); //x=0
              sum_lr0 = lr0_temp;
              
              for (int x = 1; x < x_max; x++) {
                lr0_temp = -gamma[b] - gamma[b + B] * x - log(1 + exp(-gamma[b] - gamma[b + B] * x));
                lr0_temp = lr0_temp + lgamma(phi[index_nu] + x) - lgamma(x + 1) - lgamma(phi[index_nu]) + x * log(pbgk) + phi[index_nu] * log(1 - pbgk);
                if (lr0_temp > sum_lr0) {
                  sum_lr0 = lr0_temp + log(1 + exp(sum_lr0 - lr0_temp));
                }
                else {
                  sum_lr0 = sum_lr0 + log(1 + exp(lr0_temp - sum_lr0));
                }
                
              }
              
              lpy[k] = lpy[k] + sum_lr0;
              //Rprintf("y %d %d %d = 0 and sum_lr0 is %f if the cell belongs to %d-th subtype.\n",b, i ,j ,sum_lr0, k);
            }
            
          }
          //Rprintf("lpy %d = %f for the %d-th cell.\n",k+1,lpy[k],index_n);
        }
        lpy_max = vec_max(lpy, K);
        sum_lpy = 0.0;
        for (int k = 0; k < K; k++) {
          sum_lpy = sum_lpy + exp(lpy[k] - lpy_max);
          //Rprintf("logproby[%d]=%f",k,lpy[k]);
        }
        loglike = loglike + lpy_max + log(sum_lpy);
        //}
        index_n = index_n + 1;
        //Rprintf("Finish the %d-th cell, lpy_max=%f, sum_lpy = %f, loglike = %f.\n", index_n, lpy_max, sum_lpy, loglike);
      }
    }
    
    like[0] = loglike;
    UNPROTECT(nProtect);
    return Loglikelihood;
  }
  
  SEXP Sample_X(SEXP args) {
    int nProtect = 0;
    //Rprintf("Input data.\n");
    
    // load list from R into Cpp 
    // load observed data
    int *Reads = INTEGER(getListElement(args, "Reads"));
    //[Y_11,Y_21,...,Y_N1,Y_12,Y_22,...,Y_NG]
    int B = INTEGER(getListElement(args, "BNum"))[0];
    int K = INTEGER(getListElement(args, "K"))[0];
    int* Nb = INTEGER(getListElement(args, "CNum"));
    int G = INTEGER(getListElement(args, "GNum"))[0];
    int N = 0;
    for (int b = 0; b < B; b++) {
      N = N + Nb[b];
      //Rprintf("The number of Cells in %d-th batch: %d\n", b, Nb[b]);
    }
    int iter_num;
    iter_num = INTEGER(getListElement(args, "n.iter"))[0];
    
    int* Dropout = INTEGER(getListElement(args, "z.init"));
    int* Trueread = INTEGER(getListElement(args, "x.init"));
    
    // Values of paramter
    double* alpha = REAL(getListElement(args, "alpha"));
    double* beta = REAL(getListElement(args, "beta"));
    double* nu = REAL(getListElement(args, "nu"));
    double* delta = REAL(getListElement(args, "delta"));
    double* gamma = REAL(getListElement(args, "gamma"));
    double* phi = REAL(getListElement(args, "phi"));
    int* w = INTEGER(getListElement(args, "w"));
    
    //Rprintf("Finish loading the posterior inference!\n");
    
    SEXP logmu = PROTECT(allocVector(REALSXP, N * G));
    nProtect++;
    double* log_mu = REAL(logmu);
    _update_logmu(B, Nb, N, G,
                  w, alpha, beta, nu, delta,//parameter
                  log_mu);//to be updated
    
    //Rprintf("Start sampling X and Z.\n");
    
    for (int iter = 0; iter < iter_num; iter++) {
      _update_zx(B, Nb, N, G,//dimension
                 gamma, phi, log_mu,//parameter
                 Reads, //observed data
                 Dropout, Trueread);//latent variable to be updated
      
      //Rprintf("Finish the %d-th sampling for X and Z.\n",iter);
    }
    
    UNPROTECT(nProtect);
    
  }

}
static const R_CallMethodDef CallEntries[] = {
  {"BZINBBUS",          (DL_FUNC) &BZINBBUS,          1},
  {"Cal_LogLike",       (DL_FUNC) &Cal_LogLike,       1},
  {"Sample_X",          (DL_FUNC) &Sample_X,          1},
  {NULL, NULL, 0}
};

void R_init_BUSseq(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}