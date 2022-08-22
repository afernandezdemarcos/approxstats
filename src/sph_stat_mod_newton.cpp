// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double start = 0.001;
double end = 0.25;
double by = 0.001;
arma::vec alpha_grid = arma::linspace(start, end, ((end-start)/by)+1);

//                                                           //
//        Modification of hyperspherical statistics          //
//                                                           //

//[[Rcpp::export]]
double sph_stat_PCvM_mod(double statistic, int n, double alpha, int p){
  
  double f = (1 + (+0.1130/sqrt(p) -0.5415/p) / n
                + (+0.1438/sqrt(p)) / (n*sqrt(alpha))
                + (-0.0031/sqrt(p)) / (n*alpha));
                
  return(statistic * f);
}

//[[Rcpp::export]]
double sph_stat_PAD_mod(double statistic, int n, double alpha, int p, bool PAD_2 = false){
  
  double f;
  
  if (PAD_2){
    f = (1 - 0.0751/n + 0.0692/(n*sqrt(alpha)) - 0.0014/(n*alpha));
  }else{
    f = (1 + (+0.0978/sqrt(p) -0.3596/p) / n
                  + (0 +0.1126/sqrt(p)) / (n*sqrt(alpha))
                  + (0 -0.0025/sqrt(p)) / (n*alpha));
                  
  }
  
  return(statistic * f);
}

//[[Rcpp::export]]
double sph_stat_Bakshaev_mod(double statistic, int n, double alpha, int p){
  
  double f = (1 + (+0.1189/sqrt(p) -0.5838/p) / n
                + (0 +0.1210/sqrt(p) +0.0385/p) / (n*sqrt(alpha))
                + (0 -0.0030/sqrt(p)) / (n*alpha)
  );
  
  return(statistic * f);
}

//                                                           //
//       Asymptotic quantiles of circular statistics         //
//                                                           //

// NumericVector PCvM_asymp_q = {};
// PCvM_asymp_q.attr("dim") = Dimension(, );

//                  Newton-Raphson methods                   //

//[[Rcpp::export]]
double f_quadratic(double alpha, double T_n, int n, 
                   double beta_n, double beta_n_alpha_1_2, double beta_nalpha, 
                   double T_inf2, double T_inf1, double alpha2, double alpha1){
  
  double m = (T_inf2 - T_inf1)/(alpha2 - alpha1);
  double h = (T_n - T_inf2 + m*alpha2 + beta_n*T_n/n) + (beta_n_alpha_1_2*T_n/n)/sqrt(alpha) + (beta_nalpha*T_n/n)/alpha - m*alpha;
  
  return(h);
  
}

//[[Rcpp::export]]
double de_quadratic(double alpha, double T_n, int n, 
                    double beta_n, double beta_n_alpha_1_2, double beta_nalpha, 
                    double T_inf2, double T_inf1, double alpha2, double alpha1){
  
  double m = (T_inf2 - T_inf1)/(alpha2 - alpha1);
  double der = -beta_nalpha/(alpha*alpha) -1/2*(beta_n_alpha_1_2*T_n/n)/sqrt(alpha*alpha*alpha) - m;
  
  return(der);
  
}

//[[Rcpp::export]]
double newton_raphson_quadratic(double T_n, int n, 
                                double beta_n, double beta_n_alpha_1_2, double beta_nalpha, 
                                double T_inf2, double T_inf1, double alpha2, double alpha1, double thr = 1e-7){
  
  auto fobj = [T_n, n, beta_n, beta_n_alpha_1_2, beta_nalpha, T_inf2, T_inf1, alpha2, alpha1](double alpha){
    return f_quadratic(alpha, T_n, n, beta_n, beta_n_alpha_1_2, beta_nalpha, T_inf2, T_inf1, alpha2, alpha1);
  };
  auto derobj = [T_n, n, beta_n, beta_n_alpha_1_2, beta_nalpha, T_inf2, T_inf1, alpha2, alpha1](double alpha){
    return de_quadratic(alpha, T_n, n, beta_n, beta_n_alpha_1_2, beta_nalpha, T_inf2, T_inf1, alpha2, alpha1);
  };
  
  double z;
  double a = alpha1 + (alpha2-alpha1)/2;
  int e = 0;
  do{
    e++;
    z = a - (fobj(a)/derobj(a));
    a = z;
  }while((abs(fobj(z))>thr) && (e < 7));
  
  if(a < 0){
    a = 0;
  }
  
  return(a);
  
}

//                                                           //
//       Exact distribution of hyperspherical statistics     //
//                                                           //

//[[Rcpp::export]]
double p_sph_stat_PCvM(double T_n, int n, int p, NumericVector T_inf){
  int m = alpha_grid.size();
  int i;
  double alpha;
  double T_mod;
  double sqrt_p = sqrt(p);
  double beta_n = 0.1130/sqrt_p - 0.5415/p;
  double beta_n_alpha_1_2 = 0.1438/sqrt_p;
  double beta_nalpha = -0.0031/sqrt_p;

  for(i = 0; i < m; i++){
    alpha = alpha_grid[i];
    T_mod = sph_stat_PCvM_mod(T_n, n, alpha, p);
    
    if(T_mod > T_inf[i]){
      if(i > 0){
        alpha = newton_raphson_quadratic(T_n, n, beta_n, beta_n_alpha_1_2, beta_nalpha,
                                         T_inf[i-1], T_inf[i], alpha_grid[i-1], alpha);
      }else{
        alpha = newton_raphson_quadratic(T_n, n, beta_n, beta_n_alpha_1_2, beta_nalpha,
                                         T_inf[i], T_inf[i+1], alpha, alpha_grid[i+1]);
      }
      return(alpha);
    }
  }
  return(NumericVector::get_na());
}

//[[Rcpp::export]]
double p_sph_stat_PAD(double T_n, int n, int p, NumericVector T_inf, bool PAD_2 = false){
  int m = alpha_grid.size();
  int i;
  double alpha;
  double T_mod;
  double sqrt_p = sqrt(p);
  double beta_n;
  double beta_n_alpha_1_2;
  double beta_nalpha;
  
  if(PAD_2){
    beta_n = -0.0751;
    beta_n_alpha_1_2 = 0.0692;
    beta_nalpha = -0.0014;
  }else{
    beta_n = 0.0978/sqrt_p - 0.3596/p;
    beta_n_alpha_1_2 = 0.1126/sqrt_p;
    beta_nalpha = -0.0025/sqrt_p;
  }
  
  for(i = 0; i < m; i++){
    alpha = alpha_grid[i];
    T_mod = sph_stat_PAD_mod(T_n, n, alpha, p, PAD_2 = PAD_2);
    
    if(T_mod > T_inf[i]){
      if(i > 0){
        alpha = newton_raphson_quadratic(T_n, n, beta_n, beta_n_alpha_1_2, beta_nalpha,
                                         T_inf[i-1], T_inf[i], alpha_grid[i-1], alpha);
      }else{
        alpha = newton_raphson_quadratic(T_n, n, beta_n, beta_n_alpha_1_2, beta_nalpha,
                                         T_inf[i], T_inf[i+1], alpha, alpha_grid[i+1]);
      }
      return(alpha);
    }
  }
  return(NumericVector::get_na());
}

//[[Rcpp::export]]
double p_sph_stat_Bakshaev(double T_n, int n, int p, NumericVector T_inf){
  int m = alpha_grid.size();
  int i;
  double alpha;
  double T_mod;
  double sqrt_p = sqrt(p);
  double beta_n = 0.1189/sqrt_p - 0.5838/p;
  double beta_n_alpha_1_2 = 0.1210/sqrt_p + 0.0385/p;
  double beta_nalpha = -0.0030/sqrt_p;
  
  for(i = 0; i < m; i++){
    alpha = alpha_grid[i];
    T_mod = sph_stat_Bakshaev_mod(T_n, n, alpha, p);
    
    if(T_mod > T_inf[i]){
      if(i > 0){
        alpha = newton_raphson_quadratic(T_n, n, beta_n, beta_n_alpha_1_2, beta_nalpha,
                                         T_inf[i-1], T_inf[i], alpha_grid[i-1], alpha);
      }else{
        alpha = newton_raphson_quadratic(T_n, n, beta_n, beta_n_alpha_1_2, beta_nalpha,
                                         T_inf[i], T_inf[i+1], alpha, alpha_grid[i+1]);
      }
      return(alpha);
    }
  }
  return(NumericVector::get_na());
}