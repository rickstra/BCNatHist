#include <RcppArmadillo.h>

// par[0] = A
// par[1] = B
// par[2] = delta
// par[3] = phi
// par[4] = mu
// par[5] = eta
// par[6] = beta_0
// par[7] = beta_s

//' @title CalcSymptCaseCumulative
//' @description Calculate the likelihood of being a symptomatic case, cumulative over age.
//' @export
// [[Rcpp::export]]
double CalcSymptCaseCumulative(double age, double v, arma::rowvec scr, 
                        arma::vec par, arma::colvec t, arma::colvec t_wts, 
                        double d0 = 0.5, double v0 = 0.06544985){
  t_wts *= 0.5 * age;
  t += 1;
  t *= 0.5 * age;
  
  double L = std::log(v / v0);
  double a = 1 / par[3];
  double b = a / par[4];
  
  double C = par[0] * par[1] * par[2] * std::pow(par[1] - par[0], par[2]) / 
    tgamma(a) * std::pow(b / L, a) / L;
  arma::colvec x_t = arma::exp((par[1] - par[0]) * t);
  arma::mat d = d0 * arma::exp(((1 / (age - t)) * (scr - age) + 1) * (L / 3));
  arma::mat ix = arma::conv_to<arma::mat>::from(d > d0);
  
  arma::colvec e  = arma::log(age - t) * a + 
                    par[1] * par[2] * t - 
                    (age - t) * b / L - 
              arma::sum(arma::log(1 + arma::exp(par[7] + par[6] * d)) % ix , 1); 
  
  double out = C * arma::sum(t_wts % 
    (1 - arma::exp(par[5] * (age - t) * (v0 - v) / L)) % 
    arma::exp(e) % (1 - x_t) % pow(par[1] * x_t - par[0], -par[2] - 1));
  
  return out;
}

//' @title CalcFutureScreeningSens
//' @description Calculate the probability of being detected at a future screening. 
//' Note: the screening at \code{age} should not be included in the screening history \code{scr}.
//' @export
// [[Rcpp::export]]
double CalcFutureScreeningSens(double age, arma::rowvec scr, arma::vec par,
                    arma::colvec t, arma::colvec t_wts,
                    arma::rowvec r, arma::rowvec r_wts,
                    double d0 = 0.5, double v0 = 0.06544985) {
  t_wts *= 0.5 * age;
  t += 1;
  t *= 0.5 * age;
  
  int n_scr = scr.n_elem;
  arma::colvec ones(t.n_elem, arma::fill::ones);
  
  arma::rowvec r_inv = 1 / r;
  arma::colvec x_t = arma::exp((par[1] - par[0]) * t);
  
  double a = 1 / par[3];
  double b = a / par[4];
  
  double C = -par[0] * par[1] * par[2] * std::pow(par[1] - par[0], par[2]) / 
    tgamma(a) * std::pow(b, a);
  
  arma::colvec e_t  = (par[1] * par[2]) * t - 
                      arma::log(x_t - 1) - 
                      (par[2] + 1) * arma::log(par[1] * x_t - par[0]);
  arma::rowvec e_r  = (par[5] * v0 - b + 1) * r + // +1r for gauss-laguerre
                      (a - 1) * arma::log(r);
  arma::mat e_tr = (ones * ((-par[5] * v0) * r)) % arma::exp((age - t) * r_inv);
  
  // the positive screen
  arma::mat d = d0 * arma::exp((age - t) * (r_inv / 3));
  arma::mat ix = arma::conv_to<arma::mat>::from(d > d0);
  e_tr += ((par[7] + par[6] * d) - 
    arma::log(1 + arma::exp(par[7] + par[6] * d))) % ix;
  
  // the negative screens
  if (n_scr > 0) {
    for (int i=0; i<n_scr; i++) {
      d = d0 * arma::exp((scr[i] - t) * (r_inv / 3));
      ix = arma::conv_to<arma::mat>::from(d > d0);
      e_tr -= arma::log(1 + arma::exp(par[7] + par[6] * d)) % ix;
    }
  }
  double out = C * arma::sum(arma::sum((t_wts * r_wts) %
    (arma::exp(e_t) * arma::exp(e_r)) % arma::exp(e_tr), 0));
  
  return out;
}
