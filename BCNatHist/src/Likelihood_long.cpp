#include <RcppArmadillo.h>

// par[0] = A
// par[1] = B
// par[2] = delta
// par[3] = phi
// par[4] = mu
// par[5] = eta
// par[6] = beta0
// par[7] != beta1

//' @title EvalCaseIntegralLong
//' @description Evaluate the integral for BC cases.
//' @export
// [[Rcpp::export]]
double EvalCaseIntegralLong(double age, double v, arma::rowvec scr, arma::vec par, 
                       arma::rowvec scr_par, 
                       arma::colvec t, arma::colvec t_wts, 
                       double d0 = 0.5, double v0 = 0.06544985){
  t_wts *= 0.5 * age;
  t += 1;
  t *= 0.5 * age;
  arma::colvec ones(t.n_elem, arma::fill::ones);
  
  double L = std::log(v / v0);
  double a = 1 / par[3];
  double b = a / par[4];
  
  double C = par[0] * par[1] * par[2] * std::pow(par[1] - par[0], par[2]) / 
    tgamma(a) * std::pow(b / L, a) / L;
  
  arma::colvec x_t = arma::exp((par[1] - par[0]) * t);
  
  arma::mat d = d0 * arma::exp((((1 / (age - t)) * (scr - age)) + 1) * (L / 3));
  arma::mat ix = arma::conv_to<arma::mat>::from(d > d0);
  
  arma::colvec e  = arma::log(age - t) * a + 
    par[1] * par[2] * t - 
    (age - t) * ((b + par[5] * (v - v0)) / L) - 
    arma::sum(arma::log(1 + arma::exp(ones * scr_par + par[6] * d)) % ix , 1);
  
  return C * arma::sum(t_wts % arma::exp(e) % (1 - x_t) %
                      pow(par[1] * x_t - par[0], -par[2] - 1));
  
}

//' @title EvalCensIntegralLong
//' @description Evaluate the double integral for censored.
//' @export
// [[Rcpp::export]]
double EvalCensIntegralLong(double age, arma::rowvec scr, arma::vec par, 
                            arma::rowvec scr_par,
                            arma::colvec t, arma::colvec t_wts,
                            arma::rowvec r, arma::rowvec r_wts,
                            double d0 = 0.5, double v0 = 0.06544985){
  t_wts *= 0.5 * age;
  t += 1;
  t *= 0.5 * age;

  int n_scr = scr.n_elem;
  arma::colvec ones(t.n_elem, arma::fill::ones);
  arma::rowvec r_inv = 1 / r;
  arma::mat x_t = arma::exp((par[1] - par[0]) * t);
  double a = 1 / par[3];
  double b = a / par[4];

  double C = -par[0] * par[1] * par[2] * std::pow(par[1] - par[0], par[2]) * 
    std::pow(b, a) / tgamma(a); 
  
  arma::colvec e_t  = (par[1] * par[2]) * t + arma::log(x_t - 1) -
    (par[2] + 1) * arma::log(par[1] * x_t - par[0]);
  
  arma::rowvec e_r  = (par[5] * v0 + 1) * r + // +1r for gauss-laguerre
    (a - 1) * arma::log(r) - (r * b); 
  
  arma::mat e_tr = (ones * ((-par[5] * v0) * r)) % arma::exp((age - t) * r_inv);
  if (n_scr > 0) {
    arma::mat d, ix;
    for (int i=0; i<n_scr; i++) {
      d = d0 * arma::exp((scr[i] - t) * (r_inv / 3));
      ix = arma::conv_to<arma::mat>::from(d > d0);
      e_tr -= arma::log(1 + arma::exp(scr_par[i] + par[6] * d)) % ix;
    }
  }
  return C * arma::sum(arma::sum((t_wts * r_wts) %
           (arma::exp(e_t) * arma::exp(e_r)) % arma::exp(e_tr), 0));
}

//' @title CalcScreenCaseLong
//' @description Calculate the joint likelihood of a single screen-detected case.
//' @export
// [[Rcpp::export]]
double CalcScreenCaseLong(double age, double v, arma::rowvec scr, arma::vec par, 
                          arma::rowvec scr_par,
                          arma::colvec t_pts, arma::colvec t_wts,
                          double d0 = 0.5, double v0 = 0.06544985){
  double out = scr_par[0] + par[6] * d0 * std::exp(std::log(v / v0) / 3) - 
    std::log(v) + 
    std::log(EvalCaseIntegralLong(age, v, scr, par, scr_par, t_pts, t_wts, d0, v0));
  return std::exp(out);
}

//' @title CalcSymptCaseLong
//' @description Calculate the joint likelihood of a single symptomatic case.
//' @export
// [[Rcpp::export]]
double CalcSymptCaseLong(double age, double v, arma::rowvec scr, arma::vec par, 
                         arma::rowvec scr_par,
                         arma::colvec t_pts, arma::colvec t_wts,
                         double d0 = 0.5, double v0 = 0.06544985){
  return par[5] * EvalCaseIntegralLong(age, v, scr, par, scr_par, 
                                   t_pts, t_wts, d0, v0);
}

//' @title CalcCensCaseLong
//' @description Calculate the likelihood of a single censored.
//' @export
// [[Rcpp::export]]
double CalcCensCaseLong(double age, arma::rowvec scr, arma::vec par, 
                        arma::rowvec scr_par,
                        arma::colvec t_pts, arma::colvec t_wts,
                        arma::rowvec r_pts, arma::rowvec r_wts,
                        double d0 = 0.5, double v0 = 0.06544985){
  return std::pow((par[1] - par[0]) * std::exp(par[1] * age) /
           (par[1] * std::exp((par[1] - par[0]) * age) - par[0]), par[2]) + 
           EvalCensIntegralLong(age, scr, par, scr_par, 
                                t_pts, t_wts, r_pts, r_wts, d0, v0);
}

//' @title CalcIndividualLikelihoodLong
//' @description Determine case and mode status, and calculate the corresponding joint likelihood.
//' @export
// [[Rcpp::export]]
double CalcIndividualLikelihoodLong(int is_case, int is_scr, double age, 
                                    double v, arma::rowvec scr, 
                                    double entry, arma::rowvec e_scr,
                                    arma::vec par, arma::rowvec scr_par, 
                                    arma::rowvec e_scr_par,
                                    arma::colvec t_pts, arma::colvec t_wts,
                                    arma::rowvec r_pts, arma::rowvec r_wts,
                                    double d0 = 0.5, double v0 = 0.06544985){
  double out;
  if (is_case == 0) {
    out = CalcCensCaseLong(age, scr, par, scr_par, 
                           t_pts, t_wts, r_pts, r_wts, d0, v0);
  } else if (is_scr == 1) {
    out = CalcScreenCaseLong(age, v, scr, par, scr_par, t_pts, t_wts, d0, v0);
  } else {
    out = CalcSymptCaseLong(age, v, scr, par, scr_par, t_pts, t_wts, d0, v0);
  }
  return out / CalcCensCaseLong(entry, e_scr, par, e_scr_par, 
                          t_pts, t_wts, r_pts, r_wts, d0, v0);
}