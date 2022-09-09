#include <RcppArmadillo.h>

// par[0] = A
// par[1] = B
// par[2] = delta
// par[3] = phi
// par[4] = mu
// par[5] = eta
// par[6] = beta_0
// par[7] = beta_s

// [[Rcpp::export]]
double f_CaseIntegral(double age, double v, arma::rowvec scr, 
                    arma::vec par, arma::colvec t, arma::colvec t_wts, 
                    double d0 = 0.5, double v0 = 0.06544985)
{
  double out, L, C, a, b;
  arma::colvec e, x_t;
  arma::mat d, ix;
        
  t_wts *= 0.5 * age;
  t += 1;
  t *= 0.5 * age;

  L = std::log(v / v0);
  a = 1 / par[3];
  b = a / par[4];
  
  C = par[0] * par[1] * par[2] * std::pow(par[1] - par[0], par[2]) / 
      tgamma(a) * std::pow(b / L, a) / L;
  x_t = arma::exp((par[1] - par[0]) * t);
  d = d0 * arma::exp((((1 / (age - t)) * (scr - age)) + 1) * (L / 3));
  ix = arma::conv_to<arma::mat>::from(d > d0);

  e  = arma::log(age - t) * a; 
  e += par[1] * par[2] * t; 
  e -= (age - t) * ((b + par[5] * (v - v0)) / L);
  e -= arma::sum(arma::log(1 + arma::exp(par[7] + par[6] * d)) % ix , 1); 

  out = C * arma::sum(t_wts % arma::exp(e) % (1 - x_t) %
        pow(par[1] * x_t - par[0], -par[2] - 1));

  return out;
}

// [[Rcpp::export]]
double f_CensIntegral(double age, arma::rowvec scr, arma::vec par,
                    arma::colvec t, arma::colvec t_wts,
                    arma::rowvec r, arma::rowvec r_wts,
                    double d0 = 0.5, double v0 = 0.06544985)
{
  double out, C, a, b;
  arma::colvec e_t, x_t;
  arma::rowvec e_r, r_inv;
  arma::mat e_tr, d, ix;

  arma::colvec ones(t.n_elem, arma::fill::ones);
  int n_scr = scr.n_elem; 

  t_wts *= 0.5 * age;
  t += 1;
  t *= 0.5 * age; 

  r_inv = 1 / r;
  x_t = arma::exp((par[1] - par[0]) * t);
  a = 1 / par[3];
  b = a / par[4];
  
  C = -par[0] * par[1] * par[2] * std::pow(par[1] - par[0], par[2]) * 
      std::pow(b, a) / tgamma(a); 

  e_t  = (par[1] * par[2]) * t; 
  e_t += arma::log(x_t - 1); 
  e_t -= (par[2] + 1) * arma::log(par[1] * x_t - par[0]);
  
  e_r  = (par[5] * v0 + 1) * r; // +1r for gauss-laguerre
  e_r += (a - 1) * arma::log(r); 
  e_r -= (r * b); 
  
  e_tr = (ones * ((-par[5] * v0) * r)) % arma::exp((age - t) * r_inv);
  if (n_scr > 0) {
    for (int i=0; i<n_scr; ++i) {
      d = d0 * arma::exp((scr[i] - t) * (r_inv / 3));
      ix = arma::conv_to<arma::mat>::from(d > d0);
      e_tr -= arma::log(1 + arma::exp(par[7] + par[6] * d)) % ix;
    }
  }

  out = C * arma::sum(arma::sum((t_wts * r_wts) %
        (arma::exp(e_t) * arma::exp(e_r)) % arma::exp(e_tr), 0)); 

  return out; 
}

// [[Rcpp::export]]
double f_ScreenCase(double age, double v, arma::rowvec scr, arma::vec par, 
                  arma::colvec t_pts, arma::colvec t_wts, 
                  double d0 = 0.5, double v0 = 0.06544985) 
{
  double out;
  out = par[7] + par[6] * d0 * std::exp(std::log(v / v0) / 3) - std::log(v);
  out += std::log(f_CaseIntegral(age, v, scr, par, t_pts, t_wts, d0, v0));
  out = std::exp(out);
  return out;
}

// [[Rcpp::export]]
double f_SymptCase(double age, double v, arma::rowvec scr, arma::vec par, 
                 arma::colvec t_pts, arma::colvec t_wts, 
                  double d0 = 0.5, double v0 = 0.06544985) 
{
  double out;
  out = par[5];
  out *= f_CaseIntegral(age, v, scr, par, t_pts, t_wts, d0, v0);
  return out;
}

// [[Rcpp::export]]
double f_CensCase(double age, arma::rowvec scr, arma::vec par, 
                arma::colvec t_pts, arma::colvec t_wts, 
                arma::rowvec r_pts, arma::rowvec r_wts,
                double d0 = 0.5, double v0 = 0.06544985) 
{
  double out;
  out = std::pow((par[1] - par[0]) * std::exp(par[1] * age) /
        (par[1] * std::exp((par[1] - par[0]) * age) - par[0]), par[2]);
  out += f_CensIntegral(age, scr, par, t_pts, t_wts, r_pts, r_wts, d0, v0);
  return out;
}

// [[Rcpp::export]]
double IndL(int is_case, int is_scr, double age, double v, arma::rowvec scr,
            double entry, arma::rowvec e_scr, arma::vec par,
            arma::colvec t_pts, arma::colvec t_wts,
            arma::rowvec r_pts, arma::rowvec r_wts,
            double d0 = 0.5, double v0 = 0.06544985) 
{
  double out;
  if(is_case == 0) {
    out = f_CensCase(age, scr, par, t_pts, t_wts, r_pts, r_wts, d0, v0);
  } else if (is_scr == 1) {
    out = f_ScreenCase(age, v, scr, par, t_pts, t_wts, d0, v0);
  } else {
    out = f_SymptCase(age, v, scr, par, t_pts, t_wts, d0, v0);
  }

  out /= f_CensCase(entry, e_scr, par, t_pts, t_wts, r_pts, r_wts, d0, v0);
  return out; 
}

////////////////////////////////////////////////////////////////////////////////
// Onset-Dependent Growth functions
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double f_CaseIntegral_odg(double age, double v, arma::rowvec scr, 
                          arma::vec par, arma::colvec t, arma::colvec t_wts, 
                          double d0 = 0.5, double v0 = 0.06544985)
{
  double out;
  
  t_wts *= 0.5 * age;
  t += 1;
  t *= 0.5 * age;
  
  arma::colvec e;
  arma::mat d;
  arma::mat ix;
  
  double a = 1 / par[3];
  arma::colvec b = arma::exp(-par[8] * t) / (par[3] * par[4]);
  
  double L = std::log(v / v0);
  double C = par[0] * par[1] * par[2] * std::pow(par[1] - par[0], par[2]) / 
    tgamma(a) / std::pow(L, a + 1);
  arma::colvec x_t = arma::exp((par[1] - par[0]) * t);
  
  d = d0 * arma::exp((((1 / (age - t)) * (scr - age)) + 1) * (L / 3));
  ix = arma::conv_to<arma::mat>::from(d > d0);
  
  e  = arma::log(age - t) * a; 
  e += par[1] * par[2] * t; 
  e -= ((age - t) / L) % (b + par[5] * (v - v0));
  e += arma::log(b) * a;
  e -= arma::sum(arma::log(1 + arma::exp(par[7] + par[6] * d)) % ix , 1);
  
  out = C * arma::sum(t_wts % arma::exp(e) % (1 - x_t) %
    pow(par[1] * x_t - par[0], -par[2] - 1));
  
  return out;
}

// [[Rcpp::export]]
double f_CensIntegral_odg(double age, arma::rowvec scr, arma::vec par,
                          arma::colvec t, arma::colvec t_wts,
                          arma::rowvec r, arma::rowvec r_wts,
                          double d0 = 0.5, double v0 = 0.06544985)
{
  double out;
  
  
  t_wts *= 0.5 * age;
  t += 1;
  t *= 0.5 * age;
  
  int n_scr = scr.n_elem;
  arma::colvec ones(t.n_elem, arma::fill::ones);
  arma::colvec e_t;
  arma::rowvec e_r;
  arma::mat e_tr;
  arma::mat d;
  arma::mat ix;
  
  arma::rowvec r_inv = 1 / r;
  arma::colvec x_t = arma::exp((par[1] - par[0]) * t);
  
  double a = 1 / par[3];
  arma::colvec b = arma::exp(-par[8] * t) / (par[3] * par[4]);
  
  double C = -par[0] * par[1] * par[2] * std::pow(par[1] - par[0], par[2]) / 
    tgamma(a);
  
  e_t  = (par[1] * par[2]) * t + arma::log(b) * a ;
  e_t += arma::log(x_t - 1); 
  e_t -= (par[2] + 1) * arma::log(par[1] * x_t - par[0]);
  e_r  = (par[5] * v0 + 1) * r; // +1r for gauss-laguerre
  e_r += (a - 1) * arma::log(r); 
  // 
  e_tr = -b * r; 
  e_tr+= (ones * ((-par[5] * v0) * r)) % arma::exp((age - t) * r_inv);
  
  if (n_scr > 0) {
    for (int i=0; i<n_scr; i++) {
      d = d0 * arma::exp((scr[i] - t) * (r_inv / 3));
      ix = arma::conv_to<arma::mat>::from(d > d0);
      e_tr -= arma::log(1 + arma::exp(par[7] + par[6] * d)) % ix;
    }
  }
  
  out = C * arma::sum(arma::sum((t_wts * r_wts) %
    (arma::exp(e_t) * arma::exp(e_r)) % arma::exp(e_tr), 0));

  return out;
  
}

// [[Rcpp::export]]
double f_ScreenCase_odg(double age, double v, arma::rowvec scr, arma::vec par, 
                        arma::colvec t_pts, arma::colvec t_wts, 
                        double d0 = 0.5, double v0 = 0.06544985) 
{
  double out = std::exp(par[7] + par[6] * d0 * std::exp(std::log(v / v0) / 3)) / v;
  out *= f_CaseIntegral_odg(age, v, scr, par, t_pts, t_wts, d0, v0);
  return out;
}
// [[Rcpp::export]]
double f_SymptCase_odg(double age, double v, arma::rowvec scr, arma::vec par, 
                       arma::colvec t_pts, arma::colvec t_wts, 
                       double d0 = 0.5, double v0 = 0.06544985) 
{
  double out = par[5];
  out *= f_CaseIntegral_odg(age, v, scr, par, t_pts, t_wts, d0, v0);
  return out;
}
// [[Rcpp::export]]
double f_CensCase_odg(double age, arma::rowvec scr, arma::vec par, 
                  arma::colvec t_pts, arma::colvec t_wts, 
                  arma::rowvec r_pts, arma::rowvec r_wts,
                  double d0 = 0.5, double v0 = 0.06544985) 
{
  double out = std::pow((par[1] - par[0]), par[2]) * 
    std::exp(par[1] * par[2] * age) /
      pow(par[1] * std::exp((par[1] - par[0]) * age) - par[0], par[2]);
  out += f_CensIntegral_odg(age, scr, par, t_pts, t_wts, r_pts, r_wts, d0, v0);
  return out;
}
// [[Rcpp::export]]
double IndL_odg(int is_case, int is_scr, double age, double v, arma::rowvec scr,
            double entry, arma::rowvec e_scr, arma::vec par,
            arma::colvec t_pts, arma::colvec t_wts,
            arma::rowvec r_pts, arma::rowvec r_wts,
            double d0 = 0.5, double v0 = 0.06544985) {
  double out;
  if(is_case == 0) {
    out = f_CensCase_odg(age, scr, par, t_pts, t_wts, r_pts, r_wts, d0, v0);
  } else if (is_scr == 1) {
    out = f_ScreenCase_odg(age, v, scr, par, t_pts, t_wts, d0, v0);
  } else {
    out = f_SymptCase_odg(age, v, scr, par, t_pts, t_wts, d0, v0);
  }
  
  out /= f_CensCase_odg(entry, e_scr, par, t_pts, t_wts, r_pts, r_wts, d0, v0);
  return out;
}

////////////////////////////////////////////////////////////////////////////////
// Other functions
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double f_CaseIntegral2(double age, double v, arma::rowvec scr, 
                      arma::vec par, arma::colvec t, arma::colvec t_wts, 
                      double d0 = 0.5, double v0 = 0.06544985)
{
  double out;
  
  
  d0 = std::pow(d0, par[8]);
  
  t_wts *= 0.5 * age;
  t += 1;
  t *= 0.5 * age;
  
  arma::colvec e;
  arma::mat d;
  arma::mat ix;
  
  double L = std::log(v / v0);
  double C = par[0] * par[1] * par[2] * std::pow(par[1] - par[0], par[2]) / 
    tgamma(1 / par[3]) / std::pow(par[3] * par[4] * L, 1 / par[3]) / L;
  arma::colvec x_t = arma::exp((par[1] - par[0]) * t);
  
  // d = std::log(d0) + (((1 / (age - t)) * (scr - age)) + 1) * (L / 3);
  // ix = arma::conv_to<arma::mat>::from(d > logd0);
  d = d0 * arma::exp((((1 / (age - t)) * (scr - age)) + 1) * (L * par[8] / 3));
  ix = arma::conv_to<arma::mat>::from(d > d0);
  
  e  = arma::log(age - t) / par[3]; 
  e += par[1] * par[2] * t; 
  e -= ((age - t) / L) * (1 / (par[4] * par[3]) + par[5] * (v - v0));
  // e -= arma::sum(arma::log(1 + arma::exp(par[7] + par[6] * d)) % ix , 1);
  e -= arma::sum(arma::log(1 + arma::exp(par[7] + par[6] * d)) % ix , 1);
  
  out = C * arma::sum(t_wts % arma::exp(e) % (1 - x_t) %
    pow(par[1] * x_t - par[0], -par[2] - 1));
  
  return out;
  
}


// [[Rcpp::export]]
double f_CensIntegral2(double age, arma::rowvec scr, arma::vec par,
                      arma::colvec t, arma::colvec t_wts,
                      arma::rowvec r, arma::rowvec r_wts,
                      double d0 = 0.5, double v0 = 0.06544985)
{
  double out;
  
  d0 = std::pow(d0, par[8]);
  
  t_wts *= 0.5 * age;
  t += 1;
  t *= 0.5 * age;
  
  int n_scr = scr.n_elem;
  arma::colvec ones(t.n_elem, arma::fill::ones);
  arma::colvec e_t;
  arma::rowvec e_r;
  arma::mat e_tr;
  arma::mat d;
  arma::mat ix;
  
  arma::rowvec r_inv = 1 / r;
  arma::colvec x_t = arma::exp((par[1] - par[0]) * t);
  
  double C = -par[0] * par[1] * par[2] * std::pow(par[1] - par[0], par[2]) / 
    tgamma(1 / par[3]) / std::pow(par[3], 1 / par[3]);
  
  e_t  = (par[1] * par[2]) * t - std::log(par[4]) / par[3];
  e_t += arma::log(x_t - 1); 
  e_t -= (par[2] + 1) * arma::log(par[1] * x_t - par[0]);
  e_r  = (par[5] * v0 + 1) * r; // +1r for gauss-laguerre
  e_r += (1 / par[3] - 1) * arma::log(r); 
  e_r -= (r / (par[3] * par[4])); 
  e_tr = (ones * ((-par[5] * v0) * r)) % arma::exp((age - t) * r_inv);
  if (n_scr > 0) {
    for (int i=0; i<n_scr; i++) {
      // d = std::log(d0) + (scr[i] - t) * (r_inv / 3);
      d = d0 * arma::exp((scr[i] - t) * (r_inv * (par[8] / 3)));
      // ix = arma::conv_to<arma::mat>::from(d > logd0);
      ix = arma::conv_to<arma::mat>::from(d > d0);
      // e_tr -= arma::log(1 + arma::exp(par[7] + par[6] * d)) % ix;
      e_tr -= arma::log(1 + arma::exp(par[7] + par[6] * d)) % ix;
    }
  }
  
  out = C * arma::sum(arma::sum((t_wts * r_wts) %
    (arma::exp(e_t) * arma::exp(e_r)) % arma::exp(e_tr), 0));
  
  return out;
}

// [[Rcpp::export]]
double f_ScreenCase2(double age, double v, arma::rowvec scr, arma::vec par, 
                    arma::colvec t_pts, arma::colvec t_wts, 
                    double d0 = 0.5, double v0 = 0.06544985) 
{
  if (scr.n_elem == 1) {
    scr = {-1, -1};
  }
  double out = f_CaseIntegral(age, v, scr(arma::span(1, scr.n_elem-1)), par, t_pts, t_wts, d0, v0);
  d0 = std::pow(d0, par[8]);
  out *= 1- 1 / (1 + std::exp(par[7] + par[6] * (d0 * std::exp(std::log(v / v0) * par[8] / 3))));
  out /= v;
  return out;
  
}

// [[Rcpp::export]]
double f_SymptCase2(double age, double v, arma::rowvec scr, arma::vec par, 
                   arma::colvec t_pts, arma::colvec t_wts, 
                   double d0 = 0.5, double v0 = 0.06544985) 
{
  double out = par[5];
  out *= f_CaseIntegral2(age, v, scr, par, t_pts, t_wts, d0, v0);
  return out;
}


// [[Rcpp::export]]
double f_CensCase2(double age, arma::rowvec scr, arma::vec par, 
                  arma::colvec t_pts, arma::colvec t_wts, 
                  arma::rowvec r_pts, arma::rowvec r_wts,
                  double d0 = 0.5, double v0 = 0.06544985) 
{
  double out = std::pow((par[1] - par[0]), par[2]) * 
    std::exp(par[1] * par[2] * age) /
      pow(par[1] * std::exp((par[1] - par[0]) * age) - par[0], par[2]);
  out += f_CensIntegral2(age, scr, par, t_pts, t_wts, r_pts, r_wts, d0, v0);
  return out;
}
// [[Rcpp::export]]
double IndL2(int is_case, int is_scr, double age, double v, arma::rowvec scr,
            double entry, arma::rowvec e_scr, arma::vec par,
            arma::colvec t_pts, arma::colvec t_wts,
            arma::rowvec r_pts, arma::rowvec r_wts,
            double d0 = 0.5, double v0 = 0.06544985) {
  double out;
  if(is_case == 0) {
    out = f_CensCase2(age, scr, par, t_pts, t_wts, r_pts, r_wts, d0, v0);
  } else if (is_scr == 1) {
    out = f_ScreenCase2(age, v, scr, par, t_pts, t_wts, d0, v0);
  } else {
    out = f_SymptCase2(age, v, scr, par, t_pts, t_wts, d0, v0);
  }
  
  out /= f_CensCase2(entry, e_scr, par, t_pts, t_wts, r_pts, r_wts, d0, v0);
  return out;
}

////////////////////////////////////////////////////////////////////////////
