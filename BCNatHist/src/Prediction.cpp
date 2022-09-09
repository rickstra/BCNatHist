#include <RcppArmadillo.h>

// par[0] = A
// par[1] = B
// par[2] = delta
// par[3] = phi
// par[4] = mu
// par[5] = eta
// par[6] = beta0
// par[7] != beta1

// [[Rcpp::export]]
double f_ScreenSizeInterval(double age, double v_cut, arma::rowvec scr, 
                            arma::vec par, arma::colvec t, arma::colvec t_wts, 
                            arma::rowvec v, arma::rowvec v_wts, 
                            double d0 = 0.5, double v0 = 0.06544985)
{
  double out;
  int n_scr = scr.n_elem;
  
  t_wts *= 0.5 * age;
  t += 1;
  t *= 0.5 * age;
  
  v_wts *= 0.5 * (v_cut - v0);
  v *= 0.5 * (v_cut - v0);
  v += 0.5 * (v_cut + v0);
  
  arma::rowvec L = arma::log(v / v0);
  arma::rowvec v_inv = 1 / v;
  
  arma::colvec e_t;
  arma::rowvec e_v;
  arma::mat e_tv;
  
  arma::mat d;
  arma::mat ix;
  
  arma::colvec x_t = arma::exp((par[1] - par[0]) * t);
  
  double C = par[0] * par[1] * par[2] * std::pow(par[1] - par[0], par[2]) / 
    tgamma(1 / par[3]) / std::pow(par[3] * par[4], 1 / par[3]) ;
  
  e_v = arma::exp(par[7] + par[6] * d0 * arma::exp(L / 3)) % 
    arma::pow(L,  -1 / par[3] - 1) % v_inv;
  
  e_t  = arma::log(age - t) / par[3]; 
  e_t += par[1] * par[2] * t; 
  
  e_tv = ((age - t) *  (1 / L)) * (1 / (par[4] * par[3]));
  e_tv += (age - t) * (par[5] * (v - v0) / L);
  e_tv -= ((age - t) *  (1 / L)) % (1 / (par[4] * par[3]) + par[5] * (v - v0));
  
  if (n_scr > 0) {
    for (int i=0; i<n_scr; i++) {
      d = d0 * arma::exp(((scr[i] - age)/ (age - t) + 1) * (L / 3));
      ix = arma::conv_to<arma::mat>::from(d > d0);
      e_tv -= arma::log(1 + arma::exp(par[7] + par[6] * d)) % ix;
    }
  }
  
  out = C * arma::sum(arma::sum((t_wts * v_wts) % ((arma::exp(e_t) % (1 - x_t) %
    pow(par[1] * x_t - par[0], -par[2] - 1)) * e_v) % arma::exp(e_tv), 0)) ;
  
  return out;
  
}

// [[Rcpp::export]]
double f_SymptSizeInterval(double age, double v_cut, arma::rowvec scr, 
                           arma::vec par, arma::colvec t, arma::colvec t_wts, 
                           arma::rowvec v, arma::rowvec v_wts, 
                           double d0 = 0.5, double v0 = 0.06544985)
{
  double out;
  int n_scr = scr.n_elem;
  
  t_wts *= 0.5 * age;
  t += 1;
  t *= 0.5 * age;
  
  v_wts *= 0.5 * (v_cut - v0);
  v *= 0.5 * (v_cut - v0);
  v += 0.5 * (v_cut + v0);
  
  
  arma::rowvec L = arma::log(v / v0);
  arma::rowvec v_inv = 1 / v;
  
  arma::colvec e_t;
  arma::rowvec e_v;
  arma::mat e_tv;
  
  arma::mat d;
  arma::mat ix;
  
  arma::colvec x_t = arma::exp((par[1] - par[0]) * t);
  
  double C = par[0] * par[1] * par[2] * std::pow(par[1] - par[0], par[2]) / 
    tgamma(1 / par[3]) / std::pow(par[3] * par[4], 1 / par[3]) ;
  
  e_v = par[5] * arma::pow(L,  -1 / par[3] - 1);
  
  e_t  = arma::log(age - t) / par[3]; 
  e_t += par[1] * par[2] * t; 
  
  e_tv = ((age - t) *  (1 / L)) * (1 / (par[4] * par[3]));
  e_tv += (age - t) * (par[5] * (v - v0) / L);
  e_tv -= ((age - t) *  (1 / L)) % (1 / (par[4] * par[3]) + par[5] * (v - v0));
  
  if (n_scr > 0) {
    for (int i=0; i<n_scr; i++) {
      d = d0 * arma::exp(((scr[i] - age)/ (age - t) + 1) * (L / 3));
      ix = arma::conv_to<arma::mat>::from(d > d0);
      e_tv -= arma::log(1 + arma::exp(par[7] + par[6] * d)) % ix;
    }
  }
  
  out = C * arma::sum(arma::sum((t_wts * v_wts) % ((arma::exp(e_t) % (1 - x_t) %
    pow(par[1] * x_t - par[0], -par[2] - 1)) * e_v) % arma::exp(e_tv), 0)) ;
  
  return out;
  
}

// [[Rcpp::export]]
double f_CaseIntegral4(double age, double v, arma::rowvec scr, 
                      arma::vec par, arma::rowvec scr_par, arma::colvec t, arma::colvec t_wts, 
                      double d0 = 0.5, double v0 = 0.06544985)
{
  double out;
  
  t_wts *= 0.5 * age;
  t += 1;
  t *= 0.5 * age;
  arma::colvec ones(t.n_elem, arma::fill::ones);
  
  arma::colvec e;
  arma::mat d;
  arma::mat ix;
  
  
  double L = std::log(v / v0);
  double C = par[0] * par[1] * par[2] * std::pow(par[1] - par[0], par[2]) / 
    tgamma(1 / par[3]) / std::pow(par[3] * par[4] * L, 1 / par[3]) / L;
  arma::colvec x_t = arma::exp((par[1] - par[0]) * t);
  
  d = d0 * arma::exp((((1 / (age - t)) * (scr - age)) + 1) * (L / 3));
  ix = arma::conv_to<arma::mat>::from(d > d0);
  
  e  = arma::log(age - t) / par[3]; 
  e += par[1] * par[2] * t; 
  e -= ((age - t) / L) * (1 / (par[4] * par[3]) + par[5] * (v - v0));
  e -= arma::sum(arma::log(1 + arma::exp(ones * scr_par + par[6] * d)) % ix , 1);
  
  out = C * arma::sum(t_wts % arma::exp(e) % (1 - x_t) %
    pow(par[1] * x_t - par[0], -par[2] - 1));
  
  return out;
  
}

// [[Rcpp::export]]
double f_CensIntegral4(double age, arma::rowvec scr, arma::vec par, arma::rowvec scr_par,
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
      d = d0 * arma::exp((scr[i] - t) * (r_inv / 3));
      ix = arma::conv_to<arma::mat>::from(d > d0);
      e_tr -= arma::log(1 + arma::exp(scr_par[i] + par[6] * d)) % ix;
    }
  }
  
  out = C * arma::sum(arma::sum((t_wts * r_wts) %
    (arma::exp(e_t) * arma::exp(e_r)) % arma::exp(e_tr), 0));
  
  return out;
  
}

// [[Rcpp::export]]
double f_ScreenCase4(double age, double v, arma::rowvec scr, arma::vec par, arma::rowvec scr_par,
                    arma::colvec t_pts, arma::colvec t_wts, 
                    double d0 = 0.5, double v0 = 0.06544985) 
{
  double out = std::exp(scr_par[0] + par[6] * d0 * std::exp(std::log(v / v0) / 3)) / v;
  out *= f_CaseIntegral4(age, v, scr, par, scr_par, t_pts, t_wts, d0, v0);
  return out;
}

// [[Rcpp::export]]
double f_SymptCase4(double age, double v, arma::rowvec scr, arma::vec par, arma::rowvec scr_par,
                   arma::colvec t_pts, arma::colvec t_wts, 
                   double d0 = 0.5, double v0 = 0.06544985) 
{
  double out = par[5];
  out *= f_CaseIntegral4(age, v, scr, par, scr_par, t_pts, t_wts, d0, v0);
  return out;
}

// [[Rcpp::export]]
double f_CensCase4(double age, arma::rowvec scr, arma::vec par, arma::rowvec scr_par,
                  arma::colvec t_pts, arma::colvec t_wts, 
                  arma::rowvec r_pts, arma::rowvec r_wts,
                  double d0 = 0.5, double v0 = 0.06544985) 
{
  double out = std::pow((par[1] - par[0]), par[2]) * 
    std::exp(par[1] * par[2] * age) /
      pow(par[1] * std::exp((par[1] - par[0]) * age) - par[0], par[2]);
  out += f_CensIntegral4(age, scr, par, scr_par, t_pts, t_wts, r_pts, r_wts, d0, v0);
  return out;
}

// [[Rcpp::export]]
double IndL4(int is_case, int is_scr, double age, double v, 
             arma::rowvec scr, double entry, arma::rowvec e_scr,
             arma::vec par, arma::rowvec scr_par, arma::rowvec e_scr_par,
             arma::colvec t_pts, arma::colvec t_wts,
             arma::rowvec r_pts, arma::rowvec r_wts,
             double d0 = 0.5, double v0 = 0.06544985){
  double out;
  if(is_case == 0) {
    out = f_CensCase4(age, scr, par, scr_par, t_pts, t_wts, r_pts, r_wts, d0, v0);
  } else if (is_scr == 1) {
    out = f_ScreenCase4(age, v, scr, par, scr_par, t_pts, t_wts, d0, v0);
  } else {
    out = f_SymptCase4(age, v, scr, par, scr_par, t_pts, t_wts, d0, v0);
  }
  
  out /= f_CensCase4(entry, e_scr, par, e_scr_par, t_pts, t_wts, r_pts, r_wts, d0, v0);
  return out;
}

// // [[Rcpp::export]]
// arma::mat f_SizeBounded(double x, double v, double c, double s, 
//              arma::rowvec scr, arma::vec par, 
//              arma::colvec gk_pts, arma::colvec gk_wts,
//              arma::rowvec gl_pts, arma::rowvec gl_wts,
//              double d0 = 0.5, double v0 = 0.06544985){
//   arma::mat out;
//   arma::mat t;
//   arma::mat t_w;
//   arma::mat r;
//   arma::mat r_w;
//   
//   arma::rowvec rowones(gl_pts.n_elem, arma::fill::ones);
//   arma::colvec colones(gk_pts.n_elem, arma::fill::ones);
//   arma::mat zeros(gk_pts.n_elem, gl_pts.n_elem, arma::fill::zeros);
//   arma::mat m;
//   arma::mat M;
//   
//   arma::mat x_t;
//   arma::mat r_inv;
//   arma::mat lnf_t;
//   arma::mat lnf_r;
//   arma::mat lns;
//   arma::mat lnu;
//   arma::mat d;
//   arma::mat ix;
//   
//   int n_scr = scr.n_elem;
//   double a  = 1 / par[3]; 
//   double b  = a / par[4];
//   double C  = -par[0] * par[1] * par[2] * std::pow(par[1] - par[0], par[2]) * 
//                 std::pow(b, a) / tgamma(a);
//   
//   
//   
//   if (s < 255) {
//   // define for screen and joint, t in (0, s), r in ((s-t)/L, inf)
//   t   = ((gk_pts + 1) * (s * 0.5)) * rowones;
//   t_w = gk_wts * (0.5 * s) * rowones;
//   r   = (colones * (gl_pts + 1)) * (0.5 * std::log(v / v0)) / (s - t);
//   r_w = (colones * gl_wts) * (0.5 * std::log(v / v0)) / (s - t);
//   
//   
//   x_t = arma::exp((par[1] - par[0]) * t);
//   r_inv = 1 / r;
//   
//   M = max(zeros + c, t);
//   m = min(zeros + x, t + r_inv * std::log(v / v0));
//   
//   lnf_t  = (par[1] * par[2]) * t;
//   lnf_t += arma::log(x_t - 1); 
//   lnf_t -= (par[2] + 1) * arma::log(par[1] * x_t - par[0]);
//   
//   lnf_r  = (-a - 1) * arma::log(r) - (r_inv * b); 
//   
//   // screen
//   d   = d0 * arma::exp((s - t) % (r / 3));
//   ix  = arma::conv_to<arma::mat>::from(d > d0);
//   lns = arma::log((1-1/(1 + arma::exp(par[7] + par[6] * d))) % ix);
//   if (n_scr > 0) {
//     for (int i=0; i<n_scr; i++) {
//       d    = d0 * arma::exp((scr[i] - t) % (r / 3));
//       ix   = arma::conv_to<arma::mat>::from(d > d0);
//       lns -= arma::log(1 + arma::exp(par[7] + par[6] * d)) % ix;
//     }
//   }
//   
//   // symptomatic (joint)
//   lnu = (-par[5] * v0) * r_inv % (arma::exp((m-t) % r) - arma::exp((M-t) % r));
//   
//   // put together
//   // out = arma::sum(arma::sum((arma::exp(lns) - (1-arma::exp(lnu)) % arma::exp(lns)) % 
//   //   arma::exp(lnf_t + lnf_r) % r_w % t_w, 0));
//   // out = arma::exp(lns) % arma::exp(lnf_t + lnf_r);
//   // } //else {
//   //   out = 0;
//   // }
//   // // redefine for just symptomatic, t in (0, x), r in ((x-t)/L, inf)
//   // t   = ((gk_pts + 1) * (x * 0.5)) * rowones;
//   // t_w = (gk_wts * (0.5 * x)) * rowones;
//   // r   = (colones * gl_pts) + max(c - t, zeros) / std::log(v / v0);
//   //   
//   //   r   = (colones * (gl_pts + 1)) * (0.5 * std::log(v / v0)) / (s - t);
//   //   r_w = (colones * gl_wts) * (0.5 * std::log(v / v0)) / (s - t);
//   //   
//   // 
//   // M = max(zeros + c, t);
//   // m = min(zeros + x, t + r * std::log(v / v0));
//   // 
//   // x_t = arma::exp((par[1] - par[0]) * t);
//   // r_inv = 1 / r;
//   // 
//   // lnf_t  = (par[1] * par[2]) * t;
//   // lnf_t += arma::log(x_t - 1); 
//   // lnf_t -= (par[2] + 1) * arma::log(par[1] * x_t - par[0]);
//   // 
//   // lnf_r  = (a - 1) * arma::log(r) - (r * (b - 1)); // +1r for gauss-laguerre
//   // 
//   // // just symptomatic
//   // lnu = (-par[5] * v0) * r % (arma::exp((m - t) % r_inv) - 
//   //         arma::exp((M - t) % r_inv));
//   // 
//   // 
//   // // Put together rest
//   // out += arma::sum(arma::sum((1-arma::exp(lnu)) % arma::exp(lnf_t + lnf_r) % 
//   //   r_w % t_w, 0));
//   // 
//   // don't forget the constant!
//   return C * out;
//   // return (arma::exp(lns) - (1-arma::exp(lnu)) % arma::exp(lns)) % 
//   //   arma::exp(lnf_t + lnf_r);
// }

// [[Rcpp::export]]
double f_SizeBounded_v1(double x, double v, double c, double s,
                     arma::rowvec scr, arma::vec par,
                     arma::colvec gk_pts, arma::colvec gk_wts,
                     arma::rowvec gl_pts, arma::rowvec gl_wts,
                     double d0 = 0.5, double v0 = 0.06544985){
  double out;
  arma::mat t;
  arma::mat t_w;
  arma::mat r;
  arma::mat r_w;

  arma::rowvec rowones(gl_pts.n_elem, arma::fill::ones);
  arma::colvec colones(gk_pts.n_elem, arma::fill::ones);
  arma::mat zeros(gk_pts.n_elem, gl_pts.n_elem, arma::fill::zeros);
  arma::mat m;
  arma::mat M;

  arma::mat x_t;
  arma::mat r_inv;
  arma::mat lnf_t;
  arma::mat lnf_r;
  arma::mat lns;
  arma::mat lnu;
  arma::mat d;
  arma::mat ix;

  int n_scr = scr.n_elem;
  double a  = 1 / par[3];
  double b  = a / par[4];
  double C  = -par[0] * par[1] * par[2] * std::pow(par[1] - par[0], par[2]) *
    std::pow(b, a) / tgamma(a);

  
  r   = colones * gl_pts; // + (s - t) / std::log(v / v0);
  r_w = colones * gl_wts;
  
  if (s < 255) {
    // define for screen and joint, r in (0, inf), t in (max(c-rlog(v/v0),0), s)

    t   = (gk_pts * (rowones * 0.5)) % (s - arma::max(zeros, c  - r * std::log(v / v0)));
    t  += (s + arma::max(zeros, c  - r * std::log(v / v0))) * 0.5;
    t_w = (gk_wts * (rowones * 0.5)) % (s - arma::max(zeros, c  - r * std::log(v / v0)));
    
    M = arma::max(zeros + c, t);
    m = arma::min(zeros + x, t + r * std::log(v / v0));

    x_t = arma::exp((par[1] - par[0]) * t);
    r_inv = 1 / r;

    lnf_t  = (par[1] * par[2]) * t;
    lnf_t += arma::log(x_t - 1);
    lnf_t -= (par[2] + 1) * arma::log(par[1] * x_t - par[0]);

    lnf_r  = (a - 1) * arma::log(r) - (r * (b - 1)); // +1r for gauss-laguerre

    // screen
    d   = d0 * arma::exp((s - t) % (r_inv / 3));
    ix  = arma::conv_to<arma::mat>::from(d > d0);
    lns = arma::log((1-1/(1 + arma::exp(par[7] + par[6] * d))) % ix);
    if (n_scr > 0) {
      for (int i=0; i<n_scr; i++) {
        d    = d0 * arma::exp((scr[i] - t) % (r_inv / 3));
        ix   = arma::conv_to<arma::mat>::from(d > d0);
        lns -= arma::log(1 + arma::exp(par[7] + par[6] * d)) % ix;
      }
    }

    // symptomatic (joint)
    lnu = (-par[5]*v0) * r % (arma::exp((m-t) % r_inv) - arma::exp((M-t) % r_inv));

    // put together
    out = arma::sum(arma::sum(arma::exp(lnu) % arma::exp(lns) %
      arma::exp(lnf_t + lnf_r) % r_w % t_w, 0));
  } else {
    out = 0;
  }
  // redefine for just symptomatic, r in (0, inf), t in (c-rlog(v/v0), x)
  t   = (gk_pts * (rowones * 0.5)) % (x - arma::max(zeros, c  - r * std::log(v / v0)));
  t  += (x + arma::max(zeros, c  - r * std::log(v / v0))) * 0.5;
  t_w = (gk_wts * (rowones * 0.5)) % (x - arma::max(zeros, c  - r * std::log(v / v0)));

  M = max(zeros + c, t);
  m = min(zeros + x, t + r * std::log(v / v0));

  x_t = arma::exp((par[1] - par[0]) * t);
  r_inv = 1 / r;

  lnf_t  = (par[1] * par[2]) * t;
  lnf_t += arma::log(x_t - 1);
  lnf_t -= (par[2] + 1) * arma::log(par[1] * x_t - par[0]);

  lnf_r  = (a - 1) * arma::log(r) - (r * (b - 1)); // +1r for gauss-laguerre

  // just symptomatic
  lnu = (-par[5] * v0) * r % (arma::exp((m - t) % r_inv) -
    arma::exp((M - t) % r_inv));


  // Put together rest
  out += arma::sum(arma::sum((1-arma::exp(lnu)) % arma::exp(lnf_t + lnf_r) %
    r_w % t_w, 0));

  // don't forget the constant!
  return C * out;
  // return C * (1 - arma::exp(lnu)) % arma::exp(lnf_t + lnf_r);
  // return t;
}

// [[Rcpp::export]]
double f_SizeBounded_v2(double x, double v, double c, double s,
                     arma::rowvec scr, arma::vec par,
                     arma::colvec gk_pts, arma::colvec gk_wts,
                     arma::rowvec gl_pts, arma::rowvec gl_wts, double k = 5,
                     double d0 = 0.5, double v0 = 0.06544985){
  double out;
  arma::mat t;
  arma::mat t_w;
  arma::mat r;
  arma::mat r_w;
  
  arma::rowvec rowones(gl_pts.n_elem, arma::fill::ones);
  arma::colvec colones(gk_pts.n_elem, arma::fill::ones);
  arma::mat zeros(gk_pts.n_elem, gl_pts.n_elem, arma::fill::zeros);
  arma::mat m;
  arma::mat M;
  
  arma::mat x_t;
  arma::mat r_inv;
  arma::mat lnf_t;
  arma::mat lnf_r;
  arma::mat lns;
  arma::mat lnu;
  arma::mat d;
  arma::mat ix;
  
  int n_scr = scr.n_elem;
  double a  = 1 / par[3];
  double b  = a / par[4];
  double C  = -par[0] * par[1] * par[2] * std::pow(par[1] - par[0], par[2]) *
    std::pow(b, a) / tgamma(a);
  
  r_w = colones * gl_wts;
  
  if (s < 255) {
    // define for screen and joint, t in (0, s), r in ((s-t)/L, (s-t)/L+k)
    t   = ((gk_pts + 1) * (s * 0.5)) * rowones;
    t_w = (gk_wts * (0.5 * s)) * rowones;
    
    r   = colones * (gl_pts * (0.5 * k)); // + (s - t) / std::log(v / v0);
    r  += (s - t) / std::log(v / v0);
    r  += k * 0.5;
    r_w = colones * (gl_wts * (0.5 * k));
    
    M = max(zeros + c, t);
    m = min(zeros + x, t + r * std::log(v / v0));
    
    x_t = arma::exp((par[1] - par[0]) * t);
    r_inv = 1 / r;
    
    lnf_t  = (par[1] * par[2]) * t;
    lnf_t += arma::log(x_t - 1);
    lnf_t -= (par[2] + 1) * arma::log(par[1] * x_t - par[0]);
    
    lnf_r  = (a - 1) * arma::log(r) - (r * b); // NOT +1r for gauss-laguerre
    
    // screen
    d   = d0 * arma::exp((s - t) % (r_inv / 3));
    ix  = arma::conv_to<arma::mat>::from(d > d0);
    lns = arma::log((1 - 1 / (1 + arma::exp(par[7] + par[6] * d))) % ix);
    if (n_scr > 0) {
      for (int i=0; i<n_scr; i++) {
        d    = d0 * arma::exp((scr[i] - t) % (r_inv / 3));
        ix   = arma::conv_to<arma::mat>::from(d > d0);
        lns -= arma::log(1 + arma::exp(par[7] + par[6] * d)) % ix;
      }
    }
    
    // symptomatic (joint)
    lnu = (-par[5] * v0) * r % (arma::exp((m - t) % r_inv) - arma::exp((M - t) % r_inv));
    
    // put together
    out = arma::sum(arma::sum(arma::exp(lnu) % arma::exp(lns) %
      arma::exp(lnf_t + lnf_r) % r_w % t_w, 0));
  } else {
    out = 0;
  }
  // redefine for just symptomatic, t in (0, x), r in ((x-t)/L, inf)
  t   = ((gk_pts + 1) * (x * 0.5)) * rowones;
  t_w = (gk_wts * (0.5 * x)) * rowones;
  
  r   = colones * (gl_pts * (0.5 * k)); // + (s - t) / std::log(v / v0);
  r  += (x - t) / std::log(v / v0);
  r  += k * 0.5;
  r_w = colones * (gl_wts * (0.5 * k));
  
  M = max(zeros + c, t);
  m = min(zeros + x, t + r * std::log(v / v0));
  
  x_t = arma::exp((par[1] - par[0]) * t);
  r_inv = 1 / r;
  
  lnf_t  = (par[1] * par[2]) * t;
  lnf_t += arma::log(x_t - 1);
  lnf_t -= (par[2] + 1) * arma::log(par[1] * x_t - par[0]);
  
  lnf_r  = (a - 1) * arma::log(r) - (r * b); // +1r for gauss-laguerre
  
  // just symptomatic
  lnu = (-par[5] * v0) * r % (arma::exp((m - t) % r_inv) -
    arma::exp((M - t) % r_inv));
  
  
  // Put together rest
  out += arma::sum(arma::sum((1-arma::exp(lnu)) % arma::exp(lnf_t + lnf_r) %
    r_w % t_w, 0));
  
  // don't forget the constant!
  return C * out;
  // return C * (1 - arma::exp(lnu)) % arma::exp(lnf_t + lnf_r);
}
// 
// [[Rcpp::export]]
double f_SizeBounded(double x, double v, double c, double s,
                     arma::rowvec scr, arma::vec par,
                     arma::colvec gk_pts, arma::colvec gk_wts,
                     arma::rowvec gl_pts, arma::rowvec gl_wts,
                     double d0 = 0.5, double v0 = 0.06544985){
  double out;
  arma::mat t;
  arma::mat t_w;
  arma::mat r;
  arma::mat r_w;

  arma::rowvec rowones(gl_pts.n_elem, arma::fill::ones);
  arma::colvec colones(gk_pts.n_elem, arma::fill::ones);
  arma::mat zeros(gk_pts.n_elem, gl_pts.n_elem, arma::fill::zeros);
  arma::mat m;
  arma::mat M;

  arma::mat x_t;
  arma::mat r_inv;
  arma::mat lnf_t;
  arma::mat lnf_r;
  arma::mat lns;
  arma::mat lnu;
  arma::mat d;
  arma::mat ix;

  int n_scr = scr.n_elem;
  double a  = 1 / par[3];
  double b  = a / par[4];
  double C  = -par[0] * par[1] * par[2] * std::pow(par[1] - par[0], par[2]) *
    std::pow(b, a) / tgamma(a);

  r_w = colones * gl_wts;

  if (s < 255) {
    // define for screen and joint, t in (0, s), r in ((s-t)/L, inf)
    t   = ((gk_pts + 1) * (s * 0.5)) * rowones;
    t_w = gk_wts * (0.5 * s) * rowones;
    r   = (colones * gl_pts) + (s - t) / std::log(v / v0);

    M = max(zeros + c, t);
    m = min(zeros + x, t + r * std::log(v / v0));

    x_t = arma::exp((par[1] - par[0]) * t);
    r_inv = 1 / r;

    lnf_t  = (par[1] * par[2]) * t;
    lnf_t += arma::log(x_t - 1);
    lnf_t -= (par[2] + 1) * arma::log(par[1] * x_t - par[0]);

    lnf_r  = (a - 1) * arma::log(r) - (r * (b - 1)); // +1r for gauss-laguerre

    // screen
    d   = d0 * arma::exp((s - t) % (r_inv / 3));
    ix  = arma::conv_to<arma::mat>::from(d > d0);
    lns = arma::log((1-1/(1 + arma::exp(par[7] + par[6] * d))) % ix);
    if (n_scr > 0) {
      for (int i=0; i<n_scr; i++) {
        d    = d0 * arma::exp((scr[i] - t) % (r_inv / 3));
        ix   = arma::conv_to<arma::mat>::from(d > d0);
        lns -= arma::log(1 + arma::exp(par[7] + par[6] * d)) % ix;
      }
    }

    // symptomatic (joint)
    lnu = (-par[5]*v0) * r % (arma::exp((m-t) % r_inv) - arma::exp((M-t) % r_inv));

    // put together
    out = arma::sum(arma::sum((arma::exp(lns) - (1-arma::exp(lnu)) % arma::exp(lns)) %
      arma::exp(lnf_t + lnf_r) % r_w % t_w, 0));
  } else {
    out = 0;
  }
  // redefine for just symptomatic, t in (0, x), r in ((x-t)/L, inf)
  t   = ((gk_pts + 1) * (x * 0.5)) * rowones;
  t_w = (gk_wts * (0.5 * x)) * rowones;
  r   = (colones * gl_pts) + max(c - t, zeros) / std::log(v / v0);

  M = max(zeros + c, t);
  m = min(zeros + x, t + r * std::log(v / v0));

  x_t = arma::exp((par[1] - par[0]) * t);
  r_inv = 1 / r;

  lnf_t  = (par[1] * par[2]) * t;
  lnf_t += arma::log(x_t - 1);
  lnf_t -= (par[2] + 1) * arma::log(par[1] * x_t - par[0]);

  lnf_r  = (a - 1) * arma::log(r) - (r * (b - 1)); // +1r for gauss-laguerre

  // just symptomatic
  lnu = (-par[5] * v0) * r % (arma::exp((m - t) % r_inv) -
    arma::exp((M - t) % r_inv));


  // Put together rest
  out += arma::sum(arma::sum((1-arma::exp(lnu)) % arma::exp(lnf_t + lnf_r) %
    r_w % t_w, 0));

  // don't forget the constant!
  return C * out;
  // return C * (1 - arma::exp(lnu)) % arma::exp(lnf_t + lnf_r);
}



// [[Rcpp::export]]
arma::colvec CondOnset(double age, double v, arma::rowvec scr, 
                       arma::vec par, arma::colvec t, arma::colvec t_wts, 
                       double d0 = 0.5, double v0 = 0.06544985)
{
  arma::colvec out;
  
  t_wts *= 0.5 * age;
  t += 1;
  t *= 0.5 * age;
  
  arma::colvec e;
  arma::mat d;
  arma::mat ix;
  
  double L = std::log(v / v0);
  arma::colvec x_t = arma::exp((par[1] - par[0]) * t);
  
  d = d0 * arma::exp((((1 / (age - t)) * (scr - age)) + 1) * (L / 3));
  ix = arma::conv_to<arma::mat>::from(d > d0);
  
  e  = arma::log((age - t) / par[4]) / par[3]; 
  e += par[1] * par[2] * t; 
  e -= ((age - t) / L) * (1 / (par[4] * par[3]) + par[5] * (v - v0));
  e -= arma::sum(arma::log(1 + arma::exp(par[7] + par[6] * d)) % ix , 1);
  
  out = arma::exp(e) % (1 - x_t) %
    pow(par[1] * x_t - par[0], -par[2] - 1);
  out /= arma::sum(t_wts % out);
  
  return out;
  
}

// [[Rcpp::export]]
arma::colvec f_pMVK(arma::colvec t, double A, double B, double delta){
  arma::colvec out = B * t;
  out -= arma::log(B * arma::exp((B - A) * t) - A);
  out += std::log(B - A);
  out *= delta;
  return 1 - arma::exp(out);
}

// [[Rcpp::export]]
arma::colvec f_dMVK(arma::colvec t, double A, double B, double delta){
  arma::colvec out = B * t + std::log(B - A);
  out *= delta;
  out += log(exp((B - A) * t) - 1);
  out -= log(B * exp((B - A) * t) - A) * (delta + 1);
  
  return -A * B * delta * arma::exp(out);
}

// [[Rcpp::export]]
arma::colvec RiskOnsetCPF(arma::colvec t, double age, arma::rowvec scr, 
                          arma::vec par, 
                          arma::colvec t_pts, arma::colvec t_wts, 
                          arma::rowvec r_pts, arma::rowvec r_wts,
                          double d0 = 0.5, double v0 = 0.06544985) 
{
  arma::colvec out = 1 - f_pMVK(t, par[0], par[1], par[2]);
  out /= f_CensCase(age, scr, par, t_pts, t_wts, r_pts, r_wts, d0, v0);
  return 1 - out;
}

// [[Rcpp::export]]
arma::colvec RiskOnsetPDF(arma::colvec t, double age, arma::rowvec scr, 
                          arma::vec par, 
                          arma::colvec t_pts, arma::colvec t_wts, 
                          arma::rowvec r_pts, arma::rowvec r_wts,
                          double d0 = 0.5, double v0 = 0.06544985) 
{
  arma::colvec out = f_dMVK(t, par[0], par[1], par[2]);
  out /= f_CensCase(age, scr, par, t_pts, t_wts, r_pts, r_wts, d0, v0);
  return out;
}


// [[Rcpp::export]]
double RiskSympt(double x, double age, arma::rowvec scr, 
                 arma::vec par,
                 arma::colvec t_pts, arma::colvec t_wts,
                 arma::rowvec r_pts, arma::rowvec r_wts,
                 double d0 = 0.5, double v0 = 0.06544985) {
  double out;
  
  out = f_CensCase(age + x, scr, par, t_pts, t_wts, r_pts, r_wts, d0, v0);
  out /= f_CensCase(age, scr, par, t_pts, t_wts, r_pts, r_wts, d0, v0);
  return 1 - out;
}

// [[Rcpp::export]]
double f_SensScreen(double x, double age, arma::rowvec scr, 
                    arma::vec par,
                    arma::colvec t_pts, arma::colvec t_wts,
                    arma::rowvec r_pts, arma::rowvec r_wts,
                    double d0 = 0.5, double v0 = 0.06544985) {
  double out;
  arma::rowvec r = r_pts;
  arma::colvec t;
  int n_scr = scr.n_elem;
  arma::colvec ones(t_pts.n_elem, arma::fill::ones);
  arma::colvec e_t;
  arma::rowvec e_r;
  arma::mat e_tr;
  arma::mat d;
  arma::mat ix;
  
  arma::rowvec r_inv = 1 / r;
  
  double C = -par[0] * par[1] * par[2] * std::pow(par[1] - par[0], par[2]) / 
    tgamma(1 / par[3]) / std::pow(par[3] * par[4], 1 / par[3]);
  
  // First integral from 0 to age
  t_wts *= 0.5 * age;
  t = t_pts + 1;
  t *= 0.5 * age;
  
  arma::colvec x_t = arma::exp((par[1] - par[0]) * t);
  
  e_t  = (par[1] * par[2]) * t;
  e_t += arma::log(x_t - 1); 
  e_t -= (par[2] + 1) * arma::log(par[1] * x_t - par[0]);
  e_r  = (par[5] * v0 + 1) * r; // +1r for gauss-laguerre
  e_r += (1 / par[3] - 1) * arma::log(r); 
  e_r -= (r / (par[3] * par[4])); 
  e_tr = (ones * ((-par[5] * v0) * r)) % arma::exp((age - t) * r_inv);
  
  d = d0 * arma::exp((age + x - t) * (r_inv / 3));
  ix = arma::conv_to<arma::mat>::from(d > d0);
  e_tr += arma::log(1 - 1 / (1 + arma::exp(par[7] + par[6] * d))) % ix;
  if (n_scr > 0) {
    for (int i=0; i<n_scr; i++) {
      d = d0 * arma::exp((scr[i] - t) * (r_inv / 3));
      ix = arma::conv_to<arma::mat>::from(d > d0);
      e_tr -= arma::log(1 + arma::exp(par[7] + par[6] * d)) % ix;
    }
  }
  
  
  arma::rowvec int1 = arma::sum((t_wts * r_wts) %
    (arma::exp(e_t) * arma::exp(e_r)) % arma::exp(e_tr), 0);
  
  // Second integral from age to age + x  
  
  t_wts *= x / age;
  t = t_pts * (x / 2);
  t += (2 * age + x) / 2;
  
  x_t = arma::exp((par[1] - par[0]) * t);
  
  e_t  = (par[1] * par[2]) * t;
  e_t += arma::log(x_t - 1); 
  e_t -= (par[2] + 1) * arma::log(par[1] * x_t - par[0]);
  e_r  = r; // +1r for gauss-laguerre
  e_r += (1 / par[3] - 1) * arma::log(r); 
  e_r -= (r / (par[3] * par[4])); 
  
  d = d0 * arma::exp((age + x - t) * (r_inv / 3));
  ix = arma::conv_to<arma::mat>::from(d > d0);
  e_tr = arma::log(1 - 1 / (1 + arma::exp(par[7] + par[6] * d))) % ix;
  if (n_scr > 0) {
    for (int i=0; i<n_scr; i++) {
      d = d0 * arma::exp((scr[i] - t) * (r_inv / 3));
      ix = arma::conv_to<arma::mat>::from(d > d0);
      e_tr -= arma::log(1 + arma::exp(par[7] + par[6] * d)) % ix;
    }
  }
  
  arma::rowvec int2 = arma::sum((t_wts * r_wts) %
    (arma::exp(e_t) * arma::exp(e_r)) % arma::exp(e_tr), 0);
  
  // Put together
  
  out = C * arma::sum(int1 + int2);
  
  // out /= f_CensCase(age, scr, par, t_pts, t_wts * 2 / x, r_pts, r_wts);
  return out;
}


// [[Rcpp::export]]
double SensScreen(double x, double age, arma::rowvec scr, 
                  arma::vec par,
                  arma::colvec t_pts, arma::colvec t_wts,
                  arma::rowvec r_pts, arma::rowvec r_wts,
                  double d0 = 0.5, double v0 = 0.06544985) {
  double out;
  
  out = f_SensScreen(x, age, scr, par, t_pts, t_wts, r_pts, r_wts, d0, v0);
  out /= f_CensCase(age, scr, par, t_pts, t_wts, r_pts, r_wts, d0, v0);
  
  return out;
}




// [[Rcpp::export]]
double f_RiskScreen(double age, arma::rowvec scr, arma::vec par,
                    arma::colvec t, arma::colvec t_wts,
                    arma::rowvec r, arma::rowvec r_wts,
                    double d0 = 0.5, double v0 = 0.06544985) {
  
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
  
  
  double C = -par[0] * par[1] * par[2] * std::pow(par[1] - par[0], par[2]) / 
    tgamma(1 / par[3]) / std::pow(par[3], 1 / par[3]);
  
  e_t  = (par[1] * par[2]) * t - std::log(par[4]) / par[3];
  e_t += arma::log(x_t - 1); 
  e_t -= (par[2] + 1) * arma::log(par[1] * x_t - par[0]);
  e_r  = (par[5] * v0 + 1) * r; // +1r for gauss-laguerre
  e_r += (1 / par[3] - 1) * arma::log(r); 
  e_r -= (r / (par[3] * par[4])); 
  e_tr = (ones * ((-par[5] * v0) * r)) % arma::exp((age - t) * r_inv);
  
  // the positive screen
  d = d0 * arma::exp((age - t) * (r_inv / 3));
  ix = arma::conv_to<arma::mat>::from(d > d0);
  e_tr += (par[7] + par[6] * d) % ix;
  
  // the negative screens
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
double f_RiskScreen_odg(double age, arma::rowvec scr, arma::vec par,
                        arma::colvec t, arma::colvec t_wts,
                        arma::rowvec r, arma::rowvec r_wts,
                        double d0 = 0.5, double v0 = 0.06544985) {
  
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
  
  e_t  = (par[1] * par[2]) * t + arma::log(b) * a;
  e_t += arma::log(x_t - 1); 
  e_t -= (par[2] + 1) * arma::log(par[1] * x_t - par[0]);
  e_r  = (par[5] * v0 + 1) * r; // +1r for gauss-laguerre
  e_r += (a - 1) * arma::log(r); 
  //
  e_tr = -b * r; 
  e_tr += (ones * ((-par[5] * v0) * r)) % arma::exp((age - t) * r_inv);
  
  // the positive screen
  d = d0 * arma::exp((age - t) * (r_inv / 3));
  ix = arma::conv_to<arma::mat>::from(d > d0);
  e_tr += (par[7] + par[6] * d) % ix;
  
  // the negative screens
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
