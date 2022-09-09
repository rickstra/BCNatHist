#ifndef __LIKELIHOOD__
#define __LIKELIHOOD__

  double f_CaseIntegral(double age, double v, arma::rowvec scr, 
                      arma::vec par, arma::colvec t, arma::colvec t_wts, 
                      double d0, double v0);
  double f_CensIntegral(double age, arma::rowvec scr, arma::vec par,
                        arma::colvec t, arma::colvec t_wts,
                        arma::rowvec r, arma::rowvec r_wts,
                        double d0, double v0);
  double f_ScreenCase(double age, double v, arma::rowvec scr, arma::vec par, 
                      arma::colvec t_pts, arma::colvec t_wts, 
                      double d0, double v0);
  double f_SymptCase(double age, double v, arma::rowvec scr, arma::vec par, 
                     arma::colvec t_pts, arma::colvec t_wts, 
                     double d0, double v0);  
  double f_CensCase(double age, arma::rowvec scr, arma::vec par, 
                    arma::colvec t_pts, arma::colvec t_wts, 
                    arma::rowvec r_pts, arma::rowvec r_wts,
                    double d0, double v0);
  double IndL(int is_case, int is_scr, double age, double v, arma::rowvec scr,
              double entry, arma::rowvec e_scr, arma::vec par,
              arma::colvec t_pts, arma::colvec t_wts,
              arma::rowvec r_pts, arma::rowvec r_wts,
              double d0, double v0);
  
#endif // __LIKELIHOOD__