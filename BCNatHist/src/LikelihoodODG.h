#ifndef __LIKELIHOODODG__
#define __LIKELIHOODODG__

  double EvalCaseIntegralODG(double age, double v, arma::rowvec scr, 
                      arma::vec par, arma::colvec t, arma::colvec t_wts, 
                      double d0, double v0);
  double EvalCensIntegralODG(double age, arma::rowvec scr, arma::vec par,
                        arma::colvec t, arma::colvec t_wts,
                        arma::rowvec r, arma::rowvec r_wts,
                        double d0, double v0);
  double CalcScreenCaseODG(double age, double v, arma::rowvec scr, arma::vec par, 
                      arma::colvec t_pts, arma::colvec t_wts, 
                      double d0, double v0);
  double CalcSymptCaseODG(double age, double v, arma::rowvec scr, arma::vec par, 
                     arma::colvec t_pts, arma::colvec t_wts, 
                     double d0, double v0);  
  double CalcCensCaseODG(double age, arma::rowvec scr, arma::vec par, 
                    arma::colvec t_pts, arma::colvec t_wts, 
                    arma::rowvec r_pts, arma::rowvec r_wts,
                    double d0, double v0);
  double CalcIndividualLikelihoodODG(int is_case, int is_scr, double age, double v, arma::rowvec scr,
              double entry, arma::rowvec e_scr, arma::vec par,
              arma::colvec t_pts, arma::colvec t_wts,
              arma::rowvec r_pts, arma::rowvec r_wts,
              double d0, double v0);
  
#endif // __LIKELIHOODODG__