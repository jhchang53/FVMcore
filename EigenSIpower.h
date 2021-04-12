/*
 * 	EigenSIpower.h
 * 	shift and invert solver
 */
#include "Eigen.h"

class EigenSIpower : public Eigen {
public:
  EigenSIpower();
  ~EigenSIpower();
  void setProblem(Mat Amat, Mat Bmat, double sigma);
  void guess(Vec phi);
  double solve();
private:
  double sigma;
  Mat ASmat,Bmat;
};
