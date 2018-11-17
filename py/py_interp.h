#ifndef PY_INTERP_PY
#define PY_INTREP_PY 1

#include <vector>
#include <gsl/gsl_spline.h>
#include "Python.h"

class Interp {
 public:
  Interp(const std::vector<double>& v, const size_t nrow, const int ncol_);
  ~Interp();
  //double interp(const size_t icol, const double k);
  double operator()(const size_t icol, const double k);
 private:
  const size_t ncol;
  gsl_spline** spline;
  gsl_interp_accel** acc;
};

#endif
