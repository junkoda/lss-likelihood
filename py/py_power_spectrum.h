#ifndef PY_POWER_SPECTRUM_H
#define PY_POWER_SPECTRUM_H 1

#include <gsl/gsl_spline.h>
#include "Python.h"

class PowerSpectrum {
 public:
  PowerSpectrum(PyObject* py_k, PyObject* py_P);
  ~PowerSpectrum();
  double P(const double k);
 private:
  gsl_spline* spline;
  gsl_interp_accel* acc;
};

#endif
