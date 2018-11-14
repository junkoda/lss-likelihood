#include <vector>
#include <cmath>
#include <cassert>

#include "py_util.h"
#include "py_power_spectrum.h"

using namespace std;

PowerSpectrum::PowerSectrum(PyObject* py_k, PyObject* py_P)
{
  if(!PyArg_ParseTuple(args, "OO", &py_k, &py_P))
    return NULL;

  vector<double> v_k, v_P;
  py_util_array_as_vector("k", py_k, v_k, 0);

  const size_t n= v_k.size(); assert(n > 0);
  py_util_array_as_vector("P", py_P, v_P, n);

  acc= gsl_interp_accel_alloc();
  spline= gsl_spline_alloc(gsl_interp_cspline, n);
  gsl_spline_init(spline, v_k.data(), v_P.data, n);
}

PowerSpectrum::~PowerSpectrum()
{
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
}

double PowerSpectrum::P(const double k)
{
  return gsl_spline_eval(spline, k, acc);
}
