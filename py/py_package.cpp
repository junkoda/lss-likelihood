//
// Information for Python
//
#include <iostream>
#include "Python.h"

#include "py_model.h"
#include "py_data.h"
#include "py_likelihood.h"

using namespace std;

static PyMethodDef methods[] = {
  {"_model_exp_moment", py_model_exp_moment, METH_VARARGS,
   "_model_exp_moment(a)"},
  {"_model_get_nbin", py_model_get_nbin, METH_VARARGS,
   "_model_get_nbin(_model)"},
  {"_model_nmodes", py_model_nmodes, METH_VARARGS,
   "_model_nmodes(_model)"},
  {"_model_kaiser_alloc", py_model_kaiser_alloc, METH_VARARGS,
   "_model_kaiser_alloc(), returns _Model pointer"},
  {"_model_evaluate", py_model_evaluate, METH_VARARGS,
   "_model_evaluate(_model, b, f, sigma, P0, P2, P4)"},
  {"_model_scoccimarro_alloc", py_model_scoccimarro_alloc, METH_VARARGS,
   "_model_scoccimarro_alloc(), returns _Model pointer"},
  {"_data_alloc", py_data_alloc, METH_VARARGS,
   "_data_alloc(P0, P2, P4)"},
  {"_likelihood_alloc", py_likelihood_alloc, METH_VARARGS,
   "_data_alloc(_data, _model, cov_inv)"},
  {"_likelihood_chi2", py_likelihood_chi2,  METH_VARARGS,
   "_data_alloc(_data, _model, _likelihood, params)"},   

  {NULL, NULL, 0, NULL}
};


static struct PyModuleDef module = {
  PyModuleDef_HEAD_INIT,
  "_lss_likelihood", // name of this module
  "large-scale structure data-model chi2 computation", // Doc String
  -1,
  methods
};

PyMODINIT_FUNC
PyInit__lss_likelihood(void) {
  //py_spherical_bessel_module_init();
  
  return PyModule_Create(&module);
}
