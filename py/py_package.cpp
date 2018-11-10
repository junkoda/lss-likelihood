//
// Information for Python
//
#include <iostream>
#include "Python.h"

#include "py_model.h"

using namespace std;

static PyMethodDef methods[] = {
  {"_model_exp_moment", py_model_exp_moment, METH_VARARGS,
   "_model_kaiser_eval(a)"},
  {"_model_kaiser_alloc", py_model_kaiser_alloc, METH_VARARGS,
   "_model_kaiser_alloc(), returns _Model pointer"},
  {"_model_kaiser_eval", py_model_kaiser_eval, METH_VARARGS,
   "_model_kaiser_eval(_model, b, f, sigma, P0, P2, P4)"},
  
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
