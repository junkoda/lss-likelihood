//
// Information for Python
//
#include <iostream>
#include "Python.h"

#include "py_model.h"

using namespace std;

static PyMethodDef methods[] = {
  {"_model_alloc",
   py_model_alloc, METH_VARARGS,
   "temp"},
  
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
