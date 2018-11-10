#ifndef PY_LIKELIHOOD_H
#define PY_LIKELIHOOD_H 1

#include "Python.h"

PyObject* py_likelihood_alloc(PyObject* self, PyObject* args);
PyObject* py_likelihood_chi2(PyObject* self, PyObject* args);

#endif
