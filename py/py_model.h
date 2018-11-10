#ifndef PY_MODEL_H
#define PY_MODEL_H 1

#include "Python.h"

PyObject* py_model_alloc(PyObject* self, PyObject* args);
PyObject* py_model_exp_moment(PyObject* self, PyObject* args);

PyObject* py_model_kaiser_alloc(PyObject* self, PyObject* args);
PyObject* py_model_kaiser_eval(PyObject* self, PyObject* args);
#endif
