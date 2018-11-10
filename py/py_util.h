#ifndef PY_UTIL_H
#define PY_UTIL_H 1

#include <vector>
#include "Python.h"

void py_util_array_as_vector(const char name[],
			     PyObject* py_obj,
			     std::vector<double>& v,
			     const Py_ssize_t len=0);

void py_util_vector_as_array(const char name[], std::vector<double>& v,
			     PyObject* py_obj);


#endif
