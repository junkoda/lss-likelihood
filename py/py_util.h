#ifndef PY_UTIL_H
#define PY_UTIL_H 1

#include <vector>
#include "Python.h"

void py_util_array_as_vector(const char name[],
			     PyObject* py_obj,
			     std::vector<double>& v,
			     const Py_ssize_t len=0);

size_t py_util_array2_as_vector(const char name[],
				PyObject* py_array,
				std::vector<double>& v,
				const Py_ssize_t len_expect=0,
				const Py_ssize_t ncol_expect=0);

void py_util_vector_as_array(const char name[], std::vector<double>& v,
			     PyObject* py_obj);

void py_util_sequence_as_vector(const char name[], PyObject* py_list,
				std::vector<double>& v);


#endif
