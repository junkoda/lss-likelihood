#ifndef PY_DATA_H
#define PY_DATA_H 1

#include <vector>
#include "Python.h"

class Data {
 public:
  std::vector<double> P0, P2, P4;
};

PyObject* py_data_alloc(PyObject* self, PyObject* args);

#endif
