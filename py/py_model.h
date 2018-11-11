#ifndef PY_MODEL_H
#define PY_MODEL_H 1

#include <vector>
#include "Python.h"


//
//
//
class Model {
 public:
  virtual ~Model();
  virtual void evaluate(const std::vector<double>& param,
			std::vector<double>& v_P0,
			std::vector<double>& v_P2,
			std::vector<double>& v_P4) const = 0;

};

class Kaiser : public Model {
 public:
  Kaiser(PyObject* py_k, PyObject* py_P);
  virtual void evaluate(const std::vector<double>& param,
			std::vector<double>& v_P0,
			std::vector<double>& v_P2,
			std::vector<double>& v_P4) const;

  // real-space power spectrum
  std::vector<double> v_k, v_P;
};

PyObject* py_model_alloc(PyObject* self, PyObject* args);
PyObject* py_model_exp_moment(PyObject* self, PyObject* args);

PyObject* py_model_kaiser_alloc(PyObject* self, PyObject* args);
PyObject* py_model_kaiser_evaluate(PyObject* self, PyObject* args);
#endif
