#ifndef PY_BETTER_MODEL_H
#define PY_BETTER_MODEL_H 1

#include <vector>
#include "Python.h"
#include "py_multipole.h"
#include "py_power_spectrum.h"
#include "py_model.h"

class Model1 : public Model {
 public:
  Model1(const double k_min_, const double dk_, const int nbin_,
	 const double boxsize,
	 PyObject* py_k,
	 PyObject* py_Pdd, PyObject* py_Pdt, PyObject* py_Ptt,
	 PyObject* py_AB);
  virtual void evaluate(const std::vector<double>& param,
			std::vector<double>& v_P0,
			std::vector<double>& v_P2,
			std::vector<double>& v_P4) const;
};

PyObject* py_better_model_alloc(PyObject* self, PyObject* args);

#endif
