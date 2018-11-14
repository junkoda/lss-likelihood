#ifndef PY_MODEL_H
#define PY_MODEL_H 1

#include <vector>
#include "Python.h"
#include "py_multipole.h"
#include "py_power_spectrum.h"

//
//
//
class Model {
 public:
  Model(const double k_min_, const double k_max_, const double dk_,
	const double boxsize);
  virtual ~Model();
  virtual void evaluate(const std::vector<double>& param,
			std::vector<double>& v_P0,
			std::vector<double>& v_P2,
			std::vector<double>& v_P4) const = 0;
  void nmodes(std::vector<double>& v_nmodes) const;
    
  const double k_min, k_max, dk, boxsize;
  const int nbin;
  std::vector<DiscreteWaveVector>* modes;
};

/*
class Kaiser : public Model {
 public:
  Kaiser(const double k_min_, const double k_max_, const double dk_,
	 PyObject* py_k, PyObject* py_P);
  virtual void evaluate(const std::vector<double>& param,
			std::vector<double>& v_P0,
			std::vector<double>& v_P2,
			std::vector<double>& v_P4) const;

  // real-space power spectrum
  std::vector<double> v_k, v_P;
};
*/

class Kaiser : public Model {
 public:
  Kaiser(const double k_min_, const double k_max_, const double dk_,
	 const double boxsize,
	 PyObject* py_k, PyObject* py_P);
  virtual void evaluate(const std::vector<double>& param,
			std::vector<double>& v_P0,
			std::vector<double>& v_P2,
			std::vector<double>& v_P4) const;

  // real-space power spectrum
  //std::vector<double> v_k, v_P;

  // coefficient of discrete Legendre multipoles
  vector<double> coef;
};


//PyObject* py_model_alloc(PyObject* self, PyObject* args);
PyObject* py_model_get_nbin(PyObject* self, PyObject* args);
PyObject* py_model_nmodes(PyObject* self, PyObject* args);
PyObject* py_model_exp_moment(PyObject* self, PyObject* args);

PyObject* py_model_kaiser_alloc(PyObject* self, PyObject* args);
PyObject* py_model_kaiser_evaluate(PyObject* self, PyObject* args);
#endif
