#ifndef PY_MODEL_H
#define PY_MODEL_H 1

#include <vector>
#include "Python.h"
#include "py_multipole.h"
#include "py_power_spectrum.h"

//
// Base class `Model`
//
class Model {
 public:
  Model(const double k_min_, const double dk_, const int nbin_,
	const double boxsize);
  virtual ~Model();
  virtual void evaluate(const std::vector<double>& param,
			std::vector<double>& v_P0,
			std::vector<double>& v_P2,
			std::vector<double>& v_P4) const = 0;
  void nmodes(std::vector<double>& v_nmodes) const;
    
  const double k_min, dk, boxsize;
  const int nbin;

  // Wave vectors in 3D grid
  std::vector<DiscreteWaveVector>* modes;

  // Coefficient of discrete Legendre multipoles
  vector<double> coef;
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

// Kaiser model with damping (aka dispersion model)
class Kaiser : public Model {
 public:
  Kaiser(const double k_min_, const double dk_, const int nbin_,
	 const double boxsize,
	 PyObject* py_k, PyObject* py_P);
  virtual void evaluate(const std::vector<double>& param,
			std::vector<double>& v_P0,
			std::vector<double>& v_P2,
			std::vector<double>& v_P4) const;

};

// Scoccimarro model
class Scoccimarro : public Model {
 public:
  Scoccimarro(const double k_min_, const double dk_, const int nbin_,
	      const double boxsize,
	      PyObject* py_k,
	      PyObject* py_Pdd, PyObject* py_Pdt, PyObject* py_Ptt);
  virtual void evaluate(const std::vector<double>& param,
			std::vector<double>& v_P0,
			std::vector<double>& v_P2,
			std::vector<double>& v_P4) const;
};


// Taruay model
class Taruya : public Model {
 public:
  Taruya(const double k_min_, const double dk_, const int nbin_,
	 const double boxsize,
	 PyObject* py_k,
	 PyObject* py_Pdd, PyObject* py_Pdt, PyObject* py_Ptt,
	 PyObject* py_AB);
  virtual void evaluate(const std::vector<double>& param,
			std::vector<double>& v_P0,
			std::vector<double>& v_P2,
			std::vector<double>& v_P4) const;
};


//PyObject* py_model_alloc(PyObject* self, PyObject* args);
PyObject* py_model_get_nbin(PyObject* self, PyObject* args);
PyObject* py_model_nmodes(PyObject* self, PyObject* args);
PyObject* py_model_exp_moment(PyObject* self, PyObject* args);
PyObject* py_model_evaluate(PyObject* self, PyObject* args);
  
PyObject* py_model_kaiser_alloc(PyObject* self, PyObject* args);
PyObject* py_model_scoccimarro_alloc(PyObject* self, PyObject* args);
#endif
