#include <iostream> // DEBUG
#include <vector>

#include "py_data.h"
#include "py_model.h"
#include "py_util.h"
#include "py_assert.h"
#include "py_error.h"
#include "py_likelihood.h"

using namespace std;

class Likelihood {
 public:
  Likelihood(Data const * const data_,
	     Model const * const model_);
  double chi2(const vector<double>& params) const;
  vector<double> cov_inv;
 private:
  Data const * const data;
  Model const * const model;
};

Likelihood::Likelihood(Data const * const data_,
		       Model const * const model_) :
  data(data_), model(model_)
{

}

double Likelihood::chi2(const vector<double>& params) const
{
  //
  // Likelihood evaluation
  //
  const size_t n= data->P0.size();
  assert(n == static_cast<size_t>(model->nbin));

  vector<double> model_P0(model->nbin),
    model_P2(model->nbin), model_P4(model->nbin);
  assert(cov_inv.size() == 9*n);

  model->evaluate(params, model_P0, model_P2, model_P4);

  double chi2= 0.0;

  for(size_t i=0; i<n; ++i) {
    double x0= data->P0[i] - model_P0[i];
    double x1= data->P2[i] - model_P2[i];
    double x2= data->P4[i] - model_P4[i];

    size_t i9 = 9*i;
    double c00 = cov_inv[i9    ];
    double c01 = cov_inv[i9 + 1];
    double c02 = cov_inv[i9 + 2];
    double c11 = cov_inv[i9 + 4];
    double c12 = cov_inv[i9 + 5];
    double c22 = cov_inv[i9 + 8];
    
    chi2 += c00*x0*x0 + 2.0*c01*x0*x1 + 2.0*c02*x0*x2
            + c11*x1*x1 + 2.0*c12*x1*x2 + c22*x2*x2;
  }
  
  return chi2;
}

static void py_likelihood_free(PyObject *obj)
{
  // Delete the data object
  // Called automatically by Python
  Likelihood* const li= (Likelihood*) PyCapsule_GetPointer(obj, "_Likelihood");
  assert(li);

  delete li;
}


PyObject* py_likelihood_alloc(PyObject* self, PyObject* args)
{
  // _likelihood_alloc(_data, _model, cov)
  // Args:
  //   _data: _Data pointer
  //   _model: _Model pointer
  //   _cov_inv: array of covariance matrix inverse
  
  PyObject *py_data, *py_model, *py_cov;
  if(!PyArg_ParseTuple(args, "OOO",
		       &py_data, &py_model, &py_cov))
    return NULL;

  Model const * const model=
    (Model const *) PyCapsule_GetPointer(py_model, "_Model");
  py_assert_ptr(model);

  Data const * const data=
    (Data const *) PyCapsule_GetPointer(py_data, "_Data");
  py_assert_ptr(data);

  const size_t nk= data->P0.size();

  Likelihood* const li= new Likelihood(data, model);


  // Fill Cov_inv
  Py_buffer buf;
  
  if(PyObject_GetBuffer(py_cov, &buf, PyBUF_FORMAT | PyBUF_FULL_RO) == -1)
      throw TypeError();

  char msg[128];
  if(buf.ndim != 3) {
    sprintf(msg, "Expected a 3-dimensional array for cov_inv: %d",
	    (int) buf.ndim);
    PyErr_SetString(PyExc_TypeError, msg);
    throw TypeError();
  }

  if(strcmp(buf.format, "d") != 0) {
    sprintf(msg, "Expected an array of double for cov_inv: %s", buf.format);
    PyErr_SetString(PyExc_TypeError, msg);
    throw TypeError();
  }

  if(nk != (size_t) buf.shape[0]) {
    sprintf(msg, "Expected length %d in k direction: %d",
	    (int) nk, (int) buf.shape[0]);
    PyErr_SetString(PyExc_TypeError, msg);
    throw TypeError();
  }

  if(buf.shape[1] != 3 || buf.shape[2] != 3) {
    sprintf(msg, "Expected length 3 in 1st and 2nd axis %d %d",
	    (int) buf.shape[1], (int) buf.shape[2]);
    PyErr_SetString(PyExc_TypeError, msg);
    throw TypeError();
  }


  const size_t stride_k= buf.strides[0];
  const size_t stride_i= buf.strides[1];
  const size_t stride_j= buf.strides[2];

  const size_t size= nk*3*3;
  li->cov_inv.reserve(size);

  double const * c_k= (double const*) buf.buf;
  for(size_t ik=0; ik<nk; ++ik) {
    double const * c_i= c_k;
    for(size_t i=0; i<3; ++i) {
      double const * c= c_i;
      for(size_t j=0; j<3; ++j) {
	li->cov_inv.push_back(*c);
	c = (double const *) ((char const *) c + stride_j);
      }
      c_i = (double const *) ((char const *) c_i + stride_i);      
    }
    c_k = (double const *) ((char const *) c_k + stride_k);
  }
  assert(li->cov_inv.size() == size);

  /*
  cerr << "DEBUG\n";
  for(size_t ik=0; ik<nk; ++ik) {
    for(size_t i=0; i<3; ++i) {
      for(size_t j=0; j<3; ++j) {
	cerr << li->cov_inv[9*ik + 3*i + j] << endl;
      }
    }
  }
  */

  
  return PyCapsule_New(li, "_Likelihood", py_likelihood_free);
}


PyObject* py_likelihood_chi2(PyObject* self, PyObject* args)
{
  // _likelihood_chi2(_likelihood, params)
  // Args:
  //   _likelihood: _Likelihood pointer
  //   params (list): sequence of model paramters
  
  PyObject *py_likelihood, *py_params;
  if(!PyArg_ParseTuple(args, "OO",
		       &py_likelihood, &py_params))
    return NULL;

  Likelihood* const li= (Likelihood*) PyCapsule_GetPointer(py_likelihood, "_Likelihood");
  py_assert_ptr(li);

  vector<double> params;
  
  try {
    py_util_sequence_as_vector("params", py_params, params);
  }
  catch(TypeError) {
    return NULL;
  }

  const double chi2= li->chi2(params);

  
  return Py_BuildValue("d", chi2);
}
