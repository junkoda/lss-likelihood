#include <vector>
#include <cmath>
#include <cassert>

#include "py_model.h"
#include "py_util.h"
#include "py_error.h"
#include "py_assert.h"

using namespace std;

static const double sqrt_pi = sqrt(M_PI);

//
//
//
class Model {
 public:
  virtual ~Model();
  virtual void evaluate(const double param[],
			vector<double>& v_P0,
			vector<double>& v_P2,
			vector<double>& v_P4) const = 0;

};

class Kaiser : public Model {
 public:
  Kaiser(PyObject* py_k, PyObject* py_P);
  virtual void evaluate(const double param[],
			vector<double>& v_P0,
			vector<double>& v_P2,
			vector<double>& v_P4) const;
  vector<double> v_k, v_P;
};

//
// static function
//


static void py_model_free(PyObject *obj)
{
  // Delete the power spectrum object
  // Called automatically by Python
  Model* const model= (Model*) PyCapsule_GetPointer(obj, "_Model");
  assert(model);

  delete model;
}

//
// static inline
//
static inline double exp_moment0(const double a, const double erf_a)
{
  return 0.5*sqrt_pi*erf_a/a;
}

static inline double exp_moment2(const double a,
				 const double a2,
				 const double erf_a,
				 const double exp_a)
{
  if(a > 0.015)
    return 0.25*sqrt_pi*erf_a/(a*a2) - 0.5*exp_a/a2;

  return 1.0/3.0 - a2/5.0 + a2*a2/14.0;
}

static inline double exp_moment4(const double a,
				 const double a2,
				 const double a4,
				 const double erf_a,
				 const double exp_a)
{
  if(a > 0.07)
    return (3.0/8.0)*sqrt_pi*erf_a/(a4*a) - (0.5/a2 + 0.75/a4)*exp_a;

  return 1.0/5.0 - a2/7.0 + a4/18.0 - a4*a2/66.0;
}

static inline double exp_moment6(const double a,
				 const double a2,
				 const double a4,
				 const double a6,
				 const double erf_a,
				 const double exp_a)
{
  if(a > 0.25)
    return (15.0/16.0)*sqrt_pi*erf_a/(a6*a)
      - (0.5/a2 + 5.0/(4.0*a4) + 15.0/(8.0*a6))*exp_a;

  return 1.0/7.0 - a2/9.0 + a4/22.0 - a6/78.0 + a6*a2/360.0 - a6*a4/2040.0;
}

static inline double exp_moment8(const double a,
				 const double a2,
				 const double a4,
				 const double a6,
				 const double a8,
				 const double erf_a,
				 const double exp_a)
{
  if(a > 0.3)
    return (105.0/32.0)*sqrt_pi*erf_a/(a*a8)
           - (0.5/a2 + 7.0/(4.0*a4) + 35.0/(8.0*a6) + 105.0/(16.0*a8))*exp_a;

  return 1.0/9.0 - a2/11.0 + a4/26.0 - a6/90.0 + a8/408.0 - a6*a4/2280.0;

}

//
// Python interface
//

// !!! deprecated !!!
PyObject* py_model_alloc(PyObject* self, PyObject* args)
{
  return NULL;
}

PyObject* py_model_exp_moment(PyObject* self, PyObject* args)
{
  double a;
  if(!PyArg_ParseTuple(args, "d", &a))
    return NULL;

  double a2= a*a;
  double a4= a2*a2;
  double a6= a4*a2;
  double a8= a4*a4;
  double exp_a= exp(-a2);
  double erf_a= erf(a);

  // \int_0^1 mu^n exp(-a2 mu^2)
  double fac0= exp_moment0(a, erf_a);
  double fac2= exp_moment2(a, a2, erf_a, exp_a);
  double fac4= exp_moment4(a, a2, a4, erf_a, exp_a);
  double fac6= exp_moment6(a, a2, a4, a6, erf_a, exp_a);
  double fac8= exp_moment8(a, a2, a4, a6, a8, erf_a, exp_a);


  return Py_BuildValue("ddddd", fac0, fac2, fac4, fac6, fac8);
}


PyObject* py_model_kaiser_alloc(PyObject* self, PyObject* args)
{
  // _model_kaiser_alloc(k, P)
  PyObject *py_k, *py_P;
    
  if(!PyArg_ParseTuple(args, "OO", &py_k, &py_P))
    return NULL;

  try {
    return PyCapsule_New(new Kaiser(py_k, py_P), "_Model", py_model_free);
  }
  catch(TypeError) {
    return NULL;
  }

  return NULL;
}

PyObject* py_model_kaiser_eval(PyObject* self, PyObject* args)
{
  // _model_kaiser_alloc(_model, b, f, sigma, P0, P2, P4)
  //
  // Args:
  // _model: _Model pointer
  // b, f, sigma (double): parameters
  // P0, P2, P4 (array): model prediction
  //
  double b, f, sigma;
  PyObject *py_model, *py_P0, *py_P2, *py_P4;
    
  if(!PyArg_ParseTuple(args, "OdddOOO",
		       &py_model,
		       &b, &f, &sigma,
		       &py_P0, &py_P2, &py_P4))
    return NULL;

  Model* const model=
    (Model*) PyCapsule_GetPointer(py_model, "_Model");
  py_assert_ptr(model);

  Kaiser* const kaiser= dynamic_cast<Kaiser*>(model);
  
  const size_t n= kaiser->v_k.size();
  vector<double> v_P0(n), v_P2(n), v_P4(n);
  double params[]= {b, f, sigma};
  
  kaiser->evaluate(params, v_P0, v_P2, v_P4);
  
  try {
    py_util_vector_as_array("P0", v_P0, py_P0);
    py_util_vector_as_array("P2", v_P2, py_P2);
    py_util_vector_as_array("P4", v_P4, py_P4);
  }
  catch(TypeError) {
    return NULL;
  }

  Py_RETURN_NONE;
}

//
// Model implementation
//
Model::~Model()
{

}


//
// Kaiser
//
Kaiser::Kaiser(PyObject* py_k, PyObject* py_P)
{
  py_util_array_as_vector("k", py_k, v_k);
  py_util_array_as_vector("P", py_P, v_P, v_k.size());
}

void Kaiser::evaluate(const double params[],
		      vector<double>& v_P0,
		      vector<double>& v_P2,
		      vector<double>& v_P4) const
{
  const size_t n= v_P.size();
  assert(v_P0.size() == n);
  assert(v_P2.size() == n);
  assert(v_P4.size() == n);

  const double b= params[0];
  const double f= params[1];
  const double s= params[2];

  const double b2= b*b;
  const double bf2= 2.0*b*f;
  const double ff= f*f;

  for(size_t i=0; i<n; ++i) {
    double k= v_k[i];
    double a= k*s;
    double a2= a*a;
    double a4= a2*a2;
    double a6= a4*a2;
    double a8= a4*a4;
    double exp_a= exp(-a2);
    double erf_a= erf(a);

    // \int_0^1 mu^n exp(-a2 mu^2)
    double fac0= exp_moment0(a, erf_a);
    double fac2= exp_moment2(a, a2, erf_a, exp_a);
    double fac4= exp_moment4(a, a2, a4, erf_a, exp_a);
    double fac6= exp_moment6(a, a2, a4, a6, erf_a, exp_a);
    double fac8= exp_moment8(a, a2, a4, a6, a8, erf_a, exp_a);
    
    // (2l + 1) for multipole moments
    v_P0[i]= (b2*fac0 + bf2*fac2 + ff*fac4)*v_P[i];
    v_P2[i]= 2.5*(   b2*(3.0*fac2 - fac0)
		  + bf2*(3.0*fac4 - fac2)
		  +  ff*(3.0*fac6 - fac4))*v_P[i];
    v_P4[i]= (9.0/8.0)*(   b2*(35.0*fac4 - 30.0*fac2 + 3.0*fac0)
			- bf2*(35.0*fac6 - 30.0*fac4 + 3.0*fac2)
			+  ff*(35.0*fac8 - 30.0*fac6 + 3.0*fac4))*v_P[i];
  }

}



