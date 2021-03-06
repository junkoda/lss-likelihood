//#include <iostream> // DEBUG
#include <vector>
#include <cmath>
#include <cassert>


#include "py_util.h"
#include "py_error.h"
#include "py_assert.h"
#include "py_multipole.h"
#include "py_interp.h"
#include "py_model.h"

using namespace std;

static const double sqrt_pi = sqrt(M_PI);


//
// static function
//

void py_model_free(PyObject *obj)
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

PyObject* py_model_get_nbin(PyObject* self, PyObject* args)
{
  PyObject *py_model;

  if(!PyArg_ParseTuple(args, "O",
		       &py_model))
    return NULL;

  Model* const model=
    (Model*) PyCapsule_GetPointer(py_model, "_Model");
  py_assert_ptr(model);

  return Py_BuildValue("i", model->nbin);
}

PyObject* py_model_nmodes(PyObject* self, PyObject* args)
{
  // _model_nmodes(_model)
  //
  // Args:
  // _model: _Model pointer
  //
  PyObject *py_model, *py_nmodes;

  if(!PyArg_ParseTuple(args, "OO",
		       &py_model,
		       &py_nmodes))
    return NULL;

  Model* const model=
    (Model*) PyCapsule_GetPointer(py_model, "_Model");
  py_assert_ptr(model);

  const size_t n= model->nbin;
  
  vector<double> v_nmodes(n);
  model->nmodes(v_nmodes);
  
  try {
    py_util_vector_as_array("nmodes", v_nmodes, py_nmodes);
  }
  catch(TypeError) {
    return NULL;
  }

  Py_RETURN_NONE;
}

// deprecate this??
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

PyObject* py_model_evaluate(PyObject* self, PyObject* args)
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

  const size_t n= model->nbin;
  
  vector<double> v_P0(n), v_P2(n), v_P4(n);
  vector<double> params(3, 0.0);
  params[0]= b; params[1]= f; params[2]= sigma;

  model->evaluate(params, v_P0, v_P2, v_P4);
  
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

// Kaiser
PyObject* py_model_kaiser_alloc(PyObject* self, PyObject* args)
{
  // _model_kaiser_alloc(k_min, dk, nbin, boxsize, k, P)
  double k_min, dk, boxsize;
  int nbin;
  PyObject *py_k, *py_P;
    
  if(!PyArg_ParseTuple(args, "ddidOO",
		       &k_min, &dk, &nbin, &boxsize, &py_k, &py_P))
    return NULL;

  try {
    return PyCapsule_New(new Kaiser(k_min, dk, nbin, boxsize, py_k, py_P),
			 "_Model", py_model_free);
  }
  catch(TypeError) {
    return NULL;
  }

  return NULL;
}


// Scoccimarro
PyObject* py_model_scoccimarro_alloc(PyObject* self, PyObject* args)
{
  // _model_scoccimarro_alloc(k_min, dk, nbin, boxsize, k, Pdd, Pdt, Ptt)
  double k_min, dk, boxsize;
  int nbin;
  PyObject *py_k, *py_Pdd, *py_Pdt,*py_Ptt;
    
  if(!PyArg_ParseTuple(args, "ddidOOOO",
		       &k_min, &dk, &nbin, &boxsize,
		       &py_k, &py_Pdd, &py_Pdt, &py_Ptt))
    return NULL;

  try {
    return PyCapsule_New(new Scoccimarro(k_min, dk, nbin, boxsize, py_k,
					 py_Pdd, py_Pdt, py_Ptt),
			 "_Model", py_model_free);
  }
  catch(TypeError) {
    return NULL;
  }

  return NULL;
}


// Taruya
PyObject* py_model_taruya_alloc(PyObject* self, PyObject* args)
{
  // _model_taruya_alloc(k_min, dk, nbin, boxsize, k, Pdd, Pdt, Ptt, AB)
  double k_min, dk, boxsize;
  int nbin;
  PyObject *py_k, *py_Pdd, *py_Pdt,*py_Ptt, *py_AB;
    
  if(!PyArg_ParseTuple(args, "ddidOOOOO",
		       &k_min, &dk, &nbin, &boxsize,
		       &py_k, &py_Pdd, &py_Pdt, &py_Ptt, &py_AB))
    return NULL;

  try {
    return PyCapsule_New(new Taruya(k_min, dk, nbin, boxsize, py_k,
				    py_Pdd, py_Pdt, py_Ptt, py_AB),
			 "_Model", py_model_free);
  }
  catch(TypeError) {
    return NULL;
  }

  return NULL;
}



//
// Model implementation
//
Model::Model(const double k_min_, const double dk_, const int nbin_,
	     const double boxsize_) :
  k_min(k_min_), dk(dk_), boxsize(boxsize_), nbin(nbin_)
{
  // k, mu on 3D grid points
  modes= multipole_construct_discrete_wavevectors(k_min, dk, nbin, boxsize);

  // Compute coefficients of discrete Legendre polynomials `coef`
  multipole_compute_discrete_legendre(k_min, dk, nbin,
				      boxsize, coef);
}

Model::~Model()
{
  delete [] modes;
}

void Model::nmodes(vector<double>& v_nmodes) const
{
  size_t n= static_cast<size_t>(nbin);
  assert(v_nmodes.size() == n);

  for(size_t i=0; i<n; ++i) {
    int nmodes= 0;
    
    for(vector<DiscreteWaveVector>::iterator p= modes[i].begin();
	p != modes[i].end(); ++p) {
      nmodes += p->w;
    }

    v_nmodes[i]= nmodes;
  }
}



//
// Kaiser
// continious integral
/*
Kaiser::Kaiser(const double k_min_, const double k_max_, const double dk_,
	       PyObject* py_k, PyObject* py_P) :
  Model(k_min_, k_max_, dk_)
{
  py_util_array_as_vector("k", py_k, v_k);
  py_util_array_as_vector("P", py_P, v_P, v_k.size());
}

void Kaiser::evaluate(const vector<double>& params,
		      vector<double>& v_P0,
		      vector<double>& v_P2,
		      vector<double>& v_P4) const
{
  const size_t n= v_P.size();
  assert(v_P0.size() == n);
  assert(v_P2.size() == n);
  assert(v_P4.size() == n);
  assert(params.size() == 3);

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
*/

Kaiser::Kaiser(const double k_min_, const double dk_, const int nbin_,
	       const double boxsize_,
	       PyObject* py_k, PyObject* py_P) :
  Model(k_min_, dk_, nbin_, boxsize_)
{
  PowerSpectrum ps(py_k, py_P);

  // Set realspace power spectrum to modes
  // this->modes is setup by Model::Model constructor
  for(int ibin=0; ibin<nbin; ++ibin) {
    for(vector<DiscreteWaveVector>::iterator p= modes[ibin].begin();
	p != modes[ibin].end(); ++p) {
      p->Pdd= ps.P(p->k);
    }
  }
}

void Kaiser::evaluate(const vector<double>& params,
		      vector<double>& v_P0,
		      vector<double>& v_P2,
		      vector<double>& v_P4) const
{
  size_t n= static_cast<size_t>(nbin);
  assert(v_P0.size() == n);
  assert(v_P2.size() == n);
  assert(v_P4.size() == n);
  assert(params.size() == 3);

  const double b= params[0];
  const double f= params[1];
  const double s= params[2];

  const double b2= b*b;
  const double bf2= 2.0*b*f;
  const double ff= f*f;

  for(size_t i=0; i<n; ++i) {
    int nmodes= 0;
    double P0= 0.0;
    double P2= 0.0;
    double P4= 0.0;
    
    for(vector<DiscreteWaveVector>::iterator p= modes[i].begin();
	p != modes[i].end(); ++p) {
      double k= p->k;
      double mu2= p->mu2;
      double mu4= mu2*mu2;

      nmodes += p->w;

      double P= p->w*p->Pdd*exp(-k*k*s*s*mu2);
      double l2= coef[5*i] + coef[5*i + 1]*mu2;
      double l4= coef[5*i + 2] + coef[5*i + 3]*mu2
	         + coef[5*i + 4]*mu4;
      
      P0 += (b2 + bf2*mu2 + ff*mu4)*P;
      P2 += l2*(b2 + bf2*mu2 + ff*mu4)*P;
      P4 += l4*(b2 + bf2*mu2 + ff*mu4)*P;
    }

    double fac= 1.0/static_cast<double>(nmodes);

    
    // (2l + 1) for multipole moments
    v_P0[i]= fac*P0;
    v_P2[i]= 5.0*fac*P2;
    v_P4[i]= 9.0*fac*P4;
  }
}

//
// Scoccimarro
//
Scoccimarro::Scoccimarro(const double k_min_, const double dk_, const int nbin_,
			 const double boxsize_,
			 PyObject* py_k,
			 PyObject* py_Pdd, PyObject* py_Pdt, PyObject* py_Ptt) :
  Model(k_min_, dk_, nbin_, boxsize_)
{
  PowerSpectrum pdd(py_k, py_Pdd);
  PowerSpectrum pdt(py_k, py_Pdt);
  PowerSpectrum ptt(py_k, py_Ptt);

  // Set realspace power spectrum to modes
  // this->modes is setup by Model::Model constructor
  for(int ibin=0; ibin<nbin; ++ibin) {
    for(vector<DiscreteWaveVector>::iterator p= modes[ibin].begin();
	p != modes[ibin].end(); ++p) {
      p->Pdd= pdd.P(p->k);
      p->Pdt= pdt.P(p->k);
      p->Ptt= ptt.P(p->k);
    }
  }   
}


void Scoccimarro::evaluate(const vector<double>& params,
			   vector<double>& v_P0,
			   vector<double>& v_P2,
			   vector<double>& v_P4) const
{
  size_t n= static_cast<size_t>(nbin);
  assert(v_P0.size() == n);
  assert(v_P2.size() == n);
  assert(v_P4.size() == n);
  assert(params.size() == 3);

  const double b= params[0];
  const double f= params[1];
  const double s= params[2];

  const double b2= b*b;
  const double bf2= 2.0*b*f;
  const double ff= f*f;

  for(size_t i=0; i<n; ++i) {
    int nmodes= 0;
    double P0= 0.0;
    double P2= 0.0;
    double P4= 0.0;
    
    for(vector<DiscreteWaveVector>::iterator p= modes[i].begin();
	p != modes[i].end(); ++p) {
      double k= p->k;
      double mu2= p->mu2;
      double mu4= mu2*mu2;

      nmodes += p->w;

      double Pdd= p->w*p->Pdd;
      double Pdt= p->w*p->Pdt;
      double Ptt= p->w*p->Ptt;

      double Ps= (b2*Pdd + bf2*mu2*Pdt + ff*mu4*Ptt)*exp(-k*k*s*s*mu2);
	
      double l2= coef[5*i] + coef[5*i + 1]*mu2;
      double l4= coef[5*i + 2] + coef[5*i + 3]*mu2
	         + coef[5*i + 4]*mu4;


      P0 += Ps;
      P2 += l2*Ps;
      P4 += l4*Ps;
    }

    double fac= 1.0/static_cast<double>(nmodes);

    
    // (2l + 1) for multipole moments
    v_P0[i]= fac*P0;
    v_P2[i]= 5.0*fac*P2;
    v_P4[i]= 9.0*fac*P4;
  }
}


//
// Taruya
//
Taruya::Taruya(const double k_min_, const double dk_, const int nbin_,
	       const double boxsize_,
	       PyObject* py_k,
	       PyObject* py_Pdd, PyObject* py_Pdt, PyObject* py_Ptt,
	       PyObject* py_AB) :
  Model(k_min_, dk_, nbin_, boxsize_)
{
  PowerSpectrum pdd(py_k, py_Pdd);
  PowerSpectrum pdt(py_k, py_Pdt);
  PowerSpectrum ptt(py_k, py_Ptt);

  // Set realspace power spectrum to modes
  // this->modes is setup by Model::Model constructor
  for(int ibin=0; ibin<nbin; ++ibin) {
    for(vector<DiscreteWaveVector>::iterator p= modes[ibin].begin();
	p != modes[ibin].end(); ++p) {
      p->Pdd= pdd.P(p->k);
      p->Pdt= pdt.P(p->k);
      p->Ptt= ptt.P(p->k);
    }
  }

  // Load precomputed Taruya AB terms
  vector<double> v_AB;
  size_t ncol= py_util_array2_as_vector("TaruyaAB", py_AB, v_AB);
  if(ncol != 15) throw TypeError();
  size_t nrow = v_AB.size() / ncol;

  Interp interp(v_AB, nrow, ncol);

  // Set TaruyaAB to `this->modes`
  for(int ibin=0; ibin<nbin; ++ibin) {
    for(vector<DiscreteWaveVector>::iterator p= modes[ibin].begin();
	p != modes[ibin].end(); ++p) {
      // The order of Bnab is different from original Taruya's Fortran
      // code output. Also, our Bnab does not contain sign (-1)^{a + b}
       
      const double k= p->k;
      p->A11=  interp(1,  k);
      p->A12=  interp(2,  k);
      p->A22=  interp(3,  k);
      p->A23=  interp(4,  k);
      p->A33=  interp(5,  k);
      p->B111= interp(6,  k);
      p->B211= interp(7,  k);
      p->B112= interp(8,  k);
      p->B212= interp(9,  k);
      p->B312= interp(10, k);
      p->B122= interp(11, k);
      p->B222= interp(12, k);
      p->B322= interp(13, k);
      p->B422= interp(14, k);
    }
  }
}


void Taruya::evaluate(const vector<double>& params,
			   vector<double>& v_P0,
			   vector<double>& v_P2,
			   vector<double>& v_P4) const
{
  size_t n= static_cast<size_t>(nbin);
  assert(v_P0.size() == n);
  assert(v_P2.size() == n);
  assert(v_P4.size() == n);
  assert(params.size() == 3);

  const double b= params[0];
  const double f= params[1];
  const double s= params[2];

  const double b2= b*b;
  const double bf2= 2.0*b*f;
  const double f2= f*f;
  const double f3= f*f2;
  const double f4= f2*f2;

  for(size_t i=0; i<n; ++i) {
    int nmodes= 0;
    double P0= 0.0;
    double P2= 0.0;
    double P4= 0.0;
    
    for(vector<DiscreteWaveVector>::iterator p= modes[i].begin();
	p != modes[i].end(); ++p) {
      double k= p->k;
      double mu2= p->mu2;
      double mu4= mu2*mu2;
      double mu6= mu4*mu4;
      double mu8= mu4*mu4;

      nmodes += p->w;

      double Pdd= p->w*p->Pdd;
      double Pdt= p->w*p->Pdt;
      double Ptt= p->w*p->Ptt;
      double A= p->w*(p->A11*f*mu2
		      + f2*(p->A12*mu2 + p->A22*mu4)
		      + f3*(p->A23*mu4 + p->A33*mu6));
      double B= p->w*(f2*(p->B111*mu2 + p->B211*mu4)
		      - f3*(p->B112*mu2 + p->B212*mu4 + p->B312*mu6)
		      + f4*(p->B122*mu2 + p->B222*mu4 + p->B322*mu6
			    + p->B422*mu8));
      // (-1)^{a + b} sign is here not in the file
      
      double Ps= (b2*Pdd + bf2*mu2*Pdt + f2*mu4*Ptt + A + B)*exp(-k*k*s*s*mu2);
	
      double l2= coef[5*i] + coef[5*i + 1]*mu2;
      double l4= coef[5*i + 2] + coef[5*i + 3]*mu2
	         + coef[5*i + 4]*mu4;
      
      P0 += Ps;
      P2 += l2*Ps;
      P4 += l4*Ps;
    }

    double fac= 1.0/static_cast<double>(nmodes);

    
    // (2l + 1) for multipole moments
    v_P0[i]= fac*P0;
    v_P2[i]= 5.0*fac*P2;
    v_P4[i]= 9.0*fac*P4;
  }
}
