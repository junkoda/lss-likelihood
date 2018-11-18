#include <vector>
#include <cmath>
#include <cassert>


#include "py_util.h"
#include "py_error.h"
#include "py_assert.h"
#include "py_multipole.h"
#include "py_interp.h"
#include "py_model.h"
#include "py_better_model.h"

using namespace std;

// Better-than-nothing model
PyObject* py_better_model_alloc(PyObject* self, PyObject* args)
{
  // _model_scoccimarro_alloc(imodel,
  //                          k_min, dk, nbin, boxsize, k, Pdd, Pdt, Ptt, AB)
  int imodel;
  double k_min, dk, boxsize;
  int nbin;
  PyObject *py_k, *py_Pdd, *py_Pdt,*py_Ptt, *py_AB;
    
  if(!PyArg_ParseTuple(args, "iddidOOOOO", &imodel,
		       &k_min, &dk, &nbin, &boxsize,
		       &py_k, &py_Pdd, &py_Pdt, &py_Ptt, &py_AB))
    return NULL;

  //try {
  switch(imodel) {
  case 1:
    return PyCapsule_New(new Model1(k_min, dk, nbin, boxsize, py_k,
				    py_Pdd, py_Pdt, py_Ptt, py_AB),
			 "_Model", py_model_free);
  default:
    char msg[128];
    sprintf(msg, "Unknown better_model number: %d", imodel);
    PyErr_SetString(PyExc_ValueError, msg);
    return NULL;
  }
    //}
    //catch(TypeError) {
    //return NULL;
    //}

  return NULL;
}


//
// Model1
//
Model1::Model1(const double k_min_, const double dk_, const int nbin_,
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


void Model1::evaluate(const vector<double>& params,
			   vector<double>& v_P0,
			   vector<double>& v_P2,
			   vector<double>& v_P4) const
{
  size_t n= static_cast<size_t>(nbin);
  assert(v_P0.size() == n);
  assert(v_P2.size() == n);
  assert(v_P4.size() == n);
  assert(params.size() == 4);

  const double b= params[0];
  const double f= params[1];
  const double sigma0= params[2];
  const double sigma1= params[3];

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
      double s= sigma0 + sigma1*k;
	  

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
