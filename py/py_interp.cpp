#include "py_interp.h"
#include <cstdlib>

using namespace std;

void copy_column(const vector<double>& v,
		   const size_t nrow, const size_t ncol,
		   const size_t icol, vector<double>& v_col)
{
  assert(v.size() == nrow*ncol);
  v_col.clear();
  v_col.reserve(nrow);

  for(size_t i=0; i<nrow; ++i)
    v_col.push_back(v[ncol*i + icol]);
}

Interp::Interp(const vector<double>& v, const size_t nrow, const int ncol_) :
  ncol(ncol_)
{
  spline= (gsl_spline**) malloc(sizeof(gsl_spline*)*(ncol - 1));
  assert(spline);
  acc = (gsl_interp_accel**) malloc(sizeof(gsl_interp_accel*)*(ncol - 1));
  assert(acc);

  assert(v.size() == nrow*ncol);

  vector<double> v_k, v_col;
  copy_column(v, nrow, ncol, 0, v_k);
  
  for(size_t icol=1; icol<ncol; ++icol) {
    copy_column(v, nrow, ncol, icol, v_col);

    size_t i= icol - 1;
    acc[i]= gsl_interp_accel_alloc();
    spline[i]= gsl_spline_alloc(gsl_interp_cspline, nrow);
    gsl_spline_init(spline[i], v_k.data(), v_col.data(), nrow);
  }
}

Interp::~Interp()
{
  for(size_t icol=1; icol<ncol; ++icol) {
    gsl_spline_free(spline[icol - 1]);
    gsl_interp_accel_free(acc[icol - 1]);
  }
  free(acc);
  free(spline);
}

double Interp::operator()(const size_t icol, const double k)
{
  assert(1 <= icol && icol < ncol);

  return gsl_spline_eval(spline[icol - 1], k, acc[icol - 1]);
}
