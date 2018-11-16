#ifndef PY_MULTIPOLE_H
#define PY_MULTIPOLE_H 1

#include <vector>

using namespace std;

struct DiscreteWaveVector {
  //int ikx, iky, ikz;
  double k, mu2;
  int w;
  double Pdd, Pdt, Ptt;
};

std::vector<DiscreteWaveVector>*
multipole_construct_discrete_wavevectors(
    const double k_min, const double dk, const int nbin,
    const double boxsize);

void multipole_compute_discrete_legendre(
    const double k_min, const double dk, const int nbin,
    const double boxsize, std::vector<double>& coef);

#endif

