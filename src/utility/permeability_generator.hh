#ifndef __PERMEABILITY_GENERATOR_HH__
#define __PERMEABILITY_GENERATOR_HH__

// C++ includes
#include<iostream>
#include<cstdlib>
#include<algorithm>
#include<vector>
#include<deque>
#include<valarray>
#include<string>
#include<cmath>

// C includes
#include<string.h>

template<int dim>
class EberhardPermeabilityGenerator
{
public:

  EberhardPermeabilityGenerator (Dune::FieldVector<double,dim> corr_vec, double n_variance = 1.0, double n_mean = 0.0, 
								 long num_of_modes = 1000, long sd = -1083)
  {
	space_dim = dim;
	IA = 16807;
	IM = 2147483647;
	AM = (1.0/IM);
	IQ = 127773;
	IR = 2836;
	NDIV = (1+(IM-1)/NTAB);
	EPS = 1.2e-7;
	RNMX = (1.0-EPS);

	corr_length_vec.resize(dim);
	for (int i=0; i<dim; i++) corr_length_vec[i] = corr_vec[i];
	number_of_modes = num_of_modes;
	seed = sd;
	normal_mean = n_mean;
	normal_variance = n_variance;
	alpha.resize(num_of_modes);
	q_vec.resize(dim*num_of_modes);
	norm_factor = sqrt (2.0 * n_variance / num_of_modes);

	for (long n = 0; n < number_of_modes; ++n)
	  {
		alpha[n] = 2.0 * M_PI * ran1 (&seed);
	  }
	for (long d = 0; d < space_dim; ++d)
	  {
		for (long n = 0; n < number_of_modes; ++n)
		  {
			q_vec[n*space_dim + d] = gasdev (&seed) / corr_length_vec[d];
		  }
	  }
  }

  template<typename RF>
  double eval (const Dune::FieldVector<RF,dim>& x) const
  {
	std::valarray<double> space_vec(dim);
	for (int i=0; i<dim; i++) space_vec[i] = x[i];
	double sum = 0.0, arg_of_cos = 0.0;
	for (long n = 0; n < number_of_modes; ++n, arg_of_cos = 0.0)
	  {
		for (long d = 0; d < space_dim; ++d) arg_of_cos += q_vec[n*space_dim + d] * space_vec[d];
		sum += cos (arg_of_cos + alpha[n]);
	  }
	return exp (normal_mean + norm_factor * sum);
  }

private:
  double gasdev (long *idum) // Numerical Recipes
  {
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;
	if  (iset == 0) {
	  do {
		v1=2.0*ran1(idum)-1.0;
		v2=2.0*ran1(idum)-1.0;
		rsq=v1*v1+v2*v2;
	  } while (rsq >= 1.0 || rsq == 0.0);
	  fac=sqrt(-2.0*log(rsq)/rsq);
	  gset=v1*fac;
	  iset=1;
	  return v2*fac;
	} else {
	  iset=0;
	  return gset;
	}
  }

  double ran1 (long *idum)
  {
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;
    
	if (*idum <= 0 || !iy) {
	  if (-(*idum) < 1) *idum=1;
	  else *idum = -(*idum);
	  for (j=NTAB+7;j>=0;j--) {
		k=(*idum)/IQ;
		*idum=IA*(*idum-k*IQ)-IR*k;
		if (*idum < 0) *idum += IM;
		if (j < NTAB) iv[j] = *idum;
	  }
	  iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
  }

  std::valarray<double> corr_length_vec; //vector of correlation lengths for d dimensions
  long number_of_modes;
  long seed; // must be negative
  double normal_mean; //mean of the normal field
  double normal_variance;
  double norm_factor;
  std::valarray<double> q_vec;
  std::valarray<double> alpha;
  long space_dim;
  long IA;
  long IM;
  double AM;
  long IQ;
  long IR;
  static const long NTAB=32;
  double NDIV;
  double EPS;
  double RNMX;
};

#endif
