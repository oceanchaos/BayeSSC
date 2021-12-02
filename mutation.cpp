#include "mutation.h"
//#include <except.h>

#include "cond_var.h"

#ifdef _SHOW_GAMMA_DISTR_
	#include <fstream>
	ofstream f_gamma_distr("gamma_rates.distr");
#endif

static int count_sim=0;

int
GammaRates::build(const int n, const double& a, const float& mut_ratio,
      				const float& prop_sites, const int & num_rates) {
	_size=n;
   _a=a;
   _fast_sites=(int)prop_sites*n;
   _ratio=mut_ratio;
   if (n)
   try {
   	if (_vecMut) delete[] _vecMut; _vecMut=NULL;
      _vecMut = new double[n];
      double sum_gamma_deviates=0.0;    
      double * next=&_vecMut[0];
      //Loro_27_8_98
      if (a>0.0) {
      	//Costly routine. Make sure it is called only once in the program
         if (num_rates==0) init_gamma_weights(_vecMut, n, a);
         else init_gamma_weights(_vecMut, n, a, num_rates);
         for (int i=0; i<n; ++i, ++next) {
          	//sum_gamma_deviates+=_vecMut[i];
          	sum_gamma_deviates+=*next;
         }
      }
      else {
      	for (int i=0; i<n; ++i, ++next) {
             //Assumes a two-rate mutation model
         	/*if (i<_fast_sites) _vecMut[i]=mut_ratio;
         	else _vecMut[i]=1.0;
          	sum_gamma_deviates+=_vecMut[i];
         */
         	if (i<_fast_sites) *next=mut_ratio;
         	else *next=1.0;
          	sum_gamma_deviates+=*next;
         }
      }
      #ifdef _SHOW_GAMMA_DISTR_
      f_gamma_distr << "Simulation #" << ++count_sim << " a=" << a << "   mut_ratio=" << mut_ratio << "\n";
      #endif
      //Normalization, just to be on the safe side
      next=&_vecMut[0];
      for (int i=0; i<n; ++i, ++next) {
      	//_vecMut[i]/=sum_gamma_deviates;
         *next/=sum_gamma_deviates;
         #ifdef _SHOW_GAMMA_DISTR_
         //f_gamma_distr << _vecMut[i] << "  ";
         f_gamma_distr << *next << "\n";
         #endif
      }
      #ifdef _SHOW_GAMMA_DISTR_
      f_gamma_distr << "\n";
      #endif
   }
   catch (...) {
   	if (_vecMut) delete _vecMut; _vecMut=NULL;
      _size=0;
      cout << "\nGammaRates::build(const int n, const double& a): unable to allocate memory\n";
   	return 0;
   }
   return _size;
};

int
GammaRates::get_new_rates(const int n, const double& a, const double& r,
									const double& p, const int& num_rates) {

   _fast_sites= (int) (n*p);

	if (n!=_size) {
   	build(n, a, r, p, num_rates);
      return _size;
   }

	double sum_gamma_deviates=0.0;
   double *next=&_vecMut[0];
   //Loro_27_8_98
      if (a>0.0) {
      	//Costly routine. Make sure it is called only once in the program 
         if (num_rates==0) init_gamma_weights(_vecMut, n, a);
         else init_gamma_weights(_vecMut, n, a, num_rates);
         for (int i=0; i<n; ++i, ++next) {
          	//sum_gamma_deviates+=_vecMut[i];
         	sum_gamma_deviates+=*next;
         }
      }
      else {
      	for (int i=0; i<n; ++i, ++next) {
             //Assumes a two-rate mutation model
         	//if (i<_fast_sites) _vecMut[i]=_ratio;
         	if (i<_fast_sites) *next=_ratio;
         	else *next=1.0;
          	//sum_gamma_deviates+=_vecMut[i];
          	sum_gamma_deviates+=*next;
         }
      }
   #ifdef _SHOW_GAMMA_DISTR_
   f_gamma_distr << "Simulation #" << ++count_sim << " a=" << a << "   mut_ratio=" << r << "  ";
   #endif
   //Normalization, just to be on the safe side
   next=&_vecMut[0];
   for (int i=0; i<_size; ++i, ++next) {
   	//_vecMut[i]/=sum_gamma_deviates;
      *next/=sum_gamma_deviates;
      #ifdef _SHOW_GAMMA_DISTR_
      //f_gamma_distr << _vecMut[i] << "  ";
      f_gamma_distr << *next << "  ";
      #endif
   }
   #ifdef _SHOW_GAMMA_DISTR_
   f_gamma_distr << "\n";
   #endif
   
   return _size;
}
