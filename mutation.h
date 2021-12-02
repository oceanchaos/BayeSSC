#ifndef _MUTATION_H_
#define _MUTATION_H_

#include "public.h"

class GammaRates{
   private :
   	double * _vecMut;
      double _a;
      double _ratio;
      int _fast_sites;
      int _size;
      int _rate_categories;

	public :
   	GammaRates() {
       	_size=0;
         _a=0.0;
         _fast_sites=0;
         _ratio=1.0;
         _rate_categories=0; 
         _vecMut=NULL;
      }
      GammaRates(const int& n, const double& a, const double& r, const double& p, const int& nG) {
       	_size=n;
         _vecMut=NULL;
         _a=a;
         _fast_sites= (int) p*n;
         _ratio=r;
         _rate_categories=nG;
      	if (n>0) {
      		if (!build(n,a,r,p, nG))
            	cout << "\nGammaRates(const int n, const double& a): Unable to build gamma deviate vector\n";
         }
      };
      ~GammaRates() {
      	if (_vecMut) delete[] _vecMut; _vecMut=NULL;
      }
      int build(const int n, const double& a, const float& mut_ratio,
      				const float& prop_sites, const int& num_rates);
      int get_new_rates(const int n, const double& a, const double& r,
      						const double& p, const int& nG);
      double& operator[](const int& i) {return _vecMut[i];} 
      const int& size() {return _size;}
      const double& a() {return _a;} 

};

#endif

