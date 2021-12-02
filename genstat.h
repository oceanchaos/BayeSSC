#ifndef __GENSTAT_H__
#define __GENSTAT_H__

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "deme.h"
#include "arrays.h"
#include "cstring.h"
#include "algebra.h"

typedef MY_TArrayAsVector<int> THIST;

struct summary_stats {
	int num_haps;
	int num_private_haps;
 	int seg_sites;
  	double pairws_dif;
   double hap_div;
   double nuc_div;
   double tajd;
   double fstar;
   //char hap_freq[20];
   //char mut_size[20];
   THIST mismatch;
};

class T_HAP {
  public:
	TIntVect* hap;
	int n;
	int first_deme;
	bool is_private;

	~T_HAP() {;};
   T_HAP() { n=0; first_deme=-1; is_private=1; }
   int operator==(const T_HAP& H) {
		int *site1=hap->begin(), *site2=H.hap->begin();
    	while(site1<hap->end())
    		if(*site1++!=*site2++) return 0;
      return 1;
   }
   int operator<(const T_HAP& H) {return 0; }
   void print() {
		printf("\nHAPLOTYPE: %d copies, first found in deme %d (%sPrivate)\n",n,first_deme,(is_private)?"":"Not");
		for(int *i=hap->begin();i<hap->end();i++) printf("%d",*i);
	}
};

enum e_param {deme_size, sample_age, sample_stat_grp, growth_rate, 
   mig_rate, event_time, event_src, event_sink, event_migrants, event_size,
   event_growth, event_mig_matr, pmut_rate, abstract, none};
enum e_dist {expression, uniform, normal, gammadist, exponential, geometric};
void param2str(char *x, e_param p);
char dist2str(e_dist d);
double str2double(char *x);
double rgamma(int k, double theta);

//Christian: a class for parameter approximation
class TPrior {
   public :
      double min, max;
      double curval; 
      int type; //0=none 1=int 2=long 3=float 4=double
      e_param param;
      e_dist dist;
      int num;
      symb *formula;
      
      TPrior(double mn, double mx) {
         min=mn; max=mx; type=0; param=none; dist=expression; num=0; formula=NULL;
      }
      TPrior(void) {
         min=0; max=0; type=0; param=none; dist=expression; num=0; formula=NULL;
      }
      void Print(void) { 
         char x[20]; param2str(x,param);
         printf("%s %d: %f, %c[%f,%f]", x, num, curval, dist2str(dist), min, max);
         if(formula) formula->print();
      }
      void pick(void *memloc);
      void set(void *memloc) {
         switch(type) {
            case 1: *(int *)memloc=int(curval+.5); break;
            case 2: *(long *)memloc=long(curval+.5); break;
            case 3: *(float *)memloc=float(curval); break;
            case 4: *(double *)memloc=curval; break;
            default: printf("\nERROR: Bad data type when picking value from prior."); break;
         }
      }
};
TPrior* BayesRead(void *buf, istream *ifs, e_param p, int n, int type);

typedef MY_TSArrayAsVector<T_HAP>  THapArray;

void CalcStats(THapArray *haps, summary_stats *st, const Mut_Type data_type);
void CalcInterStats(THapArray *h1, THapArray *h2, summary_stats *st, const Mut_Type data_type);

#endif
