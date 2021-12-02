#ifndef __PUBLIC_H__
#define __PUBLIC_H__

#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "cstring.h"
#include "arrays.h"

typedef double my_float; 
                        
enum Mut_Type {MICROSAT, RFLP, DNA};

extern float mtr(void);
extern float gammln(float xx);

// A routine to generate poisson deviates with mean xm 
extern float poidev(float xm);
extern void quicksort(int l, int r, double* list);

// A routine to initialize the probabilities for the Gamma distributed mutation rates.
void init_gamma_weights(double* p, int n_loci, double alph);
//Discrete gamma case, using Ziheng Yang's routines
void init_gamma_weights(double* p, int n_loci, double alph, const int num_rates);


//returns gamma deviates
extern double gamma_dev(double a);
extern double igamma_dev(int ia);
typedef MY_TArrayAsVector<my_string> StringArray;

//Arlequin output routines
//------------------------------------------------------------------------------
extern void write_Arlequin_header(const int& num_samples,const Mut_Type& data_type, ostream& os);
extern void write_Arlequin_sample_header(const int& num_pop, const int& size,ostream & os);
extern void write_Arlequin_sample_footer(ostream & os);

//PAUP* output routines
//------------------------------------------------------------------------------

enum Mut_Model {K80_GAMMA, K80_NOGAMMA, HKY85_GAMMA, HKY85_NOGAMMA/*, JK, TN93, HK85*/};

extern void write_PAUP_header(ostream& os, const char* logFile, int tot_num_nodes,
   int num_loci,int num_pop, double mut_rate, double gamma_a, double transition_rate);
extern void write_PAUP_block(ostream & os, const Mut_Model mut_mod);
extern void write_PAUP_footer(ostream & os);

extern void write_PAUP_data_header(const int ntax, const int nchar, Mut_Type datatype, ostream & os);
extern void write_PAUP_end(ostream & os);

extern void write_PAUP_end_matrix(ostream & os);

extern void write_PAUP_tree_header(ostream & os, const int num_rep);

extern void write_PAUP_replicate_number(const int n, ostream & os);
extern void write_PAUP_tree_name(ostream & os, const int num_rep);
extern void write_PAUP_trees_section_header(ostream & os);

//Ziheng Yang routines for generating discrete gamma distributions
//------------------------------------------------------------------------------
extern int Rates4Sites (double rates[],double alpha,int ncatG,int ls, int cdf, double space[]);
extern double LnGamma (double alpha);
extern double DFGamma(double x, double alpha, double beta);
extern double IncompleteGamma (double x, double alpha, double ln_gamma_alpha);
extern double PointChi2 (double prob, double v);
extern double PointNormal (double prob);
extern int DiscreteGamma (double freqK[], double rK[], double alfa, double beta, int K, int median);
extern double rndgamma (double s);
extern int MultiNomial (int n, int ncat, double prob[], int nobs[], double space[]);
extern double rndu (void);
extern int xtoy (double x[], double y[], int n);
extern int abyx (double a, double x[], int n);
//------------------------------------------------------------------------------

#endif

