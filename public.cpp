#include "public.h"
#include "MersenneTwister.h"
static MTRand rangen;

float mtr(void) { return float(rangen()); }

/* (C) Copr. 1986-92 Numerical Recipes Software Y5jc. */
//------------------------------------------------------------------------------
float gammln(float xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
/* (C) Copr. 1986-92 Numerical Recipes Software Y5jc. */
//------------------------------------------------------------------------------
//routine to initialise the probabilities for the Gamma distributed mutation rates.
// 1) we generate n_gen random values according to the gemma density function by using
//    the routine gamma_dev. Let the vector be v_gamma_dev
// 2) we sort the resulting values and make n_loci classes
// 3) Pi = v_gama_dev[ni] where ni is the center of the i'th class

void init_gamma_weights(double* p, int n_loci, double alph){

	//ofstream ofs("test_gam.txt");
	int n_gen=1000000;
   double *v_gama_dev;
   try {
   	v_gama_dev=new double[n_gen];
   }
   catch (...) {
   	cout << "Memory allocation error, exiting ...\n";
      exit(0);
   }
   double *next=&v_gama_dev[0];
   for(int i=0;i<n_gen;i++, ++next){
   	//v_gama_dev[i]=(float) gamma_dev(alph);
   	*next=(double) gamma_dev(alph);
   }
   //sort the elements
   quicksort(0,n_gen-1, v_gama_dev);
   /*
   //Check that the elements are well sorted
   ofs << "Sorted gamma vector: " << endl;
   for(int i=0;i<n_gen;i++){
   	//ofs << v_gama_dev[i] << endl;
      if(i<n_gen-1){
      	if( v_gama_dev[i] > v_gama_dev[i+1] ) ofs << "Error in sorting the vector";
      }
   }
   */
   //compute the Pi:
   int num_in_cl=n_gen/n_loci;
   int count=0;
   for(int i= num_in_cl/2; i<n_gen && count<n_loci; i+=num_in_cl){
   	p[count]=v_gama_dev[i];
      count++;
   }
   /*
   for(int i=0;i< n_loci; i++){
   	ofs << "P[" << i << "]=\t" << p[i] << endl;
   }
   */
   delete[]  v_gama_dev;
}
//------------------------------------------------------------------------------
void init_gamma_weights(double* p, int n_loci, double alph, const int num_rates){

	//ofstream ofs("test_gam.txt");
	int n_gen=n_loci;
   double *v_gama_dev, *working_space;
   try {
   	v_gama_dev=new double[n_gen];
      working_space= new double[4*num_rates];
   }
   catch (...) {
   	if (v_gama_dev) delete[] v_gama_dev;
   	if (working_space) delete[] working_space;
   	cout << "Memory allocation error, exiting ...\n";
      exit(0);
   }

   Rates4Sites (v_gama_dev, alph, num_rates , n_loci, 0, working_space);

   //sort the elements
   quicksort(0,n_gen-1, v_gama_dev);
/*
   //Computes the relative rates
   double sum_rel_rates=0.0;
   for(int i=0;i<n_gen;i++){
     	sum_rel_rates+=v_gama_dev[i];
   }

   if (sum_rel_rates>0)
   for(int i=0; i<n_gen ;++i){
   	p[i]=v_gama_dev[i]/sum_rel_rates;
   }
   else for(int i=0; i<n_gen ;++i){
   	p[i]=0.0;
   }
*/
	for(int i=0; i<n_gen ;++i){
   	p[i]=v_gama_dev[i];
   }
   /*
   for(int i=0;i< n_loci; i++){
   	ofs << "P[" << i << "]=\t" << p[i] << endl;
   }
   */
   if (v_gama_dev) delete[]  v_gama_dev;
   if (working_space) delete[] working_space;
}
//------------------------------------------------------------------------------
void quicksort(int l, int r, double* list) {
 	int i, j;
   double temp, v;
   if (r>l) {
   	v=list[r]; i=l-1; j=r;
      for (;;) {
      	while (list[++i] < v) ;
         while (list[--j] > v) ;
         if (i>=j) break;
         //swap elements;
         temp=list[i];
         list[i]=list[j];
         list[j]=temp;
      }
      //swap elements;
      temp=list[i];
      list[i]=list[r];
      list[r]=temp;
      quicksort(l, i-1, list);
      quicksort(i+1, r, list);
   }
}
//------------------------------------------------------------------------------

#define PI 3.141592654

float poidev(float xm) {
	float gammln(float xx);
	static float sq,alxm,g,oldm=(-1.0);
	float em,t,y;

	if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
			++em;
			t *= mtr();
		} while (t > g);
	} else {
		if (xm != oldm) {
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-gammln(xm+1.0);
		}
		do {
			do {
				y=tan(PI*mtr());
				em=sq*y+xm;
			} while (em < 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		} while (mtr() > t);
	}
	return em;
}
#undef PI
/* (C) Copr. 1986-92 Numerical Recipes Software Y5jc. */

/*
Rogers comments:

In answer to your question about my algorithm, I'm going to append my
entire gamma_dev function.  The comments at the top provide references
to the original sources.  The algorithm I use is supposedly the most
commonly used when alpha<1.

In case it is relevant, let me tell you about some of the trouble I've
run into generating gamma deviates with small values of alpha.  My
first gamma_dev function was in single precision.  It behaved very
strangely.  When alpha<0.1, the number of segregating sites went *up*
as alpha went *down*, which makes no sense at all.  I couldn't find
any error in the code, but I noticed that the code does things that
may stretch the limits of floating point arithmetic.  So I recompiled
using double precision for all variables within gamma_dev.  The
strange behavior went away.

The literature doesn't say much about the stability of these
algorithms when alpha is very small.  It seems that no one has ever
been interested in that case.  I'll bet that none of the commercial
statistical packages have tested their gamma deviate generator with
very small alpha values either.  Consequently, we can't test our
algorithms by comparing the quantiles of our generated values with
those generated by, say, SPSS.  The only sure way is to calculate
quantiles by direct integration of the density function.  I have done
this for alpha=0.1 and am about to compare the quantiles of my numbers
with these values.  I'll let you know what happens.

Alan

PS  Here's the code along with references.  */

/****************************************************************
Random deviates from standard gamma distribution with density
         a-1
        x    exp[ -x ]
f(x) = ----------------
         Gamma[a]

where a is the shape parameter.  The algorithm for integer a comes
from numerical recipes, 2nd edition, pp 292-293.  The algorithm for
a<1 uses code from p 213 of Statistical Computing, by Kennedy and
Gentle, 1980 edition.  This algorithm was originally published in:

Ahrens, J.H. and U. Dieter (1974), "Computer methods for sampling from
Gamma, Beta, Poisson, and Binomial Distributions".  COMPUTING
12:223-246.

The mean and variance of these values are both supposed to equal a.
My tests indicate that they do.

This algorithm has problems when a is small.  In single precision, the
problem  arises when a<0.1, roughly.  That is why I have declared
everything as double below.  Trouble is, I still don't know how small 
a can be without causing trouble.  Mean and variance are ok at least
down to a=0.01.  f(x) doesn't seem to have a series expansion around
x=0.
****************************************************************/
double gamma_dev(double a) {

  int ia;
  double u, b, p, x, y=0.0, recip_a;

  if(a <= 0) {
    cout << "\ngamma_dev: parameter must be positive\n";
    exit(1);
  }

  ia = (int) (floor(a));  /* integer part */
  a -= ia;        /* fractional part */
  if(ia > 0) {
    y = igamma_dev(ia);  /* gamma deviate w/ integer argument ia */
    if(a==0.0) return(y);
  }

  /* get gamma deviate with fractional argument "a" */
  b = (M_E + a)/M_E;
  recip_a = 1.0/a;
  for(;;) {
    u = mtr();
    p = b*u;
    if(p > 1) {
      x = -log( (b-p)/a );
      if( mtr() > pow(x, a-1)) continue;
      break;
    }
    else {
      x = pow(p, recip_a);
      if( mtr() > exp(-x)) continue;
      break;
    }
  }
  return(x+y);
}
/****************************************************************
gamma deviate for integer shape argument.  Code modified from pp
292-293 of Numerical Recipes in C, 2nd edition.
****************************************************************/
double igamma_dev(int ia)
{
  int j;
  double am,e,s,v1,v2,x,y;
  
  if (ia < 1)
  {
    cout << "\nError: arg of igamma_dev was <1\n";
    exit(1);
  }
  if (ia < 6)
  {
    x=1.0;
    for (j=0; j<ia; j++)
      x *= mtr();
    x = -log(x);
  }else
  {
    do
    {
      do
      {
	do
	{                         /* next 4 lines are equivalent */
	  v1=2.0*mtr()-1.0;       /* to y = tan(Pi * uni()).     */
	  v2=2.0*mtr()-1.0;
	}while (v1*v1+v2*v2 > 1.0);
	y=v2/v1;
	am=ia-1;
	s=sqrt(2.0*am+1.0);
	x=s*y+am;
      }while (x <= 0.0);
      e=(1.0+y*y)*exp(am*log(x/am)-s*y);
    }while (mtr() > e);
  }
  return(x);
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void
write_Arlequin_header(const int& num_samples, const Mut_Type& data_type, ostream& os) {

	os
   	<< "\n[Profile]\n"
		<< "\tTitle=\"A series of simulated samples\"\n"
     	<< "\tNbSamples=" << num_samples << "\n"
      << "\tGenotypicData=0\n"
      << "\tGameticPhase=1\n"
      << "\tRecessiveData=0\n"
      << "\tDataType=";
   switch (data_type) {
   	case DNA 		: os << "DNA\n" << "\tLocusSeparator=NONE\n"; break;
      case RFLP		: os << "RFLP\n" << "\tLocusSeparator=NONE\n"; break;
      case MICROSAT	: os << "MICROSAT\n" << "\tLocusSeparator=WHITESPACE\n"; break;
   }
   os
      << "\tMissingData='?'\n"
      << "\n[Data]\n"
      << "\t[[Samples]]\n\n";
};

//------------------------------------------------------------------------------
void
write_Arlequin_sample_header(const int& num_pop, const int& size, ostream & os) {
	os
		<< "\t\tSampleName=\"Sample " << num_pop << "\"\n"
      << "\t\tSampleSize=" << size
      << "\n\t\tSampleData= {\n";
};

void write_Arlequin_sample_footer(ostream & os) {
	os << "\n}\n";
}

void write_PAUP_header(ostream& os, const char* logFile, int tot_num_nodes,
   int num_loci,int num_pop, double mut_rate, double gamma_a, double transition_rate) {
   os << "#NEXUS\n"
   	<< "\n[Simulated data generated by program sim_coal.exe (Laurent Excoffier)"
      << "\n modified for ancient DNA by Liz Hadly's Lab, Stanford Univ."
      << "\n\tSimulation conditions:"
      << "\n\tNo. of taxa   : " << tot_num_nodes
      << "\n\tNo. of loci   : " << num_loci
      << "\n\tNo. of demes  : " << num_pop
      << "\n\tMutation rate : " << (mut_rate/num_loci) << " per locus (site) per generation"
      << "\n\tAlpha         : " << gamma_a << " (0=homogeneity)"
      << "\n\t%Transitions  : " << (100*transition_rate) << " (ts/tv=";
   if (transition_rate!=1)
   	os << (transition_rate/(1-transition_rate));
   else os << "inf";
   os << ")\n]\n\n";  

	os << "\nbegin paup; [header to log results into a log file]"
   	<< "\n\tlog start file=" << logFile << " repl=yes;"
		<< "\nend;"
      << endl;
};

void write_PAUP_block(ostream & os, const Mut_Model mut_mod) {
	os
   	<< "\nbegin paup;"
   	<< "\n\tset autoclose=yes warnreset=no increase=auto; [turns off user prompts]"
   	<< "\n\tset crit=likelihood;";
   switch (mut_mod) {
   	case 	K80_GAMMA :
      	os << "\n\tlset nst=2 tratio=est base=equal rates=gamma shape=est ncat=8;";
         break;
      case K80_NOGAMMA :
      	os << "\n\tlset nst=2 tratio=est base=equal rates=equal;";
         break;
   	case 	HKY85_GAMMA :
      	os << "\n\tlset nst=2 tratio=est rates=gamma shape=est ncat=8;";
         break;
      case HKY85_NOGAMMA :
      	os << "\n\tlset nst=2 tratio=est rates=equal;";
         break;
   }
   os	<< "\n\tlscore;"
   	<< "\n\tsavetree brlens=yes maxdec=5 append=yes;"
		<< "\nend;"
      << endl;
};
//------------------------------------------------------------------------------

void

write_PAUP_footer(ostream & os) {

	os
   	<< "\nbegin paup;  [footer to stop]"
		<< "\n\tlog stop;"
      << "\nend;"
		<< "\nquit;"
      << endl;
};

//------------------------------------------------------------------------------
extern void write_PAUP_data_header(const int ntax, const int nchar, Mut_Type datatype, ostream & os) {
	os
   	<< "\nbegin data;"
   	<< "\n\tdimensions ntax=" << ntax << " nchar=" << nchar << ";"
   	<< "\n\tformat datatype=";
   switch (datatype) {
		case DNA 		: os << "dna;"; break;
		case RFLP 		: os << "rflp;"; break;
		case MICROSAT 	: os << "microsat;"; break;
   }
   os << "\n\tmatrix" << endl;
};
//------------------------------------------------------------------------------
void
write_PAUP_tree_header(ostream & os, const int num_rep) {
	os
   	<< "\nbegin tree;"
   	<< "\n\ttree true_tree_" << num_rep << " = [&U] ";
}
//------------------------------------------------------------------------------
void
write_PAUP_end(ostream & os) {
	os << "end;" << endl;
}
//------------------------------------------------------------------------------

void

write_PAUP_end_matrix(ostream & os) {

	os << "\t;\n end;" << endl;

};

//------------------------------------------------------------------------------
void
write_PAUP_replicate_number(const int n, ostream & os) {
	os << "\n[Replicate # " << n << " ]" << endl;
};
//------------------------------------------------------------------------------
void
write_PAUP_tree_name(ostream & os, const int num_rep) {
	os
   	<< "\ttree true_tree_" << num_rep << " = [&U] ";
};
//------------------------------------------------------------------------------
void
write_PAUP_trees_section_header(ostream & os) {
	os
   	<< "#NEXUS\n"
      << "begin trees;  [Treefile generated by sim_coal.exe (Laurent Excoffier)]\n"
      << endl;
};
//------------------------------------------------------------------------------
//Ziheng Yang routines for generating discrete gamma distributions
//------------------------------------------------------------------------------

#define FOR(i,n) for(i=0; i<n; i++)
#define PointGamma(prob,alpha,beta) PointChi2(prob,2.0*(alpha))/(2.0*(beta))

//------------------------------------------------------------------------------
int Rates4Sites (double rates[],double alpha,int ncatG,int ls, int cdf,
    double space[])
{
/* Rates for sites from the gamma (ncatG=0) or discrete-gamma (ncatG>1).
   Rates are converted into the c.d.f. if cdf=1, which is useful for
   simulation under JC69-like models.

	. 	rates[] is a vector of size ls, and will have, in return the gamma
  		variates for sites.
	. 	ls means "length of sequence" in number of sites.
	. 	alpha is the parameter.  I fix the other parameter beta = alpha so that
  		the mean of the distribution is one.
	. 	ncatG is the number of categories in the discrete distribution.
  		ncatG = 0 means the continuous gamma distribution.  Use 8.
	. 	space[] is a vector of working space, with a minimum of size 3*ncatG
	. 	cdf stands for cumulative distribution function.

*/
   int h, ir,j, *counts=(int*)(space+2*ncatG);
   double *rK=space, *freqK=space+ncatG;

   if (alpha==0) 
      for (h=0; h<ls; h++) rates[h]=1;
   else {
      if (ncatG>1) {
         DiscreteGamma (freqK, rK, alpha, alpha, ncatG, 0);
         MultiNomial (ls, ncatG, freqK, counts, space+3*ncatG);
         for (ir=0,h=0; ir<ncatG; ir++) 
            for (j=0; j<counts[ir]; j++)  rates[h++]=rK[ir];
      }
      else 
         for (h=0; h<ls; h++) rates[h]=rndgamma(alpha)/alpha;
      if (cdf) {
         for (h=1; h<ls; h++) rates[h]+=rates[h-1];
         abyx (1/rates[ls-1], rates, ls);
      }
   }
   return (0);
}
//------------------------------------------------------------------------------
double LnGamma (double alpha)
{
/* returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.  
   Stirling's formula is used for the central polynomial part of the procedure.
   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
   double x=alpha, f=0, z;

   if (x<7) {
      f=1;  z=x-1;
      while (++z<7)  f*=z;
      x=z;   f=-log(f);
   }
   z = 1/(x*x);
   return  f + (x-0.5)*log(x) - x + .918938533204673 
          + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
               +.083333333333333)/x;  
}
//------------------------------------------------------------------------------
double DFGamma(double x, double alpha, double beta)
{
/* mean=alpha/beta; var=alpha/beta^2
*/
   if (alpha<=0 || beta<=0) cout << "err in DFGamma()\n";
   if (alpha>100) cout << "large alpha in DFGamma()\n";
   return pow(beta*x,alpha)/x * exp(-beta*x - LnGamma(alpha));

}
//------------------------------------------------------------------------------
double IncompleteGamma (double x, double alpha, double ln_gamma_alpha)
{
/* returns the incomplete gamma ratio I(x,alpha) where x is the upper 
           limit of the integration and alpha is the shape parameter.
   returns (-1) if in error
   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
   (1) series expansion     if (alpha>x || x<=1)
   (2) continued fraction   otherwise
   RATNEST FORTRAN by
   Bhattachrjee GP (1970) The incomplete gamma integral.  Applied Statistics,
   19: 285-287 (AS32)
*/
   int i;
   double p=alpha, g=ln_gamma_alpha;
   double accurate=1e-8, overflow=1e30;
   double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6];

   if (x==0) return (0);
   if (x<0 || p<=0) return (-1);

   factor=exp(p*log(x)-x-g);   
   if (x>1 && x>=p) goto l30;
   /* (1) series expansion */
   gin=1;  term=1;  rn=p;
 l20:
   rn++;
   term*=x/rn;   gin+=term;

   if (term > accurate) goto l20;
   gin*=factor/p;
   goto l50;
 l30:
   /* (2) continued fraction */
   a=1-p;   b=a+x+1;  term=0;
   pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
   gin=pn[2]/pn[3];
 l32:
   a++;  b+=2;  term++;   an=a*term;
   for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
   if (pn[5] == 0) goto l35;
   rn=pn[4]/pn[5];   dif=fabs(gin-rn);
   if (dif>accurate) goto l34;
   if (dif<=accurate*rn) goto l42;
 l34:
   gin=rn;
 l35:
   for (i=0; i<4; i++) pn[i]=pn[i+2];
   if (fabs(pn[4]) < overflow) goto l32;
   for (i=0; i<4; i++) pn[i]/=overflow;
   goto l32;
 l42:
   gin=1-factor*gin;

 l50:
   return (gin);
}
//------------------------------------------------------------------------------
double PointChi2 (double prob, double v)
{
/* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
   returns -1 if in error.   0.000002<prob<0.999998
   RATNEST FORTRAN by
       Best DJ & Roberts DE (1975) The percentage points of the 
       Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
   Converted into C by Ziheng Yang, Oct. 1993.
*/
   double e=.5e-6, aa=.6931471805, p=prob, g;
   double xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6;

   if (p<.000002 || p>.999998 || v<=0) return (-1);

   g = LnGamma (v/2);
   xx=v/2;   c=xx-1;
   if (v >= -1.24*log(p)) goto l1;

   ch=pow((p*xx*exp(g+xx*aa)), 1/xx);
   if (ch-e<0) return (ch);
   goto l4;
l1:
   if (v>.32) goto l3;
   ch=0.4;   a=log(1-p);
l2:
   q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch));
   t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
   ch-=(1-exp(a+g+.5*ch+c*aa)*p2/p1)/t;
   if (fabs(q/ch-1)-.01 <= 0) goto l4;
   else                       goto l2;
  
l3: 
   x=PointNormal (p);
   p1=0.222222/v;   ch=v*pow((x*sqrt(p1)+1-p1), 3.0);
   if (ch>2.2*v+6)  ch=-2*(log(1-p)-c*log(.5*ch)+g);
l4:
   q=ch;   p1=.5*ch;
   if ((t=IncompleteGamma (p1, xx, g))<0) {
      cout << "\nerr IncompleteGamma\n";
      return (-1);
   }
   p2=p-t;
   t=p2*exp(xx*aa+g+p1-c*log(ch));   
   b=t/ch;  a=0.5*t-b*c;

   s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
   s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520;
   s3=(210+a*(462+a*(707+932*a)))/2520;
   s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
   s5=(84+264*a+c*(175+606*a))/2520;
   s6=(120+c*(346+127*c))/5040;
   ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
   if (fabs(q/ch-1) > e) goto l4;

   return (ch);
}
//------------------------------------------------------------------------------
int DiscreteGamma (double freqK[], double rK[],
    double alfa, double beta, int K, int median)
{
/* discretization of gamma distribution with equal proportions in each 
   category
*/
   int i;
   double gap05=1.0/(2.0*K), t, factor=alfa/beta*K, lnga1;

   if (median) {
      for (i=0; i<K; i++) rK[i]=PointGamma((i*2.0+1)*gap05, alfa, beta);
      for (i=0,t=0; i<K; i++) t+=rK[i];
      for (i=0; i<K; i++)     rK[i]*=factor/t;
   }
   else {
      lnga1=LnGamma(alfa+1);
      for (i=0; i<K-1; i++)
         freqK[i]=PointGamma((i+1.0)/K, alfa, beta);
      for (i=0; i<K-1; i++)
         freqK[i]=IncompleteGamma(freqK[i]*beta, alfa+1, lnga1);
      rK[0] = freqK[0]*factor;
      rK[K-1] = (1-freqK[K-2])*factor;
      for (i=1; i<K-1; i++)  rK[i] = (freqK[i]-freqK[i-1])*factor;
   }
   for (i=0; i<K; i++) freqK[i]=1.0/K;

   return (0);
}
//------------------------------------------------------------------------------
double PointNormal (double prob)
{
/* returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
   returns (-9999) if in error
   Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
   Applied Statistics 22: 96-97 (AS70)

   Newer methods:
     Wichura MJ (1988) Algorithm AS 241: the percentage points of the
       normal distribution.  37: 477-484.
     Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage 
       points of the normal distribution.  26: 118-121.

*/
   double a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245;
   double a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495;
   double b2=.531103462366, b3=.103537752850, b4=.0038560700634;
   double y, z=0, p=prob, p1;

   p1 = (p<0.5 ? p : 1-p);
   if (p1<1e-20) return (-9999);

   y = sqrt (log(1/(p1*p1)));   
   z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
   return (p<0.5 ? -z : z);
}
//------------------------------------------------------------------------------
double rndgamma1 (double s);
double rndgamma2 (double s);

double rndgamma (double s)
{
/* random standard gamma (Mean=Var=s,  with shape par=s, scale par=1)
      r^(s-1)*exp(-r)
   J. Dagpunar (1988) Principles of random variate generation,
   Clarendon Press, Oxford
   calling rndgamma1() if s<1 or
           rndgamma2() if s>1 or
           exponential if s=1 
*/

   double r=0;

   if (s<=0)      cout << "jgl gamma..\n";
   else if (s<1)  r=rndgamma1 (s);
   else if (s>1)  r=rndgamma2 (s);
   else           r=-log(rndu());
   return (r);
}
//------------------------------------------------------------------------------
double rndgamma1 (double s)
{
/* random standard gamma for s<1
   switching method
*/
   double r, x=0,sml=1e-37,w;
   static double a,p,uf,ss=10,d;

   if (s!=ss) {
      a=1-s;
      p=a/(a+s*exp(-a));
      uf=p*pow(sml/a,s);
      d=a*log(a);
      ss=s;
   }
   for (;;) {
      r=rndu();
      if (r>p)        x=a-log((1-r)/(1-p)), w=a*log(x)-d;
      else if (r>uf)  x=a*pow(r/p,1/s), w=x;
      else            return (0);
      r=rndu ();
      if (1-r<=w && r>0)
         if (r*(w+1)>=1 || -log(r)<=w)  continue;
      break;
   }
   return (x);
}
//------------------------------------------------------------------------------
double rndgamma2 (double s)
{
/* random standard gamma for s>1
   Best's (1978) t distribution method
*/
   double r,d,f,g,x;
   static double b,h,ss=0;
   if (s!=ss) {
      b=s-1;
      h=sqrt(3*s-0.75);
      ss=s;
   }
   for (;;) {
      r=rndu ();
      g=r-r*r;
      f=(r-0.5)*h/sqrt(g);
      x=b+f;
      if (x <= 0) continue;
      r=rndu();
      d=64*r*r*g*g*g;
      if (d*x < x-2*f*f || log(d) < 2*(b*log(x/b)-f))  break;
   }
   return (x);
}
//------------------------------------------------------------------------------
int xtoy (double x[], double y[], int n)
{ int i; for (i=0; i<n; y[i]=x[i],i++) ;  return(0); }
//------------------------------------------------------------------------------
int abyx (double a, double x[], int n)
{ int i; for (i=0; i<n; x[i]*=a,i++) ;  return(0); }
//------------------------------------------------------------------------------
static int z_rndu=137;
static unsigned w_rndu=13757;
//------------------------------------------------------------------------------
void SetSeed (int seed)
{
   z_rndu = 170*(seed%178) + 137;
   w_rndu=seed;
}
//------------------------------------------------------------------------------
double rndu (void)
{
/* U(0,1): AS 183: Appl. Stat. 31:188-190
   Wichmann BA & Hill ID.  1982.  An efficient and portable
   pseudo-random number generator.  Appl. Stat. 31:188-190

   x, y, z are any numbers in the range 1-30000.  Integer operation up
   to 30323 required.
*/
   static int x_rndu=11, y_rndu=23;
   double r;

   x_rndu = 171*(x_rndu%177) -  2*(x_rndu/177);
   y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176);
   z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178);
   if (x_rndu<0) x_rndu+=30269;
   if (y_rndu<0) y_rndu+=30307;
   if (z_rndu<0) z_rndu+=30323;
   r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0;
   return (r-(int)r);
}
//------------------------------------------------------------------------------
int MultiNomial (int n, int ncat, double prob[], int nobs[], double space[])
{
/* sample n times from a mutinomial distribution M(ncat, prob[])
   prob[] is considered cumulative prob if (space==NULL)
*/
   int i, j, crude=(ncat>20), ncrude=5, lcrude[5];
   double r, *pcdf=(space==NULL?prob:space);

   FOR (i, ncat) nobs[i]=0;
   if (space) {
      xtoy (prob, pcdf, ncat);
      for (i=1; i<ncat; i++) pcdf[i]+=pcdf[i-1];
   }
   if (fabs(pcdf[ncat-1]-1) > 1e-5) cout << "sum P!=1 in MultiNomial\n";
   if (crude) {
      for (j=1,lcrude[0]=i=0; j<ncrude; j++)  {
         while (pcdf[i]<(double)j/ncrude) i++;
         lcrude[j]=i-1;
      }
   }
   FOR (i, n) {
      r=rndu();
      j=0;
      if (crude) {
         for (; j<ncrude; j++) if (r<(j+1.)/ncrude) break;
         j=lcrude[j];
      }
      for (; j<ncat; j++) if (r<pcdf[j]) break;
      nobs[j] ++;
   }
   return (0);
};
//------------------------------------------------------------------------------
