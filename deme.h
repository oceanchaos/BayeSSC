#ifndef __DEME_H__
#define __DEME_H__

//#include <classlib\arrays.h>
#include "arrays.h"
#include "public.h"
#include "migrmat.h"
#include "genealgy.h"
#include "genstat.h"

#define MAXLONG	0x7fffffff

const long MAX_GENERATIONS=100000000L;
const long MIN_SIZE=1;
const long MAX_SIZE=2000000000L;

//------------------------------------------------------------------------------
//Christian 06/28/04: The TSamp class. If samples are to be taken from
//the same deme at multiple time points, sample size cannot be treated as a
//simple <long>
class TSamp {
	friend ostream& operator<<(ostream& os, const TSamp& S);

    public :
      long _size;  //number of samples
      long _age;   //number of generations old these samples are
      int _deme;   //the deme this sample comes from
      int _stat_grp; //the stat group of the sample

      TSamp () {;}
      void SetVals(const long& sz, const long& ag, int d, int sg) {
         _size=sz; _age=ag; _deme=d; _stat_grp=sg;
      }
      ~TSamp() {;}
      bool operator< (TSamp t) {
         if (_age<t._age) return 1;
         return 0;
      }
};

extern ostream& operator<<(ostream& os, const TSamp& S); 

typedef MY_TSArrayAsVector<TSamp>		TSampArray;  //Christian 6/29/04

//------------------------------------------------------------------------------
class TDeme {

	friend ostream& operator<<(ostream& os , const TDeme& D);

	private :
      static TMigrationMatrix* Migrations; //A pointer to a migration matrix
      double  _deme_size;
      long  _min_deme_size, _max_deme_size;
      float _proba_coal;
      int _id;
      bool _is_migration;
      long _polym_sites;
   public :
      double _growth_rate;
      TSampArray _samples;  //Christian 6/29/04
      long  _lineages;
      long _oldest_samp;
      long _next_samp;
      TListNode OriginalNodes;  //A list of the nodes assigned to the deme (past or present)
      TListNode NodeList;    	//extant nodes only

   	TDeme() : NodeList(), OriginalNodes() {
         initialize();
      };
      TDeme(const TDeme& D)  : NodeList(), OriginalNodes() {
      	if (this!=&D) *this=D;
      };
   	TDeme(const long& ds, const TSampArray samps, const float& gr) : NodeList(), OriginalNodes()  {
         _growth_rate=gr;
         _samples=samps;
         _deme_size=(double)ds;
         _lineages=_samples.count();
         _is_migration=true; 
         Migrations=NULL;
       	_polym_sites=0L;
       	_oldest_samp=0L;
      };
      ~TDeme() {;};
      void set_growth_rate(const float& gr) { _growth_rate=gr; }
      void AddSample(const TSamp s) {
      	 _samples.Add(s);
         _lineages+=s._size;
      }
      void DumpSamples() { _samples.Flush(); _lineages=0; }
      void set_size_bounds(const long& min, const long& max) {
			_min_deme_size=min;
			_max_deme_size=max;
      }
      void set_deme_size(const long& size) {
      	_deme_size=(double)size;
   		if (size) _proba_coal=1.0/_deme_size;
         else _proba_coal=0.0;
      }
      void set_id(const int id) {_id=id;}
      void initialize() {
      	 _growth_rate=0;  //Stationary population
         _deme_size=1.0;
         _min_deme_size=MIN_SIZE;
         _max_deme_size=MAX_SIZE;
         _lineages=0;
         _is_migration=true;
         Migrations=NULL; 
       	_polym_sites=0L;
       	_oldest_samp=0L;
      }
      void set_migration_matrix(TMigrationMatrix* MM) { Migrations=MM; }
      TDeme operator=(TDeme D);
      int operator==(TDeme D) {
      	if (_id==D._id) return 1;
         return 0;
      };
      int operator<(const TDeme& D) { return _id<D._id; }
      int copy_node_pointers(const TTree& T, int curr_lin, const long& max_lin);
      //Implementation of a single migration event towards another deme (sink)
      int migrate(TDeme& sink, int x) {
   		//Add selected node to sink and then increments sink's number of lineages
   		sink.NodeList[sink._lineages++]=NodeList[x];
   		//Decrements last node's position, and replace selected node by last node
   		//CHRISTIAN: BUG!!! What if the migrating node *IS* the last node?
   		if(x!=(--_lineages)) NodeList[x]=NodeList[_lineages];
         //printf("\n-Moved #%d from %d to %d",x,_id,sink._id);
         return 1;
      };

      //Check and implements migrations
      int send_migrants(TDeme& sink, const int& i,const long& time);

      //Implements the coalescence process within population
      int coalesce_lineages(const long& time, TTree& GeneTree, int pos, bool hudson);
      float growth() {return (_growth_rate!=0.0);}

      //Implements an exponential growth of population for positive _growth_rate
      //and and exponential decline for negative _growth_rate values
      float exponential_resize() {
         _deme_size*=exp(_growth_rate);
         if ((long)_deme_size<_min_deme_size) {
         	_deme_size=_min_deme_size;
			   _growth_rate=0.0;
         }
         if ((long)_deme_size>_max_deme_size) {
         	_deme_size=_max_deme_size;
			   _growth_rate=0.0;
         }
   		if (_deme_size) _proba_coal=1.0/_deme_size;
         else _proba_coal=0.0;
         return _deme_size;
      }

      //Implements a resizing of the deme size according to a given factor
      void linear_resize(const float& factor) {
      	_deme_size*=factor;
   		if (_deme_size) _proba_coal=1.0/_deme_size;
         else _proba_coal=0.0;
      }
      void reset();
      long sample_size() {return _samples.count();}
      int  check_migrations();
      bool is_migration() {return _is_migration;}
      long deme_size() {return (long)_deme_size;}
      void print_haplotypes_for_Arlequin(ostream & os, const int& n, const Mut_Type data_type);
      void print_haplotypes_for_PAUP(ostream & os, const Mut_Type data_type);
};

extern ostream& operator<<(ostream& os , const TDeme& D);


//------------------------------------------------------------------------------
// A class specifying the conditions of a demographic historical change
// sometimes in the past, like the fusion of two demes, or a change in
// deme size

class THistoricalEvent {
	friend ostream& operator<<(ostream& os, const THistoricalEvent& HE);
	public:
	   long time;              //time in the past at which the event occurs (in generations)
      int source;    			//deme which sends of migrants
      int sink;      			//deme which receives migrants
      float migrants; 			//Proportion of migrants sent to sink from source
      float new_deme_size;    //New relative deme size for SOURCE deme
      float new_growth_rate; 	//New growth rate for SOURCE deme
      int MigMat;             //New migration matrix number
      THistoricalEvent() {};
      ~THistoricalEvent() {};
      int operator==(const THistoricalEvent& HE) {
      	if (source==HE.source &&
     		sink==HE.sink &&
            time==HE.time &&
            migrants==HE.migrants &&
            new_deme_size==HE.new_deme_size &&
            new_growth_rate==HE.new_growth_rate &&
            MigMat==HE.MigMat) return 1;
         return 0;
      }
      int operator<(const THistoricalEvent& HE) const {
      	return time<HE.time;
      }
};

extern ostream& operator<<(ostream& os, const THistoricalEvent& HE); 

//------------------------------------------------------------------------------
typedef MY_TArrayAsVector<TDeme> 			    TDemeArray;
typedef MY_TSArrayAsVector<THistoricalEvent>    TEventArray;
typedef MY_TArrayAsVector<long>                 TDemeSizeArray;
typedef MY_TArrayAsVector<float>                TGrowthArray;
typedef MY_TArrayAsVector<TMigrationMatrix>     TMigrMatArray; //Loro_16_9_99

class TDemeCollection {
	friend ostream& operator<<(ostream& os,const TDemeCollection& DC);
  	private :
      TDeme          *Demes;
      TEventArray 	*Events;
      TMigrMatArray 	*MigrMatArray;
      TMigrationMatrix 	*CoalTimes; //A matrix of mean pairwise coalescence times
      TMigrationMatrix 	*PairDiff; //A matrix of mean number of pairwise differences
      TMigrationMatrix 	*MinCoalTimes; //A matrix of minimum pairwise coalescence times
      TTree			   GeneTree;
      bool _is_migration;
      int _num_demes; 
      long _polym_sites;
      long _tot_mut;
      int _cur_mig_mat;
      int _num_stat_grps; //Christian
   public :
   	TDemeCollection() : GeneTree() {
   	   Events=NULL;
      	CoalTimes=NULL;
         MinCoalTimes=NULL;
         PairDiff=NULL;
         Demes=NULL;
         _num_demes=0; 
       	_polym_sites=0L;
       	_tot_mut=0L;
         _cur_mig_mat=0;
         _num_stat_grps=0;
      };
      ~TDemeCollection() {
      	if (CoalTimes) delete CoalTimes;
      	if (PairDiff) delete PairDiff;
      	if (MinCoalTimes) delete MinCoalTimes;
         if (Demes) delete[] Demes;
      }
      TDemeCollection& operator=(const TDemeCollection& DC);
      long count_lineages();
      int reajust_deme_size(const int i, const long& ds);
      int stat_grps() {return _num_stat_grps; } //Christian
      void set_num_demes(int n) {_num_demes=n;}
      int create_demes(TMigrMatArray *MA, TDemeSizeArray *DS, TGrowthArray *GR, TEventArray *HE, TSampArray *SS);
      int build_tree(ostream &logfile, int simnum);
      TMigrationMatrix& migr_matrix() {return (*MigrMatArray)[_cur_mig_mat];}
      MatElemType& min_coal_time(const int i, const int j) const {
      	return (*MinCoalTimes)(i,j);
      }
      MatElemType& mean_pair_diff(const int i, const int j) const {
      	return (*PairDiff)(i,j);
      }
      void print_gene_tree(ostream& os, const tree_output_type tree_type, const float & mu) {
      	GeneTree.print_nodes(os, tree_type, mu);
      }
      void print_sequences(ostream& os, const Mut_Type& mut_type) {
      	GeneTree.print_sequences(os, mut_type);
      };
      int implement_event(const int& cur_event);
      void reset(const TDemeSizeArray& SA, const TGrowthArray& GA, TSampArray SS);
      double mean_coal_time() {return GeneTree.mean_coalescence_time();}
      double sd_coal_time() {return GeneTree.sd_coalescence_time();}

      const TNode& tree_node(const int& pos) {return GeneTree[pos];}
      TDeme& operator[](const int& d) {return Demes[d];}
      const TMigrationMatrix& coal_time_mat() const {return *CoalTimes;}
      const TMigrationMatrix& mean_pair_diff_mat() const {return *PairDiff;}  
      bool is_migration() {return _is_migration;}
      int check_migrations();
      int check_growth();
      int num_demes() const {return _num_demes;} 
      const long& sprinkle_mutations(const float& mut_rate, const int& len,
				const Mut_Type& mut_type, const float& gamma_par,
                const float& mut_ratio, const float& prop_sites,
                const float& trans_rate, const int& range_constr, symb *muteq);
      void write_samples_to_Arlequin_file(ostream& os, const Mut_Type data_type);
      void write_group_section_to_Arlequin_file(ostream& os);
      void write_samples_to_PAUP_file(ostream& os, const Mut_Type data_type);
      void write_Statfile(ostream& os, const Mut_Type data_type);
};
extern ostream& operator<<(ostream& os,const TDemeCollection& DC);
#endif
