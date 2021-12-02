#ifndef _genealogy_h_              // Sentry, use file only if it's not already included.
#define _genealogy_h_

#include <fstream>
#include <stdlib.h>
#include "arrays.h"
#include "mutation.h"
#include "algebra.h"

enum tree_output_type {GENERATIONS=0, MUT_RATE=1, NUM_MUT=2};
typedef MY_TArrayAsVector<int> TIntVect;

//Declaration of two classes defined elsewhere but needed here
class TDemeCollection;
class TMigrationMatrix;

class TNode {
	friend ostream& operator<<(ostream& os, const TNode& node); //Christian: Why not here already?
	public:   
   	long time;
      TNode *desc1;
      TNode *desc2;
      TNode *ancestor;
      TIntVect* desc1_nodes;
      TIntVect* desc2_nodes;
      TIntVect* sequence;     //The sequence to mutate
      static TIntVect hits;   //The number of hits per site
      static GammaRates mut_rates; //A static vector of gamma distributed rates
      int seq_length;
      static double mut_rate;			//Global mutation rate along the tree
      long num_new_mut;       //Number of mutations as compared to the direct ancestor
      long left_mut;
      long right_mut;
      int node_number;
      int deme;        		//Remember to which deme it belongs
      int stat_group;		//...and which stat group
      static double tree_length; //Length of the tree to which the node belongs, expressed in number of mutations
      static int min_mic, max_mic; //Range constraints for microsatellite data //Loro_26_2_99
  	   TNode() {
         time=0;
         desc1=desc2=ancestor=NULL;
         desc1_nodes=desc2_nodes=NULL;
         num_desc1=num_desc2=0;
         node_count++;
         node_number=node_count;
         deme=-1;
         num_new_mut=left_mut=right_mut=0;
         seq_length=0;
         sequence=NULL;
      }      
      TNode(const TNode& N)   { if (this!=&N) *this=N; }
      ~TNode() {
      	if (desc1_nodes) {
            desc1_nodes->Flush();
         	delete desc1_nodes;
            desc1_nodes=NULL;
        }
   		if (desc2_nodes) {
            desc2_nodes->Flush();
         	delete desc2_nodes;
            desc2_nodes=NULL;
        }
   		if (sequence) {
            sequence->Flush();
         	delete sequence;
            sequence=NULL;
        }
      }
      int operator==(const TNode& N) {
      	 if (desc1==N.desc1 && desc2==N.desc2 &&
             ancestor==N.ancestor)  return 1;
         return 0;
      }
      TNode& operator=(const TNode& N);
      
      void set(int d, int t, int sg) {deme=d; time=t; stat_group=sg;}
      void set_mut_rate(const double& m) {mut_rate=m;}
      double tree_exp_length() {return tree_length;}
      void reset_tree_length() {tree_length=0;}
      int count_desc();
      int count_polym_sites();
      long get_time() {return time;};
      int build_lists_of_descendent_nodes();
   	int compute_total_coal_times_among_demes(	TDemeCollection& DemeCollec,
      														TMigrationMatrix& CoalTimes,
      														TMigrationMatrix& PairDiff,
   													  		TMigrationMatrix& MinCoalTimes);
      int print_info(ostream& os);
   	int num_desc1, num_desc2;
      void reset_node_count() {node_count=0;};
		long add_mutations(int num_mut, int len, Mut_Type mut_type, float gamma_par, float trans_rate, int range_const,symb *muteq);
		//A recursive routine to print tree topology and branch length
      void print_tree_structure(ostream& os, const tree_output_type tree_type, const float & mu);
      //Christian: One I like better
      void TempRecNodeOut(int lvl);

   protected:
   	static int node_count;
};

ostream & operator<<(ostream& os, const TNode& node);

//An array of pointers to Nodes
//It is like an indirect array of Nodes where Nodes would not be deleted
//in the destructor of the indirect array
class TListNode {
	private:
   	TNode **list;
      int _size;
	public:
		TListNode() {list=NULL; _size=0;}
		TListNode(const int size, TNode * tree);
		~TListNode() {
			if (list) delete[] list;
      }
      TNode *& operator[](const int n) const {return list[n];}
      int allocate_list(const int size);
      void de_allocate_list() {
      	if (list) delete[] list;
         list=NULL;
      }
      int size() {return _size;}
};

class TTree {
//Christian: a debugging routine to print the tree
   	friend ostream& operator<<(ostream& os, const TTree& T);

	public:
      TTree() {
      	 SampleSize=0;
      	 tree=NULL;
         delete_nodes=1;
         _mean_coalescence_time=0;
         _sd_coalescence_time=0;
      };
      TTree(int size);
      ~TTree() {
      	if (tree) {
      		if (delete_nodes) delete[] tree;
         	else delete tree;
         }
      }
      TNode& operator[](const int n) const {return tree[n];}
      TNode* operator()(const int i) const {return tree + i;}
      int allocate_tree(const int& sampsize);
      void print_nodes(ostream& os, const tree_output_type tree_type, const float & mu);
      void print_sequences(ostream& os, const Mut_Type& mut_type);
      int sample_size() const {return SampleSize;}
      void set_delete_nodes() {delete_nodes=1;}
      void set_not_delete_nodes() {delete_nodes=0;}
      float total_length() {if (tree) return tree[0].tree_length; else return 0;};
      const double& mean_coalescence_time() {return _mean_coalescence_time;}
      const double& sd_coalescence_time() {return _sd_coalescence_time;}
      void count_descendants() { tree[2*SampleSize-2].count_desc(); } //Christian: Should be in the code already?
        
   private:
      int SampleSize;
      TNode* tree;
      int delete_nodes;
      double _mean_coalescence_time;
      double _sd_coalescence_time;
      double _tree_subst_length;
      double _num_mut_on_tree;
};

extern ostream& operator<<(ostream& os, const TTree& T);

#endif  // genealogy_h sentry.

