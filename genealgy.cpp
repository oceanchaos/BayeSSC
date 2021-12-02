#include "public.h"
#include "genealgy.h"
#include "deme.h"

enum UNITS {WIDTH=6};

//Initialization of static members of TNode
TIntVect 	TNode::hits(0,10);
GammaRates 	TNode::mut_rates(0,0.0,1.0,0.0, 0);
double 		TNode::mut_rate=0.0;
double 		TNode::tree_length=0.0;
int 			TNode::node_count=0;
int 			TNode::min_mic=0;
int 			TNode::max_mic=0; 

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
TNode& TNode::operator=(const TNode& N) {
	time=N.time;
   desc1=N.desc1;
   desc2=N.desc2;
   ancestor=N.ancestor;
   num_desc1=N.num_desc1;
   num_desc2=N.num_desc2;
   node_number=N.node_number;
   deme=N.deme;
   num_new_mut=N.num_new_mut;
   left_mut=N.left_mut;
   right_mut=N.right_mut;
   seq_length=N.seq_length;
   if (desc1_nodes) 	delete desc1_nodes; 	desc1_nodes=NULL;
   if (desc2_nodes) 	delete desc2_nodes; 	desc2_nodes=NULL;
   if (sequence) 		delete sequence; 		sequence=NULL;
   if (N.desc1_nodes) {
   	try {
      	int size=N.desc1_nodes->GetItemsInContainer();
      	desc1_nodes= new TIntVect(size,10);
         for (int i=0; i<size; ++i) {
         	desc1_nodes[i]=N.desc1_nodes[i];
         }
      }
      catch (...) {
   		if (desc1_nodes) 	delete desc1_nodes; 	desc1_nodes=NULL;
         cout << "\nTNode::operator=(): unable to allocate memory\n";
         return *this;
      }
   }
   
   if (N.desc2_nodes) {
   	try {
      	int size=N.desc2_nodes->GetItemsInContainer();
      	desc2_nodes= new TIntVect(size,10);
         for (int i=0; i<size; ++i) {
         	desc2_nodes[i]=N.desc2_nodes[i];
         }
      }
      catch (...) {
   		if (desc2_nodes) 	delete desc2_nodes; 	desc2_nodes=NULL;
         cout << "\nTNode::operator=(): unable to allocate memory\n";
         return *this;
      }
   }

   if (N.sequence) {
   	try {
      	int size=N.sequence->GetItemsInContainer();
      	sequence= new TIntVect(size,10);
         for (int i=0; i<size; ++i) {
         	sequence[i]=N.sequence[i];
         }
      }
      catch (...) {
   		if (sequence) 	delete sequence; 	sequence=NULL;
         cout << "\nTNode::operator=(): unable to allocate memory\n";
         return *this;
      }
   }
   return *this;
}
//----------------------------------------------------------------------------
// A recursive method to count the number of nodes below a given node,
// including himself
int TNode::count_desc() {
	num_desc1=num_desc2=0;
	if (desc1) num_desc1=desc1->count_desc();
	if (desc2) num_desc2=desc2->count_desc();
   //Counts only the basal tips
   if (!desc1 && !desc2) return num_desc1+num_desc2+1;
   return num_desc1+num_desc2;
};
//----------------------------------------------------------------------------
int TNode::count_polym_sites() {
	int polym_loci=0, size=hits.GetItemsInContainer();
	for (int i=0; i<size; ++i)
   	if (hits[i]) ++polym_loci;
   return polym_loci;
};
//----------------------------------------------------------------------------
//A recursion procedure to establish the list of descendent nodes to the left
//and to the right of the current node
int
TNode::build_lists_of_descendent_nodes() {

	//Creates empty list of descendent nodes

		try {
   		if (!desc1_nodes) desc1_nodes= new TIntVect(1,10);
         else desc1_nodes->Flush();
      	if (!desc2_nodes) desc2_nodes= new TIntVect(1,10);
         else desc2_nodes->Flush();
   	}
   	catch (...) {
   		if (desc1_nodes) delete desc1_nodes; desc1_nodes=NULL;
      	if (desc2_nodes) delete desc2_nodes; desc2_nodes=NULL;
      	cout 	<< "TNode:: build_list_of_descendent_nodes() : Unable to allocate memory"
      			<< endl;
      	return 0;
   	}

   //Begin recursion
   if (desc1) {
   	//Adds the descendent node list of descent nodes to the current one
   	if (!desc1->build_lists_of_descendent_nodes()) return 0;
      int num_nodes=desc1->desc1_nodes->GetItemsInContainer();
      for (int i=0; i<num_nodes; ++i) {
      	desc1_nodes->Add((*desc1->desc1_nodes)[i]);
      }
      num_nodes=desc1->desc2_nodes->GetItemsInContainer();
      for (int i=0; i<num_nodes; ++i) {
      	desc1_nodes->Add((*desc1->desc2_nodes)[i]);
      }
   }
   else {
   	//Adds its id in its own list of descendent nodes !
   	desc1_nodes->Add(node_number);
   }
   if (desc2) {
   	//Adds the descendent node list of descent nodes to the current one
   	if (!desc2->build_lists_of_descendent_nodes()) return 0;
      int num_nodes=desc2->desc1_nodes->GetItemsInContainer();
      for (int i=0; i<num_nodes; ++i) {
      	desc2_nodes->Add((*desc2->desc1_nodes)[i]);
      }
      num_nodes=desc2->desc2_nodes->GetItemsInContainer();
      for (int i=0; i<num_nodes; ++i) {
      	desc2_nodes->Add((*desc2->desc2_nodes)[i]);
      }
   }
   return 1;
};
//----------------------------------------------------------------------------
ostream & operator<<(ostream& os, const TNode& node) {
   bool showseq=0;
	if(node.desc1 || node.desc2) {
		os << "\nInternal node No " << node.node_number
	     	<< "\n----------------"
    	   << "\nTime since present       : " << node.time
	   	<< "\nNo. of left descendents  : " << node.num_desc1
         << "\nNo. of right descendents : " << node.num_desc2
	      	/*<< "\nLeft descendent address  : " << node.desc1
         << "\nRight descendent address : " << node.desc2
         << "\nAncestor address         : " << node.ancestor*/ << endl;
	} else
		os << "\nTerminal node No " << node.node_number << " (Age="
  		   << node.time << ", Deme=" << node.deme << ", StatGroup="
           << node.stat_group << ")";
   if(showseq) {
     	if(!node.sequence) printf("NO SEQUENCE!");
     	else {
			int *c=node.sequence->begin();
     		for(int i=0;i<40;i++) os<<*c++;
   	}
	}
   return os;
};

// Print node information recursively
int TNode::print_info(ostream& os) {
	os << *this;
	if (desc1) desc1->print_info(os);
	if (desc2) desc2->print_info(os);
    return 1;
};

//A recursive procedure to add mutations to the tree starting from the root
long TNode::add_mutations(int num_mut, int len, Mut_Type mut_type, float gamma_par, float trans_rate, int range_const, symb *muteq) {
	//cout << *this;
   long tot_mut=0L, desc_mut;
   double drand_num;
   float subst_length;
   seq_length=len;
   //Create new sequence and propagate the ancestor's mutations WHAT?!?!?!
   if (!sequence)  {
   	try { sequence= new TIntVect(len, 0);
   	} catch (...) {
   		if (sequence) delete sequence; sequence=NULL;
      	seq_length=0;
      	return 0;
   	}
   }
   if (ancestor) {
       	int *cursite=sequence->begin(), *cursite_anc=ancestor->sequence->begin();
       	for (int i=0; i<len; ++i) *cursite++=*cursite_anc++; //Inheritance of ancestor's mutations
   } else {    //Loro_15_9_98 Create the ancestor's genetic data
       	int* cursite=sequence->begin();
       	if (mut_type==DNA) for (int i=0; i<len; ++i) *cursite++=int(4*mtr()); // 0-3
        	else for (int i=0; i<len; ++i) *cursite++=0;
   }

   //This is the number of new mutations to generate as compared to parent node
   num_new_mut=num_mut;

   //Generate those mutations
   for (int i=0, pos; i<num_mut; ++i) {
		//Find the position of the site to be hit by a mutation
      if (fabs(gamma_par)<1e-7) { //Close to zero, even mutation rates are assumed among sites
        	pos=int(mtr()*len);
      } else {
      	drand_num=mtr();
      	double tot_prob=0.0;
			//Explore from the end of the sequence because it is where mutations are
			//more likely to occur
			for (pos=len-1; pos>-1; --pos) {
				tot_prob+=mut_rates[pos];
				if (drand_num<tot_prob) break;
			}
      }
      ++hits[pos];
		int& site=(*sequence)[pos];
		switch (mut_type) {
		  	case MICROSAT:	//Loro_26_2_99
			  	if (range_const) { //There are bouncing walls at min_mic and max_mic
            	if (site==min_mic) ++site;
               else if (site==max_mic) --site;
               else (mtr()<0.5) ? site-- : site++; //otherwise they are free to move randomly
				} else (mtr()<0.5) ? site-- : site++;
		 		break;
         case DNA:  drand_num=mtr();
				//Here we implement transition bias- 01: AG;   23: CT
            switch (site) {
            case 0: 	if (drand_num<trans_rate) site=1;
            			else site=(mtr()<0.5) ? 2:3;
                     break;
            case 1: if (drand_num<trans_rate) site=0;
            			else site=(mtr()<0.5) ? 2:3;
                     break;
            case 2: if (drand_num<trans_rate) site=3;
            			else site=(mtr()<0.5) ? 0:1;
                     break;
            case 3: if (drand_num<trans_rate) site=2;
            			else site=(mtr()<0.5) ? 0:1;
                     break;
            }
				break;
         case RFLP: site=!site; break;
         default: break;
      }
   }

   //Continue recursion
   if (desc1) {
        desc_mut=0;
		subst_length=mut_rate*(time-desc1->time);
        if(muteq) for(int t=time; t>desc1->time; t--) {
			double td=(double)t;
			//printf("\nt=%f mu=%f",td,muteq->solve(&td));
			//system("PAUSE");
			if(mtr()<muteq->solve(&td)) desc_mut++;
		} else desc_mut= (int) poidev(subst_length);
      tree_length+=subst_length;
      tot_mut+=desc_mut+desc1->add_mutations(desc_mut, len, mut_type,gamma_par,trans_rate, range_const,muteq);
   }
   if (desc2) {
       desc_mut=0;
       subst_length=mut_rate*(time-desc2->time);
       if(muteq) for(int t=time; t>desc1->time; t--) {
			double td=(double)t;
			if(mtr()<muteq->solve(&td)) desc_mut++;
		} else desc_mut= (int) poidev(subst_length);
       tree_length+=subst_length;
       tot_mut+=desc_mut+desc2->add_mutations(desc_mut, len, mut_type, gamma_par,trans_rate, range_const,muteq);
   }
   return tot_mut;
};
//----------------------------------------------------------------------------
//A recursive routine to print tree topology and branch length
//Loro 27.8.98
void
TNode::print_tree_structure(ostream& os, const tree_output_type tree_type, const float & mu) {
	if (desc1 && desc2) os << "(";
   if (desc1) desc1->print_tree_structure(os, tree_type, mu);
   if (desc1 && desc2) {
   	os << ", ";
   }
   else {
   	os << (deme+1) << "." << node_number << "." << stat_group;
   }
   if (desc2) desc2->print_tree_structure(os, tree_type, mu);
   if (desc1 && desc2) os << ")";
   if (ancestor) {
   	switch (tree_type)  {
      	case GENERATIONS : 	os << ":" << (ancestor->time-time);
                              break;
      	case MUT_RATE    : 	os << ":" << (ancestor->time-time)*mu;
                              break;
      	case NUM_MUT     : 	os << ":" << num_new_mut;
                              break;
      }
   }
   else os << ";\n";
};
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
int
TListNode::allocate_list(const int size) {
	if (list) delete[] list;
   try {
		list = new TNode*[size];
      _size=size;
   }
   catch (...) {
   	if (list) delete[] list;
      list=NULL;
      _size=0;
      cout << "TListNode::allocate_list : unable to allocate memory" << endl;
      return 0;
   }
   return 1;
}
//----------------------------------------------------------------------------
TListNode::TListNode(const int size, TNode * tree) {
	if (tree)
	try {
		list = new TNode*[size];
      for (int i=0; i<size; ++i) {
      	tree[i].time=0;
      	tree[i].desc1=tree[i].desc2=NULL;
         list[i]=tree + i;
      }
      _size=size;
   }
   catch (...) {
   	if (list) delete[] list;
      list=NULL;
      _size=0;
      cout << "TListNode::TListNode(const int size, TNode * tree) : unable to allocate memory" << endl;
      return;
   }
};
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

TTree::TTree(int size) {
	SampleSize=size;
   delete_nodes=1;
   _mean_coalescence_time=0;
   _sd_coalescence_time=0;
   try {
   	tree=new TNode[2*SampleSize];
   }
   catch (...) {
   	if (tree) delete[] tree; tree=NULL;
      cout << "TTree::TTree(int size) : unable to allocate memory" << endl; 
   }
};
//----------------------------------------------------------------------------
int
TTree::allocate_tree(const int& sampsize) {                           
	SampleSize=sampsize;
	int array_size=2*SampleSize-1;
   //Create tree structure if needed
   if (tree) delete[] tree;
   try {
   	tree=new TNode[array_size];
   }
   catch (...) {
   	if (tree) delete[] tree;
      return 0;
   }
   return 1;
};

//----------------------------------------------------------------------------
void TTree::print_nodes(ostream& os, const tree_output_type tree_type, const float & mu) {
	//Root the tree
	TNode& root=tree[2*SampleSize-2];
   root.print_tree_structure(os, tree_type, mu);
}

//----------------------------------------------------------------------------
void TTree::print_sequences(ostream& os, const Mut_Type& mut_type) {
   int len=tree[0].seq_length;
   char DNA_letters[4]={'A','G','C','T'};
   if (len) {
   	//Prints first the number of hits per site
   	os << "Hits\t";
		for (int j=0; j<len; ++j) {
      	os << tree[0].hits[j] << " ";
      }
      os << "\n\n";

      //Then prints the sequences
      for (int i=0; i<SampleSize; ++i) {
   		os << (i+1) << "\t\t";
      	if (tree[i].sequence) {
         	int* cursite=tree[i].sequence->begin();
   			for (int j=0; j<len; ++j) {
            	if (mut_type==DNA)
   					//os << DNA_letters[(*tree[i].sequence)[j]];
                  os << DNA_letters[*cursite++];
               else
   					//os << (*tree[i].sequence)[j];
                  //Loro_1_3_99 : In order to avoid negative numbers
               	if (mut_type==MICROSAT) os << (10000+*cursite++) << " ";
               	else os << *cursite++;
            }
      		os << "\n";
      	}
   	}
   }
   os << "\n";
};

extern ostream& operator<<(ostream& os, const TTree& T) {
    TNode& root=T.tree[2*T.sample_size()-2];
    cout<<"\n\nTREE STRUCTURE";
    root.TempRecNodeOut(0);
}

void TNode::TempRecNodeOut(int lvl) {
    cout<<"\n";
    for(int i=0;i<lvl;i++)
        cout << "  ";
    cout<< "Node " << node_number << ", Age:" << time <<", Desc:[" << num_desc1 << ","<<num_desc2<<"]";
    if(desc1) desc1->TempRecNodeOut(lvl+1);
    if(desc2) desc2->TempRecNodeOut(lvl+1);
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


