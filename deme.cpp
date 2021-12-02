////////////////////////////////////////////////////////////////////////////////
//
// A class implementing a deme within a population, exchanging migrants with
// other demes
//
//
////////////////////////////////////////////////////////////////////////////////


#include "deme.h"
#include "arrays.h"
#include "cond_var.h"

#ifdef _SHOW_MUT_DISTR_
	#include <fstream>
	ofstream f_site_hits("mut_hits_distr.sum");
	static int count_simul=0;
#endif

//Debugging internal function
void SearchLine(TNode *anc, TNode *x, int *hits, int deme, int line, int depth);

//------------------------------------------------------------------------------

//Initialization of static members of TDeme
TMigrationMatrix* TDeme::Migrations=NULL;

//==============================================================================
//------------------------------------------------------------------------------
TDeme TDeme::operator=(TDeme D) {
   _growth_rate=D._growth_rate;
   _samples=D._samples;
   _deme_size=D._deme_size;
   _lineages=D._lineages;
   _id=D._id;
   _min_deme_size=D._min_deme_size;
   _max_deme_size=D._max_deme_size;
   _proba_coal=D._proba_coal;
   Migrations=D.Migrations;
   _is_migration=D._is_migration;
   NodeList.de_allocate_list();
   int size=D.NodeList.size();
   if (size) {                
      NodeList.allocate_list(size);
      for (int i=0; i<size; ++i)
         NodeList[i]=D.NodeList[i];
   }
   size=D.OriginalNodes.size();
   if (size) {
      OriginalNodes.allocate_list(size);
      for (int i=0; i<size; ++i)
         OriginalNodes[i]=D.OriginalNodes[i];
   }   
   _polym_sites=D._polym_sites;
   return *this;
}

//------------------------------------------------------------------------------
int TDeme::copy_node_pointers(const TTree& T, int base, const long& max_lin) {
   int curr_lin=base;
   _lineages=0;
   _next_samp=0;
   //This node list is dynamic and will change according to migrations and coalescences
   // (but only allocate memory the first time)
   if(!NodeList.size()) NodeList.allocate_list(max_lin); 
   //Keeps a copy of the original nodes
   if(!OriginalNodes.size()) OriginalNodes.allocate_list(_samples.count());
   for(int i=0; i<_samples.GetItemsInContainer(); i++) {
      for(int j=0; j<_samples[i]._size; j++, curr_lin++) {
         T(curr_lin)->set(_id, _samples[i]._age, _samples[i]._stat_grp);
         //Christian 08: only include modern samples originally
        	if(_samples[i]._age==0) NodeList[_lineages++]=T(curr_lin); 
        	OriginalNodes[curr_lin-base]=T(curr_lin);
      }
      if(_samples[i]._age>_oldest_samp) _oldest_samp=_samples[i]._age;
      if(!_next_samp || (_next_samp>_samples[i]._age && _samples[i]._age)) _next_samp=_samples[i]._age;
   }
   if(!_next_samp) _next_samp=-1;
   for (int i=curr_lin-base; i<max_lin; i++) NodeList[i]=NULL; //not originalnodes!!!
   return 1;
}
//------------------------------------------------------------------------------
//Check and implements migrations
int TDeme::send_migrants(TDeme& sink, const int& i, const long& time) {
   bool trace=0; int n=0;
   if (!_lineages) return 0;
   for (int j=0; j<_lineages; j++) {
      //Christian: only extant lineages migrate, but if there is migration, allow what
      // used to be the last lineage a chance to migrate
      while (mtr()<=(*Migrations)(_id,i) && NodeList[j]->get_time()<=time && j<_lineages) {
      	migrate(sink,j);
         if(trace) printf("\nTime=%d Migrate %d->%d, Line %d (%d.%d)",time,_id,sink._id,j, NodeList[j]->deme, NodeList[j]->node_number);
         n++;
      }
   }
   return 1;
};
//------------------------------------------------------------------------------
//Implements the coalescence process within population
//Christian: if 'hudson' is true then the next coalescent time may be appropriately chosen
//from an expontential distribution. The return value is the waiting time
//to the next coalescence. Otherwise, the return value is the number of lineages
//that coalesced in the current 'time'th generation
int TDeme::coalesce_lineages(const long& time, TTree& GeneTree, int pos, bool hudson) {
   bool trace=0;
   int n_coals=0;
   //Check for old samples
   if(time==_next_samp) {
		if(trace) { cout<<"\nAdding old samples from pool of " << _samples.count() << endl; }
		for(int i=0; i<_samples.count(); i++) {
			if(trace) cout<<" "<<i;
			if(OriginalNodes[i]->time==time) { if(trace) cout<<"*";
				NodeList[_lineages++]=OriginalNodes[i]; }
		}
		_next_samp=_oldest_samp+1;
		for(int i=0; i<_samples.GetItemsInContainer(); i++)
			if(_samples[i]._age>time && _samples[i]._age<_next_samp)
				_next_samp=_samples[i]._age;
		if(_next_samp==_oldest_samp+1) _next_samp=-1;
		if(trace) printf("\nGen%d Deme%d n=%d k=%d, next t=%d",time, _id, int(_deme_size),_lineages,_next_samp);
	}
   if (_lineages==1 && hudson) {
		printf("ERROR! Tree has coalesced without reaching the root!");
		printf("\n(%d lineages left in deme%d at Gen%d)",_lineages,_id,time);
		system("PAUSE");
		exit(1);
	}
   //If there is only one EXTANT lineage left, do not even try to coalesce it...
   if (_lineages<2) return 0;
   //if(trace) printf("\nGen%d Deme%d n=%d k=%d, next t=%d",time, _id, int(_deme_size),_lineages,_next_samp);

	if (_deme_size<=0) {
     	cout << "Deme " << _id << " has zero population!\n"; return 0; }
   //prob = k(k-1)/(2N) BUGBUGBUG! Only asymptotically true
   double proba_coal;
   int skipgen=0;
   if(!hudson) proba_coal=_proba_coal*0.5*_lineages*(_lineages-1);
   else {
     proba_coal=1.0;
     skipgen=(int) (-2.0*_deme_size*log(mtr())/float(_lineages*(_lineages-1)));
   }
    
   if (proba_coal==0.0) return 0;
   //Coalesce two lineages with probability proba
   while (mtr()<proba_coal && _lineages>1) {
   	if(trace) cout<< " p="<<proba_coal;
      //Pick two EXTANT lineages at random for a coalescence
      int pick1=int(mtr()*_lineages),pick2=int(mtr()*_lineages);
		pick2=(pick2+(pick1==pick2))%_lineages; //if same, pick next; if last, wrap
		if(NodeList[pick1]->time>time || NodeList[pick2]->time>time) { 
			cout<<"ERROR: tried to coalesce a non-existant node!"; 
			system("PAUSE"); return 0; 
		}
      if(trace) cout<<"\n   Coalescing " << NodeList[pick1]->node_number << " and " 
					<< NodeList[pick2]->node_number << endl; 
      TNode *desc1=NodeList[pick1], *desc2=NodeList[pick2], *ancestor=GeneTree(pos++);

      //Update node information
      ancestor->desc1=desc1;
      ancestor->desc2=desc2;
      ancestor->time=time+skipgen;
      desc1->ancestor=ancestor;
      desc2->ancestor=ancestor;

	   //Replace one coalescing node with the ancestor
      //and swap the second with the last node after updating the number
      //of remaining lineages 
      NodeList[pick1]=ancestor;
      NodeList[pick2]=NodeList[--_lineages];
      n_coals++;
      if(!hudson) proba_coal=_proba_coal*0.5*_lineages*(_lineages-1);
      else proba_coal=0.0;
   	if(trace) for(int i=0; i<_lineages; i++) printf("%d~%d ",NodeList[i]->node_number,NodeList[i]->time);
   } if(trace) printf("!");
   if(hudson) {
      if(trace) printf("\n*Skipping %d generations to %d",skipgen,time+skipgen);
      return skipgen;
   }
   return n_coals;
};

//------------------------------------------------------------------------------
void TDeme::reset() {
   //Restoring initial number of lineages
	_lineages=_samples.count();

   int size=NodeList.size();
   //Copying original nodes into dynamic list
	for (int i=0; i<_samples.count(); ++i) 
      NodeList[i]=OriginalNodes[i];
   //Filling up the remaining spaces in node list with empty pointers
   for (int i=_samples.count(); i<size; ++i) NodeList[i]=NULL; 
   _min_deme_size=MIN_SIZE;
   _max_deme_size=MAX_SIZE;
   _proba_coal=1.0/int(_deme_size);
}
//------------------------------------------------------------------------------
//A simple procedure to check whether there are migrations from this population
int TDeme::check_migrations(){
   if (Migrations) {
      _is_migration=false;
   	int size=Migrations->size();
      for (int i=0; i<size; ++i)
      if (_id!=i) {
      	if ((*Migrations)(_id,i)>0.0) {
         	_is_migration=true;
            return 1;
         }
      }
      return 0;
   }
   return 0;
};
//------------------------------------------------------------------------------
void TDeme::print_haplotypes_for_Arlequin(ostream & os, const int& n, const Mut_Type data_type) {
   char DNA_letters[4]={'A','G','C','T'};
   if(!OriginalNodes[0]) return;
   int len=OriginalNodes[0]->seq_length;

   os << "#Number of mutation hits per site\n"
   	<< "#Sites ";
   for (int j=0; j<len; ++j) {
      	os
            << setw(4)
            << setiosflags(ios::right)
         	<< j << " ";
      }      
   os << "\n";

   #ifdef _SHOW_MUT_DISTR_
   	f_site_hits << "Simulation #" << ++count_simul << "\t";
   #endif

   os << "#Hits  ";
   for (int j=0; j<len; ++j) {
      	os 
            << setw(4)
            << setiosflags(ios::right)
            << OriginalNodes[0]->hits[j] << " ";

         #ifdef _SHOW_MUT_DISTR_
      	f_site_hits
            << setiosflags(ios::right)
            << OriginalNodes[0]->hits[j] << "\t";
         #endif
      }
   os << "\n";

   #ifdef _SHOW_MUT_DISTR_
   f_site_hits  << "\n";
   #endif

   for (int i=0; i<_samples.count(); ++i) {
   	//Write id and frequency
   	os << n << "_" << i << " 1 ";
      TNode* CurNode=OriginalNodes[i];
      if (CurNode->sequence) {
      	int* cursite=CurNode->sequence->begin();
      	for (int j=0; j<len; ++j) {
         	if (data_type==DNA) os << DNA_letters[*cursite++];
            else
            	if (data_type==MICROSAT) os << (10000+*cursite++) << " ";
            	else os << *cursite++;
         }
         os << "\n";
      }
   }
};
//------------------------------------------------------------------------------
void TDeme::print_haplotypes_for_PAUP(ostream & os, const Mut_Type data_type) {
    bool trace=0;
    char DNA_letters[4]={'A','G','C','T'};
    int len=OriginalNodes[0]->seq_length;
    int n=_samples.count();
    if(trace) {
      cout<<"\nPrinting "<< n << " haplotypes from deme "<<_id;
      for(int i=0;i<n;i++) cout<<"\n"<<i<<"\t"<<OriginalNodes[i];
    }
    for (int i=0; i<n; ++i) {
        TNode* CurNode=OriginalNodes[i];
        if(trace) cout<<*CurNode;
        if (CurNode->sequence) { //Write id
            int* cursite=CurNode->sequence->begin();
            os << (CurNode->deme+1) << "." << CurNode->node_number << "."
					<< CurNode->time << "." << CurNode->stat_group << "\n";
            for (int j=0; j<len; ++j) {
  	            if (data_type==DNA) os << DNA_letters[*cursite++];
               else if (data_type==MICROSAT) os << (10000+*cursite++) << " ";
               else os << *cursite++;
            }
            os << "\n";
        }
    }
};
//------------------------------------------------------------------------------
//Christian 6/29/04 Oneliner, no breaks
ostream& operator<<(ostream& os, const TSamp& S) {
   char buf[80];
	sprintf(buf,"Size:%3d Age:%6d Deme:%2d StatGrp:%2d",S._size,S._age,S._deme,S._stat_grp);
	os<<buf;
   return os;
}

//------------------------------------------------------------------------------
ostream& operator<<(ostream& os, const TDeme& D) {
//Christian: modified to include more than one sample per deme
   int n=D._samples.GetItemsInContainer();
   os << "\n--DEME "<< D._id << "--"
	   << "\nDeme size   : " << (int)D._deme_size
   	<< "\nSample Grps : " << n;
	for(int i=0;i<n;i++) os << "\n - " << D._samples[i];
   os << "\nGrowth rate : " << D._growth_rate
      << "\n";
   return os;
}

//------------------------------------------------------------------------------
ostream&
operator<<(ostream& os, const THistoricalEvent& HE) {
	os	  << "\nTime             : " << HE.time
   	  << "\nSource           : " << HE.source
        << "\nSink             : " << HE.sink
        << "\nMigrants         : " << HE.migrants
        << "\nNew size         : " << HE.new_deme_size
        << "\nNew growth rate  : " << HE.new_growth_rate
        << "\nNew migr. matrix : " << HE.MigMat
        << endl;
   return os;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//==============================================================================
//------------------------------------------------------------------------------
long TDemeCollection::count_lineages() {
   long count=0;
   for (int i=0; i<_num_demes; ++i)
   	count+=Demes[i]._samples.count();
   return count;
};
//------------------------------------------------------------------------------
 TDemeCollection& TDemeCollection::operator=(const TDemeCollection& DC) {
   _is_migration=DC._is_migration;
   _num_demes=DC._num_demes;
   _cur_mig_mat=DC._cur_mig_mat;
   if (Demes) delete[] Demes;
   if (DC.Demes) {
      Demes= new TDeme[DC._num_demes];
      for (int i=0; i<DC._num_demes; ++i)
         Demes[i]=DC.Demes[i];
   }
   Events=DC.Events;
   MigrMatArray=DC.MigrMatArray;
   GeneTree=DC.GeneTree;
   if (CoalTimes) {
   	delete CoalTimes;
      CoalTimes=NULL;
   }
   if (DC.CoalTimes) {
      	try {
         	CoalTimes= new TMigrationMatrix();
        }
        catch (...) {
         	cout 	<< "TDeme::operator=(const TDeme& D): unable to allocate memory"
                  << endl;
            if (CoalTimes) {
      			delete CoalTimes;
      			CoalTimes=NULL;
      		}
            return *this;
        }
   	   *CoalTimes=*DC.CoalTimes;
   }
   if (PairDiff) {
   	delete PairDiff;
      PairDiff=NULL;
   }
   if (DC.PairDiff) {
      	try {
         	PairDiff= new TMigrationMatrix();
         }
         catch (...) {
         	cout 	<< "TDeme::operator=(const TDeme& D): unable to allocate memory"
                  << endl;
            if (PairDiff) {
      			delete PairDiff;
      			PairDiff=NULL;
      		}
            return *this;
         }
   	*PairDiff=*DC.PairDiff;
   }
   if (MinCoalTimes) {
   	delete MinCoalTimes;
      MinCoalTimes=NULL;
   }
   if (DC.MinCoalTimes) {
      	try {
         	MinCoalTimes= new TMigrationMatrix();
         }
         catch (...) {
         	cout 	<< "TDeme::operator=(const TDeme& D): unable to allocate memory"
                  << endl;
   			if (MinCoalTimes) {
      			delete MinCoalTimes;
      			MinCoalTimes=NULL;
      		}
            return *this;
         }
   	*MinCoalTimes=*DC.MinCoalTimes;
   }
   _polym_sites=DC._polym_sites;
   _tot_mut=DC._tot_mut;
   return *this;
}

//Check for migrations and update coalescence probabilities within each population
int TDemeCollection::check_migrations(){
   _is_migration=false;
   //Examine all populations in turns to know which
   //one experience migrations. !! Do not stop the checking after first migration found
   for (int i=0; i<_num_demes; ++i)  {
   	if (Demes[i].check_migrations()) {
      	_is_migration=true;
      }
   }
   if (_is_migration) return 1;
   return 0;
};

int TDemeCollection::check_growth(){
   for (int i=0; i<_num_demes; ++i) 
   	if (Demes[i].growth()!=0.0)
         return 1;
   return 0;
};

int TDemeCollection::implement_event(const int& cur_event) {
   THistoricalEvent& CurEvent=(*Events)[cur_event];
   //Step 1: Send migrants from source to sink, fixed Christian 7/28/04
   // Check also that the sink population has a deme size >0
   bool trace=0;
   if(trace) {cout<<CurEvent; system("PAUSE"); }
   if (CurEvent.migrants>0.0 && CurEvent.source!=CurEvent.sink) {
		if(trace) printf("The %d lineages in %d have a %f percent chance to go to %d.",Demes[CurEvent.source]._lineages,CurEvent.source,CurEvent.migrants*100.0,CurEvent.sink);
      for (int i=0; i<Demes[CurEvent.source]._lineages; ++i) {
         while(mtr()<CurEvent.migrants && Demes[CurEvent.source].NodeList[i]->time<=CurEvent.time) {
            Demes[CurEvent.source].migrate(Demes[CurEvent.sink],i);
            if(i>=Demes[CurEvent.source]._lineages) break;
         }
      }
   }
   //Step 2: Resize sink deme size
   Demes[CurEvent.sink].linear_resize(CurEvent.new_deme_size);
   //Step 3: Readjust sink growth rate
   Demes[CurEvent.sink].set_growth_rate(CurEvent.new_growth_rate);
   //Step 4: Set new migration matrix; fixed Christian 7/28/04
   if (_cur_mig_mat!=CurEvent.MigMat) {
      for(int i=0; i<_num_demes; ++i)
         Demes[i].set_migration_matrix(&(*MigrMatArray)[CurEvent.MigMat]);
      _cur_mig_mat=CurEvent.MigMat;
   }
   return 1;
};
//------------------------------------------------------------------------------
int TDemeCollection::build_tree(ostream &logfile, int simnum) {
   bool trace=0;
   long oldest=0,num_lineages=count_lineages();
   long pos=2*GeneTree.sample_size()-num_lineages;
   int i,j,coal,num_active_demes=0, num_events=Events->GetItemsInContainer(), cur_event=0;
   int next_event[num_events];

   //initialize next_event list
   if(num_events) {
   	next_event[0]=0;
	   for(i=1; i<num_events; i++) {
			j=i-1;
			while(j>=0 && (*Events)[i].time<(*Events)[next_event[j]].time)
				next_event[j+1]=next_event[j--];
			next_event[j+1]=i;
		}
	}
	if(trace) for(i=0;i<num_events;i++) cout<<(*Events)[i];
   //Loro_21_9_99 :
   //By default take the first migration matrix of the array; zero migration matrix
   // if no migration matrix is provided
   TDeme aDeme;
   //Copies the static pointer to the migration matrix in each deme
   aDeme.set_migration_matrix(&(*MigrMatArray)[0]);
   _cur_mig_mat=0;

   //Check that there are migrations in the genealogical process
   bool implement_migration=check_migrations(), implement_growth=check_growth();
   for(i=0; i<_num_demes; i++) {
		if(oldest<Demes[i]._oldest_samp) oldest=Demes[i]._oldest_samp;
		if(Demes[i]._lineages) num_active_demes++;
	}   		
   if(trace) {
      printf("\n\n**BUILDING TREE**\nInitial Conditions:\n");
      for(i=0; i<_num_demes; i++) cout<<Demes[i];
   }
   
   for (long time=0, j, k; num_lineages>1; ++time) {
      if(time<0) {printf("\n\nERROR: TIME = %d",time); system("PAUSE"); return(0); }
      if(logfile) logfile<<"\n"<<simnum<<"\t"<<time;
      //if(time % 10000==0) cout<< "\n - Time="<<time <<" Lineages=" << num_lineages;
      //Step 1: Check for possible historical events
      if(num_events && cur_event<num_events) {
			while ((*Events)[next_event[cur_event]].time==time) {
	      	implement_event(next_event[cur_event++]);
	      	implement_migration=check_migrations();
	      	implement_growth=check_growth();
	         num_active_demes=0;
			   for(i=0; i<_num_demes; i++) if(Demes[i]._lineages) num_active_demes++;
			   if(cur_event==num_events) break;
				if(trace && cur_event && (*Events)[next_event[cur_event-1]].time==time) {
					printf("\nTime %d Event(s) %d\nNew state:",time,next_event[cur_event-1]);
	      		for(i=0; i<_num_demes; i++) printf("\n-Deme%d: n=%d k=%d r=%f next samp=%d",
						i,Demes[i].deme_size(),Demes[i]._lineages,Demes[i]._growth_rate,Demes[i]._next_samp);
	      		//system("PAUSE");
				}
	      } // end while event
		} //end if events

      //Step 2: Migration round
      if (implement_migration) {
      	for (j=0; j<_num_demes; ++j) {
         	if (Demes[j].is_migration()) {
      			for (k=0; k<_num_demes; ++k) {
         			if (j!=k) Demes[j].send_migrants(Demes[k],k,time);
                }//end for
         	} //end if
      	} //end for
      } //end if

		if(!implement_migration && cur_event>=num_events && num_active_demes>1) {
			printf("\nERROR at time %d: lineages in two unconnected demes cannot coalesce. Please adjust input file.",time);
			system("PAUSE"); return 0;
		}
   	//Step 3: Coalescence round; Christian: speed things up if the simulation
   	//has reached a simple point by picking the next coalescent time from an
   	//exponential distribution
      for (j=0; j<_num_demes; ++j) {
         if(cur_event<num_events || implement_migration || implement_growth || oldest>time) {
				Demes[j]._lineages;
         	coal=Demes[j].coalesce_lineages(time, GeneTree, pos, 0);
            if (coal) {
            	num_lineages-=coal;
					if(trace) printf("\nTime %d; %d COALS in deme %d (%d lines left)",time,coal,j,num_lineages);
               pos+=coal;
            }
         } else {
				if(Demes[j]._lineages) {
         		time+=Demes[j].coalesce_lineages(time, GeneTree, pos++, 1)-1;
            	num_lineages--;
					if(trace) printf("\nTime %d; 1 COAL in deme %d (%d lines left)",time,j,num_lineages);
				}
         }
         if(logfile) logfile << "\t" << Demes[j].deme_size() << "\t" << Demes[j]._lineages;
      } //end for Demes

   	//Step 4: Resize deme size after growth for next generation
      if (implement_growth) {
        //if(trace) printf("\nTime %d Resize: ",time);
        for (j=0; j<_num_demes; ++j) {
      	 if (Demes[j].growth()) {
            if (!Demes[j].exponential_resize()) {
	               cout << "TDemeCollection::build_tree(): size of deme "
               	     << j << "is zero !" << endl;
                   return 0;
            } //end if resize error
          } //end if Demes.growth()
          //if(trace) printf("%d ",Demes[j].deme_size());
        } //end for Demes
      } //end if implement_growth
      
      //Christian: Step 5: Verify that the original lineages are either still in the
      //deme collection, or descendents of nodes in the deme collection.
/* Take it on faith for now...
      TNode *x, *y;
      int n_fnd;
      for (j=0; j<_num_demes; ++j) {
         //if(trace) printf("\n\nOriginal Nodes from Deme %d\n",j);
         for(i=0;i<Demes[j].sample_size();i++) {
            x=Demes[j].OriginalNodes[i];
            //if(trace) printf("\n%d. %X",i,x);
            n_fnd=0;
            //search for it among the current lineages + their desc in all demes
            for(int k=0;k<_num_demes; ++k) {
               for(int l=0; l<Demes[k]._lineages; l++) {
                  y=Demes[k].NodeList[l];
                  SearchLine(y,x,&n_fnd, k, l, 0);
               }
            }
            if (n_fnd!=1) { printf("**ERROR** Node %d in deme %d has vanished!", i,j); system("PAUSE"); return 0;}
         }
      } // */
      if(cur_event<num_events || implement_migration || implement_growth || oldest>time) {
			if(time>=MAX_GENERATIONS) { //We'll simulate a maximum of MAX_GENERATIONS generations
				cout<<"\n **FAILED TO COALESCE IN " <<MAX_GENERATIONS<<" GENERATIONS.";
				//system("PAUSE");
				break;
			}
		}
    } //end for generations
    if (num_lineages>1) return 0;
	 return 1;
};

void SearchLine(TNode *anc, TNode *x, int *hits, int deme, int line, int depth) {
   if (x==anc) {
      (*hits)++;
      /*if(depth) printf(" Gen%d desc of Line %d in Deme %d", depth, line, deme);
      else printf(" Line %d in Deme %d",line, deme);*/
   }
   if (anc->desc1) SearchLine(anc->desc1,x,hits,deme,line,depth+1);
   if (anc->desc2) SearchLine(anc->desc2,x,hits,deme,line,depth+1);
}

int TDemeCollection::create_demes(TMigrMatArray *MA, TDemeSizeArray *DS, TGrowthArray *GR,
                                  TEventArray *HE, TSampArray *SS) {
   MigrMatArray=MA;
   if(_num_demes!=(*MA)[0].size()) { 
      printf("ERROR: Migration matrix does not match number of demes!"); return 0; }
   if(!Demes) Demes=new TDeme[_num_demes];
   for (int i=0; i<_num_demes; ++i) {
      Demes[i].set_migration_matrix(&(*MigrMatArray)[0]); //same migmatr for all demes
      Demes[i].set_id(i); //different ids
   }
   if (DS->GetItemsInContainer()!=_num_demes) {
      printf("ERROR: Number of deme sizes does not match number of demes!"); return 0; }
   for (int i=0; i<_num_demes; ++i)
   	Demes[i].set_deme_size((*DS)[i]);
   if (GR->GetItemsInContainer()!=_num_demes) {
      printf("ERROR: Number of growth rates does not match number of demes!"); return 0; }
   for (int i=0; i<_num_demes; ++i) 
   	Demes[i].set_growth_rate((*GR)[i]);
	Events=HE;

   //Flush and reload samples
   for(int i=0; i<_num_demes; ++i) Demes[i].DumpSamples();
   int size=SS->GetItemsInContainer();
   TSamp s;
   for (int i=0; i<size; ++i) {
      s=(*SS)[i];
      if(s._deme>=_num_demes) {
         printf("\nERROR: Sample %d has invalid deme number (%d)",i,s._deme); return 0; }
      Demes[s._deme].AddSample(s);
      if(s._stat_grp>_num_stat_grps) _num_stat_grps=s._stat_grp;
   }
   //Count total initial number of lineages (sum of samples, ancient and modern)
   long num_lineages=count_lineages();
   if (!num_lineages) { cout << "\nERROR: No samples"; return 0; }
   if(!GeneTree.sample_size()) GeneTree.allocate_tree(num_lineages); //don't allocate twice
   for (int i=0, current_lineage=0; i<_num_demes; ++i) {
      Demes[i].copy_node_pointers(GeneTree, current_lineage, num_lineages);
      current_lineage+=Demes[i]._samples.count();
   }
   return _num_demes;
}

//------------------------------------------------------------------------------
//A procedure to sprinkle mutations once the tree branch lengths are known
const long& TDemeCollection::sprinkle_mutations(const float& mut_rate, const int& len,
						const Mut_Type& mut_type, const float& gamma_par,
                  const float& mut_ratio, const float& prop_sites,
                  const float& trans_rate, const int& range_constr, symb *muteq) {
  	//if(1) { GeneTree.count_descendants(); cout<<GeneTree; }
   //Get the root node
	TNode& root=GeneTree[2*GeneTree.sample_size()-2];
	//Set tree length to zero
   root.reset_tree_length();
   //Initialize hits per site to zero
   int size=root.hits.GetItemsInContainer(), to_add=len-size;
   //Loro_29_8_98 //Resize sequence length if needed
	for (int i=0; i<to_add; ++i) root.hits.Add(0);
   for (int i=0; i<len; ++i) root.hits[i]=0;
   root.set_mut_rate(mut_rate);
   //root.count_desc();
   //root.TempRecNodeOut(0);  //print tree
   _tot_mut=root.add_mutations(0, len, mut_type, gamma_par, trans_rate, range_constr, muteq);
   return _tot_mut;
}
//------------------------------------------------------------------------------
extern ostream& operator<<(ostream& os,const TDemeCollection& DC) {
    int i=0, numdemes=DC.num_demes();
    os << "\nDEMES"
       << "\n======================";
	for (; i<numdemes; ++i) {
   	    os << DC.Demes[i] << "\n---------------------------------------\n";
   }
   int size=DC.MigrMatArray->GetItemsInContainer();
   for (int i=0; i<size; ++i) {
   	os << "\nMIGRATION MATRIX " << i
   		<< "\n=================\n"
      	<< (*DC.MigrMatArray)[i]
      	<< "\n";
   }
   os << "\nHISTORICAL EVENTS"
   	<< "\n=================\n";

   for (int i=0, numevents=DC.Events->GetItemsInContainer(); i<numevents; ++i) {
   	os << "\nEvent #" << (i+1)
      	<< "\n----------\n\n"
         << (*DC.Events)[i]
         << "\n---------------------------------------\n";
	}
   return os;
};

//------------------------------------------------------------------------------
void
TDemeCollection::write_samples_to_Arlequin_file(ostream& os, const Mut_Type data_type) {
   for (int i=0; i<_num_demes; ++i)
   if (Demes[i].sample_size())  //Loro_20_9_99
   {
   	write_Arlequin_sample_header(i+1, Demes[i].sample_size(), os);
      Demes[i].print_haplotypes_for_Arlequin(os, i, data_type);
      write_Arlequin_sample_footer(os);
   }

}
//------------------------------------------------------------------------------
void
TDemeCollection::write_group_section_to_Arlequin_file(ostream& os) {
	os
   	<< "\n[[Structure]]\n"
      << "\n\tStructureName=\"Simulated data\""
      << "\n\tNbGroups=1"
      << "\n\tGroup={";
   for (int i=0; i<_num_demes; ++i) {
   	if (Demes[i].sample_size())   //Loro_20_9_99
   	os << "\n\t   \"Sample " << (i+1) << "\"";
   }
   os << "\n\t}" << endl;
};
//------------------------------------------------------------------------------
void
TDemeCollection::write_samples_to_PAUP_file(ostream& os, const Mut_Type data_type) {
   for (int i=0; i<_num_demes; ++i)
		Demes[i].print_haplotypes_for_PAUP(os, data_type);
}
//------------------------------------------------------------------------------
//Christian: The stat file!
void TDemeCollection::write_Statfile(ostream& os, const Mut_Type data_type) {
	bool trace=0;
 	THapArray myhaps[_num_stat_grps+2];
	summary_stats gs;
	long ind,mrca=GeneTree[2*GeneTree.sample_size()-2].time;
	TNode *x; T_HAP curhapc, curhap;

	//Gather haptypes, make a combined group in myhaps[0], and also stratified into myhaps[statgroup+1]
	for(int curdeme=0; curdeme<_num_demes; curdeme++) {
		if(trace) printf("\n -Gathering haplotypes for deme %d\n   - Nodes:",curdeme);
		for(int curnode=0; curnode<Demes[curdeme].OriginalNodes.size(); curnode++) {
			if(trace) printf("\n%d: ",curnode);
			x=Demes[curdeme].OriginalNodes[curnode];
         curhapc.hap=x->sequence;
         //combined group
         ind=myhaps[0].Find(curhapc); //search all sequences: myhaps(0)
         if(ind==INT_MAX) {
				curhapc.n=1;
           	curhapc.is_private=1;
           	curhapc.first_deme=curdeme;
		   	myhaps[0].Add(curhapc);
		   	if(trace) printf(" New in deme %d",curdeme);
         } else {
				myhaps[0][ind].n++;
				if(myhaps[0][ind].first_deme!=curdeme) myhaps[0][ind].is_private=0;
				if(trace) printf(" = %d in deme %d",ind,curdeme);
			}
			//individual stat group
			curhap.hap=x->sequence;
			ind=myhaps[x->stat_group+1].Find(curhap); //search just this sequence's statgroup
			if(ind==INT_MAX) {
	 		 	curhap.n=1;
	 		  	curhap.is_private=1;
	         curhap.first_deme=curdeme;
			   myhaps[x->stat_group+1].Add(curhap);
			   if(trace) printf(". New in SG %d",x->stat_group);
			} else {
	         myhaps[x->stat_group+1][ind].n++;
				if(myhaps[x->stat_group+1][ind].first_deme!=curdeme) myhaps[x->stat_group+1][ind].is_private=0;
				if(trace) printf(". =%d in SG %d",x->stat_group);
			}
		}  //end for each node
	} //end for each deme
	if(trace) {
		for(int i=0;i<=_num_stat_grps;i++) {
			printf("\nGroup %d:",i);
			for(int j=0;j<myhaps[i+1].GetItemsInContainer();j++) printf(" %d",myhaps[i+1][j].n);
		}
		printf("\nCombined:");
		for(int j=0;j<myhaps[0].GetItemsInContainer();j++) printf(" %d",myhaps[0][j].n);
	}
	if(trace) printf("\n\nCalculating stats: Group");
   for(int curgroup=1; curgroup<=_num_stat_grps+1; curgroup++) {
		if(trace) printf(" %d",curgroup);
      CalcStats(&(myhaps[curgroup]), &gs, data_type);
      int ld=gs.mismatch.GetItemsInContainer();
      if(data_type==DNA)
          os << ","<<gs.num_haps << "," << gs.num_private_haps << "," << gs.seg_sites<<","
             << gs.pairws_dif<<","<< gs.hap_div <<","<<gs.nuc_div<<","<<gs.tajd <<","<<gs.fstar<<",{";
      else
         os << ","<<gs.num_haps << "," << gs.num_private_haps << "," << gs.hap_div<<","
         	<< gs.pairws_dif <<",{";
      for(int i=0; i<ld; i++) os<<gs.mismatch[i]<<" ";
      os << "},";

      for(int curgroupy=curgroup+1; curgroupy<=_num_stat_grps+1; curgroupy++) {
         CalcInterStats(&(myhaps[curgroup]),&(myhaps[curgroupy]),&gs,data_type);
         if(data_type==DNA)
            os << ","<<gs.num_haps << "," << gs.num_private_haps << "," << gs.pairws_dif<<","
               << gs.hap_div <<","<<gs.nuc_div<<","<<gs.tajd <<",";
      	else os << ","<<gs.num_haps << "," << gs.num_private_haps << "," << gs.pairws_dif <<","
				<< gs.tajd <<",";
      }
      myhaps[curgroup].Flush();
   }
	//combined stats
	if(trace) printf(" All");
	CalcStats(&(myhaps[0]),&gs, data_type);
	if(data_type==DNA) {
	   os << ","<<gs.num_haps << "," << gs.num_private_haps << "," << gs.seg_sites<<","
         << gs.pairws_dif<<","<< gs.hap_div <<","<<gs.nuc_div<<","<<gs.tajd <<","<<gs.fstar<<",{";
   } else os<<","<<gs.num_haps<<","<<gs.num_private_haps<<","<<gs.hap_div<<","<< gs.pairws_dif <<",{";
   for(int i=0;i<gs.mismatch.GetItemsInContainer();i++) os<<gs.mismatch[i]<<" ";
   os << "}," << mrca;
}
