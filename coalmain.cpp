/*//////////////////////////////////////////////////////////////////////////////
Copyright Laurent Excoffier, Genetics and Biometry Lab, University of Geneva
E-mail: laurent.excoffier@anthro.unige.ch
URL: http://anthropologie.unige.ch/~laurent

The use of the following source code in only allowed if it is kept intact
with its original banner
//////////////////////////////////////////////////////////////////////////////*/

#include "genealgy.h"
#include "public.h"
#include "time.h"
#include "migrmat.h"
#include "deme.h"
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include <unistd.h>
#include "coalmain.h"
#include "genstat.h"
#include "cond_var.h"
#include "algebra.h"

using namespace std;

class TInput {
   public:
	 int num_pop;
    int num_samp;
    int num_loci;
    Mut_Type data_type;
    double mut_rate;
    double transition_rate;
    double gamma_a;
    double mut_ratio;
    double prop_sites;
    double num_rate_categories;
    TNode first_node;
    int range_constraint;
    TDemeCollection myDemes;
    TEventArray myEvents;
    TDemeSizeArray *myDemeSizes;
 	 TGrowthArray *myGrowthRates;
    TSampArray *mySamples; //Christian 6/29/04
    TMigrMatArray *myMigrationMatrices;  //Loro_16_9_99
    vector<TPrior> myPriors;
    symb *muteq;

   TInput(void) {
    	myDemeSizes=new TDemeSizeArray(10,10);
    	myGrowthRates=new TGrowthArray(10,10);
   	mySamples=new TSampArray(10,10); //Christian 6/29/04
      myMigrationMatrices=new TMigrMatArray(1,10);  //Loro_16_9_99
      muteq=NULL;
   }
   ~TInput(void) {
      delete myDemeSizes;
      delete mySamples;
      delete myMigrationMatrices;
      delete myGrowthRates;
   }
   
   void *LookUp(TPrior *tp) {
      int id[3],x; float *y;
      switch(tp->param) {
         case deme_size: return &((*myDemeSizes)[tp->num]); break;
         case sample_age: return &((*mySamples)[tp->num]._age); break;
         case sample_stat_grp: return &((*mySamples)[tp->num]._stat_grp); break;
         case growth_rate: return &((*myGrowthRates)[tp->num]); break;
         case mig_rate: 
            x=tp->num;
            for(int i=2;i>=0;i--) {
               id[i]=x%num_pop;
               x=(x-id[i])/num_pop;
            }
            y=&((*myMigrationMatrices)[id[0]](id[1],id[2]));
            return y; break;
         case event_time: return &myEvents[tp->num].time; break;
         case event_src: return &myEvents[tp->num].source; break;
         case event_sink: return &myEvents[tp->num].sink; break;
         case event_migrants: return &myEvents[tp->num].migrants; break;
         case event_size: return &myEvents[tp->num].new_deme_size; break;
         case event_growth: return &myEvents[tp->num].new_growth_rate; break;
         case event_mig_matr: return &myEvents[tp->num].MigMat; break;
         case pmut_rate: return &mut_rate; break;
         case abstract: return &tp->curval; break;
         default: return NULL; break;
      }
   }
   
   void write_Statfile(ostream& os) {
      myDemes.write_Statfile(os,data_type);
      if(int np=myPriors.size()) {
         os<<","; //empty column under "PRIORS"
         for(int curprior=0; curprior<np; curprior++) {
            switch(myPriors[curprior].type) {
               case 1: os<<","<<int(myPriors[curprior].curval+.5); break;
               case 2: os<<","<<long(myPriors[curprior].curval+.5); break;
               case 3: os<<","<<float(myPriors[curprior].curval); break;
               case 4: os<<","<<myPriors[curprior].curval; break;
            }
         }
      }
      os<<"\n";
   }
   void write_StatfileHeader(ostream& os) {
      if(data_type==DNA) {
        	for(int i=0;i<=myDemes.stat_grps();i++) {
       		os << "GROUP "<<i<<",Haptypes,PrivHaps,SegSites,PairDiffs,HapDiver,NucltdDiv,TajimasD,F*,MismatDist,";
               for(int j=i+1; j<=myDemes.stat_grps();j++)
    		      os <<i<< " VS "<<j<<",PrivTo"<<i<<",PrivTo"<<j<<",PairDiffs,MeanDiv(Hs-bar),PoolDiv(Ht),Fst,";
           }
           os << "COMBINED,Haptypes,PrivHaps,SegSites,PairDiffs,HapDiver,NucltdDiv,TajimasD,F*,MismatDist,MRCA";
      } else {
         for(int i=0;i<=myDemes.stat_grps();i++) {
            os << "SAMPLE "<<i+1<<",Alleles,PrivAlleles,Heterozyg,AllelicVar,RepFreq,";
            for(int j=i+1; j<=myDemes.stat_grps();j++)
               os <<i<< " VS "<<j<<",PrivTo"<<i<<",PrivTo"<<j<<",G^2,Rst,";
         }
         os << "COMBINED,Alleles,PrivAlleles,Heterozyg,AllelicVar,RepFreq,MRCA";
      }
      if(int np=myPriors.size()) {
         char x[20]; 
         os<<",PRIORS";
         for(int i=0; i<np; i++) {
            param2str(x,myPriors[i].param);
            if(myPriors[i].param!=mig_rate)
               os << "," << x << " " << myPriors[i].num;
            else {
               int id[3], j;
               j=myPriors[i].num;
               for(int k=2;k>=0;k--) {
                  id[k]=j%num_pop;
                  j=(j-id[k])/num_pop;
               }
               os << ","<<x<<" "<<id[0]<<" ("<<id[1]<<"->"<<id[2]<<")";
            }
         }
      }
      os<<"\n";
   }
   
   void write_ArlFile(ostream& os) {
      int i;
      os << "#Arlequin input file written by SerialSimCoal.exe\n\n"
			<< "#Simulation parameters:\n"
         << "#======================\n"
         << "\n#Deme sizes\n";

      for (i=0; i<num_pop; ++i)
			os << "#Deme #" << i << "\t" << (*myDemeSizes)[i] << "\n";
      os << "\n#Sample sizes\n";
		for (i=0; i<num_samp; ++i) 
			os << "#Deme #" << (*mySamples)[i]._deme << "\t" << (*mySamples)[i]._size << "\n";
      os << "\n#Growth rates\n";
		for (i=0; i<num_pop; ++i)
			os << "#Deme #" << i << "\t" << (*myGrowthRates)[i] << "\n";
      for (i=0; i<myMigrationMatrices->GetItemsInContainer(); ++i) 
      	os << "\n#Migration matrix " << i << "\n"
				<< (*myMigrationMatrices)[i];
      os << "\n\n#Historical events\n";
		for (i=0; i<myEvents.GetItemsInContainer(); ++i)
			os << "#Event #" << i << "\n" << myEvents[i] << "\n";
		os << "\n#Data type							 : ";
      switch (data_type) {
      	case DNA: os << "DNA (Transition rate="<< setprecision(3)
                             << (transition_rate*100.0) << "%)\n"; break;
         case RFLP: os << "RFLP\n"; break;
         case MICROSAT:
            if (range_constraint)
               os << "MICROSAT (Range constraint : " << range_constraint << " ["
						<< (10000+first_node.min_mic) << ";" << (10000+first_node.max_mic) << "])\n";
            else os << "MICROSAT (No range constraint)\n"; 
            break;
      }
      os
      	<< "\n#Mutation rate/generation  : " << data_type
         << "\n#Number of loci to simulate: " << num_loci
       	<< "\n#Gamma parameter           : " << gamma_a << "\n\n";

      //Loro_20_9_99: Get the number of demes with non-empty samples
      int num_real_samples=0;
      for (int d=0; d<num_pop; ++d)
      	if (myDemes[d].sample_size()) ++num_real_samples;

      write_Arlequin_header(num_real_samples, data_type, os);
	   myDemes.write_samples_to_Arlequin_file(os, data_type);
      myDemes.write_group_section_to_Arlequin_file(os);
   }
    
    void RepickPriors(void) {
		double v[myPriors.size()];
        for(int i=0; i<myPriors.size(); i++) {//set all priors first
            if(!myPriors[i].formula) myPriors[i].pick(LookUp(&(myPriors[i])));
            v[i]=myPriors[i].curval;
		} 
        for(int i=0; i<myPriors.size(); i++) //then calculate formulae
            if(myPriors[i].formula) {
                myPriors[i].curval=myPriors[i].formula->solve(v);
                void* loc=LookUp(&(myPriors[i]));
                myPriors[i].set(loc);
                //printf("\n- Prior %d: %f",i+1,(float)(*myPriors)[i].curval); //KILL ME
            }
   }
};

//Internal functions
bool LoadFile(string path, TInput *t);
void GetFiles(int argc, char* argv[], vector<string> *infile_name, int *num_rand_samples, int *num_files, int *flags);
string remove_extension(string s) { return(s.erase(s.find('.'))); }

int main(int argc, char *argv[]) {
   bool trace=0;
   int flags;

   vector<string> infile_name;
   int num_rand_samples=1, num_files=0, i;

   printf("+----------------------------+\n|   BAYSIAN SERIAL SIMCOAL   |\n| by the Hadly Lab, Stanford |\n|       via SimCoal 1.0      |\n+----------------------------+");
   
   GetFiles(argc, argv, &infile_name, &num_rand_samples, &num_files, &flags); //Handle the input
   int _PRINT_ARLEQUIN_OUTPUT_=flags&0x1, _PRINT_SEQ_PAUP_=flags&0x2, _OUT_TRACE_=flags&0x4, _NEXUS_=flags&0x8;
   int _PRINT_TREE_PAUP_=flags&0x16;
   if(infile_name.size()==0) return 1;
   if(_NEXUS_) {
		printf("\nUsing NEXUS input, NO SIMULATIONS WILL BE RUN!");
		//read nexus file
		//make stat file
		exit(0);
	}
   printf("\nARLEQUIN output %s",_PRINT_ARLEQUIN_OUTPUT_ ? "full" : "suppressed");
   printf("\nPAUP output %s",_PRINT_SEQ_PAUP_ ? "full" : "suppressed");
   printf("\nDeme tracer is %s",_OUT_TRACE_ ? "enabled" : "disabled");
	if(_PRINT_TREE_PAUP_) printf("\nMaking tree files"); else printf("\nTree files are suppressed");

   //Reads and generate genealogies for the required number of files
   for (int nf=0; nf<num_files; ++nf) {
      printf("\n-----------\nNow reading %s:\n",infile_name[nf].c_str());
      //Creating input and output file names
		string outfile_name, arl_batch_file_name,
	             arl_generic_name, paup_batch_file_name,
	             paup_generic_name, Stat_file_name,
	             paup_name, log_file_name,
	             paup_true_trees_name, paup_mut_trees_name, temp_in;

		outfile_name		= infile_name[nf]+".gen";
		arl_batch_file_name = infile_name[nf]+".arb";
		arl_generic_name	= infile_name[nf];
		paup_generic_name	= infile_name[nf];
		paup_batch_file_name= infile_name[nf]+".bat";
		paup_name			= infile_name[nf]+".paup";
		paup_true_trees_name= infile_name[nf]+"_true_trees.trees";
		paup_mut_trees_name	= infile_name[nf]+"_mut_trees.trees";
		log_file_name 		= infile_name[nf]+".log";
		Stat_file_name		= infile_name[nf]+"_stat.csv";
		temp_in				= infile_name[nf]+".par";
		
		ofstream ArlBatch,PaupBatch,PaupFile,LogFile,PaupTrueTreesFile,PaupMutTreesFile;
		//Opening result files
		ofstream ofs(outfile_name.c_str());	
 		if (!ofs) { cout << "\nUnable to open output file (" << outfile_name << ")\n"; return 1; }
		if(_PRINT_ARLEQUIN_OUTPUT_) ArlBatch.open(arl_batch_file_name.c_str());
		if(_PRINT_SEQ_PAUP_) PaupBatch.open(paup_batch_file_name.c_str());
      if(_OUT_TRACE_) LogFile.open(log_file_name.c_str());

      TInput t;
      //Christian: Read all the input
      if(!LoadFile(temp_in,&t)) return 0;
      t.RepickPriors();
      cout<<"\n\nThis file contains "<<t.myPriors.size() << " priors";
      if(trace) for(int i=0; i<t.myPriors.size(); i++) {printf("\nPrior #%d\n",i+1); t.myPriors[i].Print(); }
      if(trace) system("PAUSE");
      cout<<"\nInitializing Demes";
      //Initialization of the demes
      if (!t.myDemes.create_demes(t.myMigrationMatrices, t.myDemeSizes, t.myGrowthRates, 
         &t.myEvents, t.mySamples)) return 0;
		cout<<t.myDemes;
      //Output initial conditions
   	ofs 	<< "\nDemes initial conditions"
         	<< "\n========================\n"
         	<< t.myDemes << endl;
      if(_PRINT_SEQ_PAUP_) PaupFile.open(paup_name.c_str());
      //Christian: stat file
      ofstream StatInfile(Stat_file_name.c_str());
      if(!StatInfile) { cout<<"Could not open statfile!"; return(1); }
      t.write_StatfileHeader(StatInfile);
      if(LogFile) {
			LogFile<<"Sim\tTime";
			for(int i=0; i<t.myDemes.num_demes(); i++) LogFile<<"\tDeme"<<i<<"\tk";
		}

   	cout << "\nBuilding " << num_rand_samples << " genealogies ...\n";
   	if(trace) system("PAUSE");

      if(_PRINT_TREE_PAUP_) {
         PaupTrueTreesFile.open(paup_true_trees_name.c_str());
         PaupMutTreesFile.open(paup_mut_trees_name.c_str());
         write_PAUP_trees_section_header(PaupTrueTreesFile);
         write_PAUP_trees_section_header(PaupMutTreesFile);
      }
      int tot_num_nodes=0;
      for (i=0; i<t.num_samp; ++i) tot_num_nodes+=(*t.mySamples)[i]._size ;
   
      if(_PRINT_SEQ_PAUP_) 
         write_PAUP_header(PaupFile,log_file_name.c_str(),tot_num_nodes,t.num_loci,
            t.num_pop,t.mut_rate,t.gamma_a,t.transition_rate);

      float mut_rate_per_locus=0.0;
      if (t.num_loci) mut_rate_per_locus=t.mut_rate/t.num_loci; //was data_type
      Mut_Model mut_mod=(fabs(t.gamma_a)>1e-7) ? K80_GAMMA : K80_NOGAMMA;
      
      //Build unequal mutation rates if needed
      if (!t.gamma_a==0.0 && t.data_type==DNA) 
      	t.first_node.mut_rates.get_new_rates(t.num_loci, t.gamma_a, t.mut_ratio, t.prop_sites, (int)t.num_rate_categories);
	
      my_float tree_mut_length=0.0, tree_mut_length2=0.0, tree_exp_length=0.0, tree_exp_length2=0.0;
      int tot_cases=0;
      while(++tot_cases<num_rand_samples) {
   	   if(trace) cout<<"\nRandom Sample " << tot_cases <<"\n - Building tree...";
   		if (t.myDemes.build_tree(LogFile,tot_cases)) { //HERE WE GO!!
         	//Add mutations
         	if(trace) cout<<"done.\n - Sprinkling Mutations..."; 
      		long mut=t.myDemes.sprinkle_mutations(t.mut_rate, t.num_loci, t.data_type,
                    t.gamma_a, t.mut_ratio, t.prop_sites, t.transition_rate, t.range_constraint, t.muteq);
      		if(trace) cout<<"done.";
            tree_mut_length+=mut;
            tree_mut_length2+=mut*mut;
            tree_exp_length+=t.first_node.tree_exp_length();
            tree_exp_length2+=t.first_node.tree_exp_length()*t.first_node.tree_exp_length();

         if(_PRINT_SEQ_PAUP_) {
   			//Writing data to PAUP file
   			if(trace) cout<<"\n - Writing to PAUP file...";
            write_PAUP_replicate_number((tot_cases+1), PaupFile);
            write_PAUP_data_header(tot_num_nodes, t.num_loci, t.data_type, PaupFile);
            t.myDemes.write_samples_to_PAUP_file(PaupFile, t.data_type);
            write_PAUP_end_matrix(PaupFile);
			   write_PAUP_tree_header(PaupFile, (tot_cases+1));
            t.myDemes.print_gene_tree(PaupFile, MUT_RATE, mut_rate_per_locus);
            write_PAUP_end(PaupFile);
            if (t.num_loci) {
           	   PaupFile 
                  << "\n\t[Tree length = " << (t.first_node.tree_exp_length()/t.num_loci)
                  << "\n\tNumber of sites hit by mutations = "  << t.first_node.count_polym_sites()
                  << "]\n";
            }
            write_PAUP_block(PaupFile,mut_mod);
            if(trace) cout<<"done.";
         }

         if(_PRINT_TREE_PAUP_) {
   			if(trace) cout<<"\n - Writing to PAUP Tree...";
   			write_PAUP_tree_name(PaupTrueTreesFile, (tot_cases+1));
            write_PAUP_tree_name(PaupMutTreesFile, (tot_cases+1));
            t.myDemes.print_gene_tree(PaupTrueTreesFile, GENERATIONS, mut_rate_per_locus);
            t.myDemes.print_gene_tree(PaupMutTreesFile, NUM_MUT, mut_rate_per_locus);
            if(trace) cout<<"done.";
         }
         if (!((tot_cases+1)%10))
      		cout << "\nGenealogy # " << (tot_cases+1) << "/" << num_rand_samples;

         
   		if(_PRINT_SEQ_PAUP_) {
		   if(trace) cout<<"\n - Writing to PAUP Sequence...";
		   //Build paup file name               
            char PaupOut[150], PaupLog[150];
            sprintf(PaupOut,"%s_%d.pau",paup_generic_name.c_str(),tot_cases);
            sprintf(PaupLog,"%s_%d.log",paup_generic_name.c_str(),tot_cases);
      		ofstream PaupOutFile(PaupOut);
            PaupBatch << "paup " << PaupOut << "\n";
            write_PAUP_header(PaupOutFile,PaupLog,tot_num_nodes,t.num_loci,
               t.num_pop,t.mut_rate,t.gamma_a,t.transition_rate);
            write_PAUP_replicate_number((tot_cases+1), PaupOutFile);
            write_PAUP_data_header(tot_num_nodes, t.num_loci, t.data_type,PaupOutFile);
            t.myDemes.write_samples_to_PAUP_file(PaupOutFile, t.data_type);
				write_PAUP_end_matrix(PaupOutFile);
				write_PAUP_tree_header(PaupOutFile, (tot_cases+1));
            t.myDemes.print_gene_tree(PaupOutFile, MUT_RATE, mut_rate_per_locus);
            write_PAUP_end(PaupOutFile);
            	PaupOutFile
            		<< "\n\t[Tree length = " << (t.first_node.tree_exp_length()/t.num_loci)
                  << "\n\tNumber of sites hit by mutations = "  << t.first_node.count_polym_sites()
                  << "]\n";
				write_PAUP_block(PaupOutFile,mut_mod);
				write_PAUP_footer(PaupOutFile);  
            if(trace) cout<<"done.";
         }

         if(_PRINT_ARLEQUIN_OUTPUT_) {
			if(trace) cout<<"Writing ARLEQUIN file...";
            //Build Arlequin file name
            char ArlOut[150];
            sprintf(ArlOut,"%s_%d.arp",paup_generic_name.c_str(),tot_cases);
            ofstream ArlInFile(ArlOut);
            ArlBatch << ArlOut << "\n";
            t.write_ArlFile(ArlInFile);
            ArlInFile.close();
			if(trace) cout<<"done.";
       	}

			if(trace) cout<<"\n - Writing to Statfile...";
      	t.write_Statfile(StatInfile);
   	   if(trace) cout<<"done.";
        } //end if !tree

      	if(trace) cout<<"\n - Repicking from Priors...";
      	t.RepickPriors();
        if(trace) cout<<"done.";
        t.myDemes.create_demes(t.myMigrationMatrices, t.myDemeSizes, t.myGrowthRates, 
            &t.myEvents, t.mySamples);
   	    if(trace) cout<<t.myDemes;
   	    if(trace) system("PAUSE"); 
      } //end random coal sample
      if (t.num_loci) {
         tree_exp_length/=t.num_loci;
         tree_exp_length2/=t.num_loci*t.num_loci;
      } else {
         tree_exp_length=0;
         tree_exp_length2=0;
      }

    if(_PRINT_SEQ_PAUP_) write_PAUP_footer(PaupFile);   
    if(_PRINT_TREE_PAUP_) {
      write_PAUP_end(PaupTrueTreesFile);
      write_PAUP_end(PaupMutTreesFile);
      PaupTrueTreesFile
        << "\nMean length of the trees : "
        << (tree_exp_length/num_rand_samples);
      if (num_rand_samples>1)
      	PaupTrueTreesFile
            << " +- "
         	<< setiosflags(ios::showpoint | ios::fixed | ios::right)
         	<< setw(10)
         	<< setprecision(5)
            << ( sqrt((tree_exp_length2 - tree_exp_length*tree_exp_length/num_rand_samples)/
            			 (num_rand_samples-1)) )
            << endl;

      PaupMutTreesFile
        << "\nMean length of the trees : "
        << (tree_mut_length/num_rand_samples);
      if (num_rand_samples>1)
      	PaupMutTreesFile
            << " +- "
         	<< setiosflags(ios::showpoint | ios::fixed | ios::right)
         	<< setw(10)
         	<< setprecision(5)
            << ( sqrt((tree_mut_length2 - tree_mut_length*tree_mut_length/num_rand_samples)/
            			 (num_rand_samples-1)) )
            << endl;
    }
        
   	cout << "\n\ndone in " <<  1000.0*float(clock())/float(CLOCKS_PER_SEC) << "ms\n" << endl;
	
    //Print to output file
    ofs   << "\nFinal results\n=============\n";

    ofs 	<< "\nData type : ";
    switch (t.data_type) {
      	case DNA 		: ofs << "DNA"; break;
      	case RFLP 		: ofs << "RFLP"; break;
      	case MICROSAT 	: ofs << "MICROSAT"; break;
    }

      if (t.data_type==DNA) {
      	if (t.gamma_a>0.0) {
         	ofs	<< "\nHeterogeneity of mutation rates: Gamma value = "
            	<< setiosflags(ios::showpoint | ios::fixed | ios::right)
         		<< setw(10)
         		<< setprecision(5)
         		<< t.gamma_a;
         }
         else ofs << "\nNo heterogeneity of mutation rates";
         ofs    << "\nProportion of transitions : "
            	<< setiosflags(ios::showpoint | ios::fixed | ios::right)
         		<< setw(10)
         		<< setprecision(5) << (100*t.transition_rate) << " (ts/tv=";
         if (t.transition_rate!=1) ofs << (t.transition_rate/(1-t.transition_rate)) << ")";
      	 else ofs << "inf";
      }
            
      ofs   << "\nNumber of loci : " << t.num_loci;
      ofs	<< "\nTotal mutation rate per generation : " << t.data_type;
      ofs	<< "\nMutation rate per locus per generation : " << t.data_type/t.num_loci;
	
      ofs
      	<< "\nMean length of the trees (subst. rate per site) : " 
         	<< setiosflags(ios::showpoint | ios::fixed | ios::right)
         << (tree_exp_length/num_rand_samples);
      if (num_rand_samples>1)
      	ofs << " +- "
         	<< setiosflags(ios::showpoint | ios::fixed | ios::right)
         	<< setw(10)
         	<< setprecision(5)
            << ( sqrt((tree_exp_length2 - tree_exp_length*tree_exp_length/num_rand_samples)/
            			 (num_rand_samples-1)) )
            << endl;

      ofs
      	<< "\nMean length of the trees (total number of mutations) : "
         << (tree_mut_length/num_rand_samples);
      if (num_rand_samples>1)
      	ofs << " +- "
         	<< setiosflags(ios::showpoint | ios::fixed | ios::right)
         	<< setw(10)
         	<< setprecision(5)
            << ( sqrt((tree_mut_length2 - tree_mut_length*tree_mut_length/num_rand_samples)/
            			 (num_rand_samples-1)) )
            << endl;
      ofs   << "\n\nTotal time for " << tot_cases
         	<< " simulations: " << 1000*float(clock())/float(CLOCKS_PER_SEC) << " ms"
         	<< endl;

      //Do not forget to reset the static member of TNode
      TNode resetNode;
      resetNode.reset_node_count();
      if(trace) system("PAUSE");
	}
   return 1;  
}

//------------------------------------------------------------------------------
int TMeanMat::compute_mean(const int& n) {
	_num_updates=n;
	if (!_num_updates) return 0;
	for (int i=0; i<_size; ++i) {
   	for (int j=0; j<_size; ++j) {
      	_data[i][j]/=_num_updates;
      }
   }
   return 1;
};

//Compute the s.d. from the values present in the matrix, assuming that
//those are the sum of square values
int
TMeanMat::compute_sd(const int& n, const TMeanMat& mean) {
   if (n<2) return 0;
	_num_updates=n   ;
	if (!_num_updates) return 0;
	for (int i=0; i<_size; ++i) {
   	for (int j=0; j<_size; ++j) {
      	double s2=_data[i][j],
         		 m=mean._data[i][j]; //assumes that mean contains mean values
      	_data[i][j]=sqrt( (s2-n*m*m)/(n-1) );
      }
   }
   return 1;
};

int
TMeanMat::update_with(const TMigrationMatrix& MM) {
	if (!MM.size() || _size!=MM.size()) return 0;
   for (int i=0; i<_size; ++i) {
   	for (int j=0; j<_size; ++j) {
      	_data[i][j]+=MM(i,j);
      }
   }return 1;
};

int TMeanMat::update_with_square(const TMigrationMatrix& MM) {
	if (!MM.size() || _size!=MM.size()) return 0;
   for (int i=0; i<_size; ++i) {
   	for (int j=0; j<_size; ++j) {
      	_data[i][j]+=MM(i,j)*MM(i,j);
      }
   }return 1;
};

bool LoadFile(char *path, TInput *t) {
   //Initialize input file
   ifstream ifs;
   ifs.open(path);
  	if (!ifs) {
      cout << "Unable to open input file (2)(" << path << ")" << endl;
      return 0;
   }
   bool is_anc_file=0;

   my_string line(20);
   line.read_to_delim(ifs);   //Reads a blank line
   cout<<endl<<line<<endl;

   TPrior *tp=BayesRead(&(t->num_pop),&ifs,none,0,1);        //Reads the number of populations to simulate
   if(tp) { cout<<"ERROR: Cannot use a variable number of populations!"; return 0; }
   t->myDemes.set_num_demes(t->num_pop);
   cout<<"# populations:"<<t->num_pop<<endl;

/*Christian 6/29/04: If the rest of the first line contains 'with ancient', then
  not all the samples are the same age, and appropriate action will be taken below.
  If not, the program is still backwards compatible */
   line.read_to_delim(ifs);
   if(line.contains("with ancient"))
   {   is_anc_file=1;
       cout<<"Contains ancient sampling information\n";
   } else cout<<"No ancient samples\n";

  	cout << "\nDeme sizes\n";
  	line.read_to_delim(ifs);   //Reads a blank line
	long cur_size; int i; 
  	for (i=0; i<t->num_pop; ++i) {
     	tp=BayesRead(&cur_size,&ifs,deme_size,i,2);
      if(tp) t->myPriors.push_back(*tp);
     	line.read_to_delim(ifs);
  		t->myDemeSizes->Add(cur_size);
  		cout << "Deme " << i << "\t" << (*t->myDemeSizes)[i] << "\n";
  	}

   //Reads the sample sizes
  	line.read_line(ifs);  //Reads a blank line
/*Christian 6/29/04: If the file contains ancient dates, then the format here
  is 'samp age deme stat_grp'. The number of sample groups is on the first line. Otherwise
  one sample group + stat group per deme is assumed */
   cout << "\nSamples\n----------";
   if(is_anc_file) { ifs>>t->num_samp; line.read_to_delim(ifs);
   } else t->num_samp=t->num_pop; //no priors
   TSamp cursamp;
  	for (i=0; i<t->num_samp; ++i) {
      if(is_anc_file) {
         tp=BayesRead(&cursamp._size,&ifs,none,0,2); if(tp) {printf("\nERROR: Variable sample sizes not allowed"); return 0;}
         tp=BayesRead(&cursamp._age,&ifs,sample_age,i,2); if(tp) t->myPriors.push_back(*tp);
         tp=BayesRead(&cursamp._deme,&ifs,none,0,1); if(tp) {printf("\nERROR: Variable home deme for samples not allowed"); return 0;}
         tp=BayesRead(&cursamp._stat_grp,&ifs,sample_stat_grp,i,1); if(tp) t->myPriors.push_back(*tp);
      } else {
         tp=BayesRead(&cur_size,&ifs,none,0,2); if(tp) {printf("\nERROR: Variable sample sizes not allowed"); return 0;}
         cursamp.SetVals(cur_size,0,i,i);
      }
      line.read_to_delim(ifs);
 		t->mySamples->Add(cursamp);
		cout << "\nSample " << i << " - " << (*t->mySamples)[i];
	}
//Christian 6/29/04: End modifications

   line.read_line(ifs);  //Reads a blank line
   cout << "\nGrowth rates\n";
   float cur_growth;
   for (i=0; i<t->num_pop; ++i) {
      tp=BayesRead(&cur_growth,&ifs,growth_rate,i,3); if(tp) t->myPriors.push_back(*tp);
      line.read_to_delim(ifs);
      t->myGrowthRates->Add(cur_growth);
      cout << "Deme" << i << ":\t" << (*t->myGrowthRates)[i] << "\n";
   }

   //Reading how many migration matrices to read
   int numMigMat=0; //Default
   line.read_line(ifs);  //Reads a blank line
   ifs >> numMigMat; //no priors
   cout<< numMigMat << " migration matrices";
   line.read_line(ifs);

   //Loro_16_9_99
   if (numMigMat) {
      for (int m=0; m<numMigMat; ++m) {
         //Reading migration matrix
         TMigrationMatrix myMigrationRates(t->num_pop);
   		line.read_line(ifs);  //Reads matrix comment line
   		for (i=0; i<t->num_pop; ++i) {
   			for (int j=0; j<t->num_pop; ++j) {
      			tp=BayesRead(&(myMigrationRates(i,j)),&ifs,mig_rate,m*t->num_pop*t->num_pop+i*t->num_pop+j,3);
               cout<<myMigrationRates(i,j)<<" ";
				}
         	line.read_to_delim(ifs); //Reads until the end of line
         	cout<<"\nPriors="<<t->myPriors.size()<<"\n";
         }
         t->myMigrationMatrices->Add(myMigrationRates);
   		cout << "\nMigration matrix\n";
   		cout << (*t->myMigrationMatrices)[m];
      	//line.read_to_delim(ifs); //Reads until the end of line
   	}
   } else { //Then a model without migration is assumed
   	TMigrationMatrix myMigrationRates(t->num_pop); //filled with zeroes
      t->myMigrationMatrices->Add(myMigrationRates);
   }

	int num_events;
   cout << "\nHistorical events (";
   line.read_line(ifs); //Reads comment line
	ifs >> num_events;   //Reads the number of historical events to read, no priors
	line.read_line(ifs); //Reads the string " historical events\n"
	cout<<num_events<<")\n";
	for (i=0; i<num_events; ++i) {
		THistoricalEvent curevent;
		tp=BayesRead(&curevent.time,&ifs, event_time,i,2); if(tp) t->myPriors.push_back(*tp);
   	tp=BayesRead(&curevent.source,&ifs, event_src,i,1); if(tp) t->myPriors.push_back(*tp);
      tp=BayesRead(&curevent.sink,&ifs, event_sink,i,1); if(tp) t->myPriors.push_back(*tp);
      tp=BayesRead(&curevent.migrants,&ifs, event_migrants,i,3); if(tp) t->myPriors.push_back(*tp);
      tp=BayesRead(&curevent.new_deme_size,&ifs, event_size,i,3); if(tp) t->myPriors.push_back(*tp);
      tp=BayesRead(&curevent.new_growth_rate,&ifs, event_growth,i,3); if(tp) t->myPriors.push_back(*tp);
      tp=BayesRead(&curevent.MigMat,&ifs, event_mig_matr,i,1); if(tp) t->myPriors.push_back(*tp);

		cout<<curevent;
		t->myEvents.Add(curevent,false);  //BUG: Sorts by date, which means priors point to the wrong event if events not chronological in parfile!
		//for(int j=0;j<=i; j++) cout << "Event " << j << "\n" << t->myEvents[j] << endl; system("PAUSE"); //KILL ME!

		line.read_to_delim(ifs); //Reads until the end of line
		//cout<< curevent;
	}
	//for(i=0;i<num_events; i++) cout << "Event " << i << "\n" << t->myEvents[i] << endl; system("PAUSE"); //KILL ME!

	//Read mutation rate
   line.read_line(ifs); //Read comments
   double currate;
   tp=BayesRead(&currate, &ifs, pmut_rate,0,4);
	if(tp) { t->myPriors.push_back(*tp); t->muteq=tp->formula; tp->Print(); }
	t->mut_rate=currate;
   printf("\nMutation rate: %f / generation; %s",currate,(t->muteq)?"variable":"constant");
   line.read_to_delim(ifs); //Reads until the end of line

   //Read number of loci
   line.read_line(ifs);
   ifs >> t->num_loci; //no priors
   cout << "\nNumber of loci to simulate : " << t->num_loci;
   line.read_to_delim(ifs); //Reads until the end of line

   //Reads data type
   t->transition_rate=0.0;
   //Loro_26_2_99
   t->range_constraint=0; //A value of zero assumes no range constraint
   line.read_line(ifs);  //Read comment line
   bool par=line.read_to_delim(ifs, ' ');
   line.to_upper();
   //Loro_26_2_99!
   if (line.contains("MICROSAT")) {
   	t->data_type=MICROSAT;
   	if(par) ifs >> t->range_constraint; //Extract range constraint for microsat data
   	line.read_to_delim(ifs);
   } else if (line.contains("DNA")) {
   	t->data_type=DNA;
      if(par) ifs >> t->transition_rate; //Extract transition rate for DNA
   	line.read_to_delim(ifs);
   } else if (line.contains("RFLP")) t->data_type=RFLP;
   // Loro_28_2_99 : Be careful: only one node allowed to be created before building the tree
   t->first_node.reset_node_count(); //Do not forget to reset the node count...

   if (t->data_type==MICROSAT && t->range_constraint) {//Fix minimum and maximum size for microsats
      double y=2.0, x=t->range_constraint;
      if (fmod(x,y)==0.0) { //Even number
      	t->first_node.min_mic=1-(t->range_constraint/2);
         t->first_node.max_mic=t->range_constraint/2;
      } else {
      	t->first_node.min_mic=(1-t->range_constraint)/2;
         t->first_node.max_mic=(t->range_constraint+1)/2; //bug?!
      }
   }

   cout << "\nData type : ";
   switch (t->data_type) {
    	 case DNA: cout << "DNA\n" << "Transition rate : " << t->transition_rate << "\n";
 			 	break;
      case RFLP: cout << "RFLP\n"; break;
      case MICROSAT:
         	if (t->range_constraint)
         		cout << "MICROSAT"
          			 << "\nRange constraint : " << t->range_constraint
                   << " (" << (10000+t->first_node.min_mic) << ";"
                   << (10000+t->first_node.max_mic) << ")\n";
         	else
         		cout << "MICROSAT" << "\n No range constraint\n";
            break;
       default : break;
   }

    //Loro_26_2_99
   t->gamma_a=0.0; t->mut_ratio=1.0; t->prop_sites=1.0; t->num_rate_categories=0;
   if (t->data_type==DNA) {
   	line.read_line(ifs); //Reads gamma parameter comment line
   	ifs >> t->gamma_a;
   	cout << "\nGamma parameter : " << t->gamma_a;
   	//If gamma parameter is negative, then it means that we have a
   	//two mutation rate model. So read the mutation ratio and the proportion
   	//of sites with the high mutation rate
   	if (t->gamma_a<0.0) ifs >> t->mut_ratio >> t->prop_sites;
   	else if (t->gamma_a>0.0) { ifs >> t->num_rate_categories;
   		cout << " " <<t->num_rate_categories << " (" << t->prop_sites*100 << "% at x"<<t->mut_ratio << ")";
		}
   }
   
	//Abstract parameters
	i=0;
   while(++i) {
      line.read_to_delim(ifs);
      tp=BayesRead(NULL,&ifs,abstract,i-1,3);
      if(!tp) break;
      t->myPriors.push_back(*tp);
   }
   cout<<"\nAbstract parameters: " << i-1;
   return 1;
}

void GetFiles(int argc, char* argv[], vector<string> *infile_name, int *num_rand_samples, int *num_files, int *flags) {
   string usage,batchFile; char *x;
   usage ="\nBayeSSC.exe usage:\n";
   usage+="\n   BayeSSC.exe -a -p -t -b --trees [batch file] -f [file] [# sims]";
   usage+="\n   -a: Output files for ARLEQUIN";
   usage+="\n   -b [batch file]: A file containing a list of .par files to process";
   usage+="\n   -f [file] [# sims]: A .par input file and the number of simulations to run";
   usage+="\n   -n [file] calculate statistics for a .nexus file";
   usage+="\n   -p: Output files for PAUP";
   usage+="\n   -t: Trace the size of demes and # of lineages in a .log file";
   usage+="\n   --trees: Supress tree files (will be created by default)";
   ifstream ifs_batch;

	int i=1;
	*num_files=0;
	while(i<argc) {
		if(strcmp(argv[i],"-a")==0) *flags+=1;
		if(strcmp(argv[i],"-p")==0) *flags+=2;
		if(strcmp(argv[i],"-t")==0) *flags+=4;
		if(strcmp(argv[i],"--trees")==0) *flags+=16;
		if(strcmp(argv[i],"-n")==0) { //nexus
            *flags+=8;
            batchFile=argv[++i];
            ifs_batch.open(batchFile.c_str());
            if (!ifs_batch) {
        		cout << "Unable to open nexus file (" << batchFile << ")\n";
        		cout << usage; return;
			}
			infile_name->push_back(batchFile);
		}
		if(strcmp(argv[i],"-b")==0) { //a file containing the name of the files to be processed
      	batchFile=argv[++i];
      	ifs_batch.open(batchFile.c_str());
      	if (!ifs_batch) {
        		cout << "Unable to open input file (" << batchFile << ")\n";
        		cout << usage; return;
      	} else {
     			string line; //Number of random samples to draw must be on first line
      		ifs_batch >> *num_rand_samples;
        		getline(ifs_batch,line);
      		while (ifs_batch.good()) {
     	    		getline(ifs_batch,line);
            	if (line!="") {
          	   	++(*num_files);
                	infile_name->push_back(remove_extension(line));
            	}
      		}
      		ifs_batch.close();
      	}
		}
		if(strcmp(argv[i],"-f")==0) { //a file to process and the number of simulations
      	*num_files=1;
      	batchFile=argv[++i];
   		ifs_batch.open((remove_extension(batchFile)+".par").c_str());
      	if (!ifs_batch) {
   	  		cout << "\nUnable to open input file: " << getcwd(NULL,0) << "\\" << remove_extension(batchFile)+".par" << endl;
        		cout << usage; return;
      	}
	      *num_rand_samples=atoi(argv[++i]);
   	   ifs_batch.close();
      	infile_name->push_back(remove_extension(batchFile));
         if (!(*num_rand_samples)) cout << usage;
		}
		x=argv[i]+1;
		if(argv[i][0]=='-' && strspn(x,"aptbf-")==0) cout << "\nERROR: Unrecognized flag: " << argv[i] << "\n" << usage;
		i++;
	}
   if(*num_files==0) {
      *num_files=1;
      string curName,tf;
      cout << "\n\nInput file name (*.par): ";
      cin >> curName;
      tf=remove_extension(curName)+".par";
   	ifs_batch.open(tf.c_str());
      while(!ifs_batch.is_open()) {
			cout<<"\n No file named " << getcwd(NULL,0) << "\\" << remove_extension(curName) << ".par\n\nInput file name: ";
      	cin >> curName;
   		ifs_batch.open(tf.c_str());
		}
      infile_name->push_back(remove_extension(curName));
      ifs_batch.close();
		
      cout << "\n\nNo. of random samples to generate : ";
      cin >> *num_rand_samples;
   }
}
