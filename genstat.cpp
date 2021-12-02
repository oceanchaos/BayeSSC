//Christian: Genetic statistics calculator; Simcoal compatable version
#include "genstat.h"
#include <math.h>

//Algebra functions
int tokenize(char *s, char *sep, char **toked, bool include_tok);
int numtokens(char *s, char *sep, bool echo);
symb *CreateTree(char **toked, int vars);

int SeqDif(int *seq1, int *seq2, int n) {
	if(!seq1 || !seq2) return 0;
	int ndif=0;
	for(int i=0; i<n; i++)
		if( *(seq1+i) != *(seq2+i)) ndif++;
	return ndif;
}

void CalcStats(THapArray *haps, summary_stats *st, const Mut_Type data_type) {
	float a1=0, a2=0, c2, e1S, e2S, vf,uf, n=0;
	float num_bases=haps[0][0].hap->end()-haps[0][0].hap->begin();
   int ss,nu_sing=0,nu=0,bp_freq[4],minmis=*(haps[0][0].hap->begin());
	int maxmis=(data_type==DNA)?0:minmis;
   st->num_haps=haps->GetItemsInContainer();
   st->num_private_haps=0;
   st->seg_sites=0;
   st->mismatch.Flush();
   if(st->num_haps<2 && data_type==DNA) {
		st->hap_div=st->nuc_div=st->tajd=0;
		st->mismatch.Add(haps[0][0].n*(haps[0][0].n-1)/2);
		return;
	}
	for(int i=0;i<st->num_haps;i++) {
        n+=float(haps[0][i].n);
        if(haps[0][i].is_private) st->num_private_haps++;
   }
   st->pairws_dif=0.0;
   //Calculate diversity, ...
   st->hap_div=1.0;
   st->mismatch.Add(0);
	for(int i=0; i<st->num_haps; i++) {
		st->hap_div-=pow((float)(haps[0][i].n/n),float(2.0));
		//...pairwise dif, mismatch...
		if(data_type==DNA) {
			for(int j=i+1; j<st->num_haps; j++) {
				ss=SeqDif(haps[0][i].hap->begin(),haps[0][j].hap->begin(),int(num_bases));
         	st->pairws_dif+=ss*haps[0][i].n*haps[0][j].n;
         	if(ss>maxmis) {
					for(int k=maxmis; k<ss; k++) st->mismatch.Add(0);
					maxmis=ss;
				}
         	st->mismatch[ss]+=haps[0][i].n*haps[0][j].n;
			}
       	st->mismatch[0]+=haps[0][i].n*(haps[0][i].n-1)/2;
      } else { 
         a1+=haps[0][i].n*pow(float(*(haps[0][i].hap->begin())),2);
         a2+=haps[0][i].n*float(*(haps[0][i].hap->begin()));
         //keep track of the repeat# range
			if(minmis>*(haps[0][i].hap->begin())) minmis=*(haps[0][i].hap->begin());
			if(maxmis<*(haps[0][i].hap->begin())) maxmis=*(haps[0][i].hap->begin());
      }
	}
	if(data_type==DNA) {
      if(n>1.0) {
			//..segsites and singletons...
			for(int i=0; i<int(num_bases); i++) {
            for(int j=0; j<4; j++) bp_freq[j]=0; //zero out bps
            bp_freq[*(haps[0][0].hap->begin()+i)]+=haps[0][0].n;
				for(int j=1; j<st->num_haps; j++) bp_freq[*(haps[0][j].hap->begin()+i)]+=haps[0][j].n;
				ss=-1;
				for(int j=0; j<4; j++) {
					nu_sing+=(bp_freq[j]==1); //count singletons
					ss+=(bp_freq[j]>0); //and total mutations (simple-minded method)
				}
				nu+=ss-1;
				if(ss>0) st->seg_sites++;
			}
         st->pairws_dif*=2/(n*(n-1));
         st->nuc_div=st->pairws_dif/num_bases;
         //...and Tajima's D.
         a1=a2=0.0;
       	for(double i=1.0; i<n; i++) { a1+=1.0/i; a2+=pow(i,-2.0); } //a1=an, a2=bn in F
         e1S=st->seg_sites*((n+1)/(3*(n-1))-1/a1)/a1;
         c2=2*(pow(n,2)+n+3)/(9*n*(n-1))-(n+2)/(a1*n)+a2/pow(a1,2);
         e2S=c2*st->seg_sites*(st->seg_sites-1)/(pow(a1,2)+a2);
         st->tajd=(st->pairws_dif-st->seg_sites/a1)/sqrt(e1S+e2S);
         //..and Fu & Li's F* (1993)
         vf=2*(n*a1-2*n-2)/((n-1)*(n-2)) + (n-2)/pow(n-1,2)+2/(n-1)*(1.5-(2*(a1+1/n)-3)/(n-2)-1/n); //dn
         vf=(vf+2*(n*n+n+3)/(9*n*(n-1))-2/(n-1)*(4*a2-6+8/n))/(pow(a1,2)+a2);
         uf=(n/(n-1)+(n+1)/(3*n-3)-4/(n*(n-1))+(2*n+2)/pow(n-1,2)*((a1+1/n)-2*n/(n+1)))/a1-vf;
         st->fstar= (st->pairws_dif-(n-1)/n*nu_sing)/sqrt(uf*nu+vf*nu*nu);
      } else { st->seg_sites=0; st->tajd=0.0; st->nuc_div=0.0;}
   } else { //store allelic variance in pairws_dif for now, only 1st locus
      st->pairws_dif=(a1-pow(a2,2)/n)/(n-1);
      for(int i=minmis; i<maxmis; i++) st->mismatch.Add(0);
      for(int i=0; i<st->num_haps; i++) { //save repcount freq of first locus as mismatch distrib
			st->mismatch[*(haps[0][i].hap->begin())-minmis]+=haps[0][i].n;
		}
   }
}

void CalcInterStats(THapArray *h1, THapArray *h2, summary_stats *st, const Mut_Type data_type) {
	int ind, n1=0, n2=0, ss;
	float num_bases=h1[0][0].hap->end()-h1[0][0].hap->begin();
	float het1=1.0, het2=1.0, pdif1=0.0, pdif2=0.0;
	st->nuc_div=1.0;
	st->num_haps=0; 
	st->pairws_dif=0;
	for(int i=0; i<h1->GetItemsInContainer(); i++) n1+=h1[0][i].n;
	st->num_private_haps=0; 
	for(int i=0; i<h2->GetItemsInContainer(); i++) n2+=h2[0][i].n;
	if(data_type==DNA) {
      for(int i=0; i<h1->GetItemsInContainer(); i++) {
         het1-=pow(float(h1[0][i].n)/float(n1),2);
         ind=h2->Find(h1[0][i]);
         if(ind==INT_MAX) { 
            st->num_haps++; //put priv1 in num_haps
            st->nuc_div-=pow(.5*float(h1[0][i].n)/float(n1),2.0);   //pooled heterozygosity (Ht) in nuc_div
         } else
            st->nuc_div-=pow(float(h1[0][i].n)/float(n1)+float(h2[0][ind].n)/float(n2),2)/4.0;
         for(int j=i+1; j<h1->GetItemsInContainer(); j++) { //pairwise difs for ONLY h1
         	ss=SeqDif(h1[0][i].hap->begin(),h1[0][j].hap->begin(),int(num_bases));
         	pdif1+=float(ss*h1[0][i].n*h1[0][j].n);
			}

         for(int j=0; j<h2->GetItemsInContainer(); j++) {
            if(i==0) {
		         for(int k=j+1; k<h2->GetItemsInContainer(); k++) { //pairwise difs for ONLY h2
      		   	ss=SeqDif(h2[0][j].hap->begin(),h2[0][k].hap->begin(),int(num_bases));
         			pdif2+=float(ss*h2[0][j].n*h2[0][k].n);
					}
               ind=h1->Find(h2[0][j]);
               if(ind==INT_MAX) {
                  st->num_private_haps++; //put priv2 in num_private_haps
                  st->nuc_div-=pow(.5*float(h2[0][j].n)/float(n2),2.0);   //continue on with (Ht)
               }
               het2-=pow(float(h2[0][j].n)/float(n2),2);
            }
            ss=SeqDif(h1[0][i].hap->begin(),h2[0][j].hap->begin(),int(num_bases));
            st->pairws_dif+=float(ss*h1[0][i].n*h2[0][j].n);
         }
      }
      st->pairws_dif*=1.0/(float(n1*n2));
      st->hap_div=(het1+het2)/2;  //put mean expected het (Hs-bar) in hap_div
      float mpairws_dif_wn=float(pdif1)/float(n1*(n1-1))+float(pdif2)/float(n2*(n2-1));
      st->tajd=1.0-mpairws_dif_wn/st->pairws_dif; //Fst stored in tajd = (bw-wn)/bw
   } else { //msat and RFLP
      st->num_private_haps=h2->GetItemsInContainer();
      float a91=0, a92=0, Sw, Sb=0, Sbar=0; //for Rst
      for(int i=0; i<h1->GetItemsInContainer(); i++) {
         ind=h2->Find(h1[0][i]);
         //put priv1 in num_haps; 2 in num_private_haps; G^2 in pairws_dif
         if (ind==INT_MAX) {
				st->num_haps++;
				st->pairws_dif+=0.6931472*h1[0][i].n;  //obs*log(obs/pred)=obs*log2 if ai2=0
			} else {
				st->num_private_haps--;
				float pred=n1*(float(h1[0][i].n)/float(n1)+float(h2[0][ind].n)/float(n2))/2.0;
				st->pairws_dif+=(h1[0][i].n==0)?0:h1[0][i].n*log(h1[0][i].n/pred);
			}
         for(int j=i+1; j<h1->GetItemsInContainer(); j++)
				a91+=h1[0][i].n*pow(*(h1[0][i].hap->begin())-*(h1[0][j].hap->begin()),2.0); //SSD(grp1)
			for(int j=0; j<h2->GetItemsInContainer(); j++)
			   Sb+=h1[0][i].n*h2[0][j].n*pow(*(h1[0][i].hap->begin())-*(h2[0][j].hap->begin()),2.0); //SSD(bwn grps)
		}
      for(int i=0; i<h2->GetItemsInContainer(); i++)
         for(int j=i+1; j<h2->GetItemsInContainer(); j++)
				a92+=h2[0][i].n*pow(*(h2[0][i].hap->begin())-*(h2[0][j].hap->begin()),2.0); //SSD(grp2)
		Sw=float(a91)/float(n1*(n1-1))+float(a92)/float(n2*(n2-1));
		Sb*=1/float(n1*n2);
		float nbar=float(n1*n2)/2.0;
		Sbar=(Sw*(nbar-1)+Sb*nbar)/(2*nbar-1);
		st->tajd=(Sbar>0)?(Sbar-Sw)/Sbar:0; //Rst stored in tajd: brutal calc! only for 1st locus
		st->pairws_dif*=2.0; //G^2=2*sum(obs*log(obs/pred)); ALL LOCI
	} //end if DNA or msat
}

TPrior* BayesRead(void *buf, istream *ifs, e_param p, int n, int type) {
   char c=ifs->peek(), str[150], *x,*y, temp[150];
   while(c==' ' || c=='\t' || c=='\n' || c=='\r') { ifs->get(c); c=ifs->peek(); }
   if(c=='{') {
      TPrior *t=new TPrior;
      (*ifs) >> str;
      x=strstr(str+1,":");
      if(!x) { //expression
         x=strtok(str+1,"}");
         if(!x) {printf("\nERROR: No '}' found in formula!"); return NULL;}
         while((y=strstr(x,"-ln("))!=NULL) {sprintf(temp,"~1*%s",y+1); strcpy(x,temp); }
         while((y=strstr(x,"-exp("))!=NULL) {sprintf(temp,"~1*%s",y+1); strcpy(x,temp); }
         while((y=strstr(x,"-("))!=NULL) {sprintf(temp,"~1*%s",y+1); strcpy(x,temp); }
         y=x-1;
         while((y=strstr(y+1,"-"))!=0)
            if(y==x || *(y-1)=='(' || *(y-1)=='+' || *(y-1)=='-' || *(y-1)=='*' || *(y-1)=='/' || *(y-1)=='^' )
                *y='~';
         int vars=numtokens(x,"+-*/^()",1);
         char **extok=(char**)(new double[2*vars-1]);
         tokenize(x,"+-*/^()",extok,1);
         t->formula=CreateTree(extok,vars);
         delete [] extok;
      } else { //prior distribution
         x=strtok(str+1,":"); //"[distribution type]:"
         switch(*x) {
				case 'U': t->dist=uniform; break;
				case 'G': t->dist=gammadist; break;
				case 'E': t->dist=exponential; break;
				case 'N': t->dist=normal; break;
				case 'M': t->dist=geometric; break;
				default: printf("\nERROR: What kind of distribution is \"%s\"?",x); return NULL; break;
			}
			if(t->dist==exponential || t->dist==geometric) x=strtok(NULL,"}"); //just one parameter
			else x=strtok(NULL,",");
         if(!x) {printf("\nERROR: Could not read prior distribution!"); return NULL;}
         t->min=str2double(x);
         if(t->dist==exponential || t->dist==geometric) t->max==0.0; //just one parameter
			else {
				x=strtok(NULL,"}");
         	if(!x) {printf("\nERROR: Could not read prior distribution!"); return NULL;}
         	t->max=str2double(x);
			}
			if(buf) switch(type) {
            case 1: *(int*)buf=(int)t->min; break;
            case 2: *(long*)buf=(long)t->min; break;
            case 3: *(float*)buf=(float)t->min; break;
            case 4: *(double*)buf=t->min; break;
			}
      }
      t->param=p; t->num=n; t->type=type; t->curval=t->min;
      return t;
   } else {
      if(buf) {
         switch(type) {
            case 1: (*ifs) >> *(int*)buf; break;
            case 2: (*ifs) >> *(long*)buf; break;
            case 3: (*ifs) >> *(float*)buf; break;
            case 4: (*ifs) >> *(double*)buf; break;
         }
      }
      return NULL;
   }
}

double str2double(char *x) {
	double y=0,z=1.0,v;
	int dec=0;
	if(!x) {
		printf("\nWARNING: Tried to convert NULL string to float");
		return 0.0;
	}
	while(*x) {
		if(*x=='-' || *x=='~') z=-1.0;
		else if (*x=='.') dec=1;
		else if (*x>='0' && *x<='9') {
			v=(double)((*x)-'0');
			if(dec==0) y=y*10+v;
			else y=y+v/pow(float(10),float(dec++));
		}
		else if(*x==' ' || *x=='\t') y*=1.0;
		else if(*x=='T') y=-1;
		else {
			printf("\nWARNING: Tried to convert non-numeric string");
			break;
		}
		x++;
	}
	return z*y;
}

//Returns the cut string in toked
int tokenize(char *s, char *sep, char **toked, bool include_tok) {
   int n=0; 
   toked[0]=s;
   while(*s) {
      if(strchr(sep,*s)) {
         char *x=toked[n];
         toked[n]=new char[s-x+1];
         if(s-x>0) strncpy(toked[n],x,s-x);
         toked[n++][s-x]='\0';
         if(include_tok) {
            toked[n]=new char[2];
            sprintf(toked[n++],"%c\0",*s);
         }
         toked[n]=s+1;
      }
      s++;
   }
   char *x=toked[n];
   toked[n]=new char[s-x+1];
   if(*x) strcpy(toked[n++],x);
   else *(toked[n])='\0';
   return n;
}

void param2str(char *x, e_param p) {
   switch(p) {
      case deme_size: strcpy(x,"Deme Size"); break;
      case sample_age: strcpy(x,"Sample Age"); break;
      case sample_stat_grp: strcpy(x,"Sample StatGrp"); break;
      case growth_rate: strcpy(x,"Growth Rate"); break;
      case mig_rate: strcpy(x,"MigrMatrix"); break;
      case event_time: strcpy(x,"Event Time"); break;
      case event_src: strcpy(x,"Event Source"); break;
      case event_sink: strcpy(x,"Event Sink"); break;
      case event_migrants: strcpy(x,"Event Migrants"); break;
      case event_size: strcpy(x,"Event Size"); break;
      case event_growth: strcpy(x,"Event Growth"); break;
      case event_mig_matr: strcpy(x,"Event MigrMatrix"); break;
      case pmut_rate: strcpy(x,"Mutation Rate"); break;
      case abstract: strcpy(x,"Abstract"); break;
      default: strcpy(x,"[UNKNOWN TYPE]"); break;
   }
}

char dist2str(e_dist d) {
	switch(d) {
		case expression: return '#'; break;
		case uniform: return 'U'; break;
		case normal: return 'N'; break;
		case gammadist: return 'G'; break;
		case exponential: return 'E'; break;
		case geometric: return 'M'; break;
		default: return '!'; break;
	}
}

//counts the number of strings in s that can be made by cutting at characters in sep
int numtokens(char *s, char *sep, bool echo) {
   int n=1; //no cuts=1 string, 1 cut=2 str
   while(*s) {
      if(strchr(sep,*s)) n++;
      s++;
   }
   return n;
}

//splices the expression string into a reverse-polish evaluation tree
symb *CreateTree(char **toked, int vars) { 
   symb *temp, *root=NULL;
   int plvl=0;
   for(int i=0;i<2*vars-1;i++) {
      temp=new symb;
      temp->parens=plvl;
      if (strlen(toked[i])==0) continue;
      else if(strcmp(toked[i],"+")==0) temp->tp=PLUS;
      else if(strcmp(toked[i],"-")==0) temp->tp=MINUS;
      else if(strcmp(toked[i],"*")==0) temp->tp=TIMES;
      else if(strcmp(toked[i],"/")==0) temp->tp=DIVIDE;
      else if(strcmp(toked[i],"^")==0) temp->tp=POWER;
      else if(strcmp(toked[i],"ln")==0) { temp->tp=FUNCTION; temp->val=1.0; }
      else if(strcmp(toked[i],"exp")==0) { temp->tp=FUNCTION; temp->val=2.0; }
      else if(strcmp(toked[i],"(")==0) { plvl++; continue; }
      else if(strcmp(toked[i],")")==0) { plvl--; continue; }
      else if(*(toked[i])=='[') {
         temp->tp=PRIOR;
         char *x=toked[i]+1;
         while(*x!=']') *(x-1)=*(x++);
         *(x-1)='\0';
         temp->val=str2double(toked[i]);
      } else { 
         temp->tp=VALUE; 
         temp->val=str2double(toked[i]);
      }
      delete toked[i];
      if(!root) {
         root=temp;
         if(root->tp==TIMES || root->tp==DIVIDE || root->tp==POWER) {
            printf("Syntax error 1: Cannot start an expression with any of these: */^"); return 0;
         } else if (root->tp==PLUS || root->tp==MINUS) {//handles -3 or +2
            temp=new symb;
            temp->parens=plvl;
            root->l=temp;
         }
      } else if (temp->tp!=VALUE && temp->tp!=PRIOR) {
         symb *t2=root, *t3=NULL;
         while(1) {
            if(!t2) break;
            else if(t2->parens<temp->parens) { t3=t2; t2=t2->r; }
            else if(t2->parens==temp->parens && t2->tp < temp->tp) { t3=t2; t2=t2->r; }
            else break;
         }
         if(t2) { 
            temp->l=t2;
            (t3)? t3->r=temp : root=temp;
         } else { //root
            if(temp->tp==PLUS || temp->tp==MINUS || temp->tp==FUNCTION) { 
               t3->r=temp;
               t2=new symb;
               t2->parens=plvl;
               temp->l=t2;
            } else {
               printf("Syntax error 2");
               return 0;
            }
         }
      } else if (root!=temp) { 
         symb *t2=root;
         while(t2->r) t2=t2->r;
         t2->r=temp;
      } else {
         printf("Syntax error 3");
         return 0;
      }
   }
   return root;
}

double rgamma(int k, double theta) {
	double a=0.0;
	for(int i=0; i<k; i++) a+=-log(mtr());
	return(theta*a);
}

void TPrior::pick(void *memloc) {
   if(!type || param==none || !memloc) return;
   float y;
   if(!formula) {
      y=mtr();
      switch(dist) {
			case uniform: curval=y*(max-min)+min; break;
			case normal: curval=sqrt(-2*log(y))*cos(6.283185*mtr())*max+min; break;
			case gammadist: curval=rgamma((int)min,max); break;
			case exponential: curval=-log(y)*min; break;
			case geometric: curval=ceil(log(y)/log(min)); break;
		}
      set(memloc);
   } else printf("\nERROR: Trying to evaluate a formula as a prior.");
}
