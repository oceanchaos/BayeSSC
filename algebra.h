#ifndef _ALGEBRA_H_
#define _ALGEBRA_H_

#include <math.h>

enum op_type {PLUS, MINUS, TIMES, DIVIDE, POWER, VALUE, PRIOR, FUNCTION};

class symb {
   public:
      op_type tp;
      double val; //if op_type=FUNCTION, 1=ln() 2=exp
      int parens;
      symb *r, *l;

      symb() {r=NULL; l=NULL; val=0; parens=0; tp=VALUE; }
      ~symb() {
         if(r) delete r;
         if(l) delete l;
         printf("x");
      }
      void print(int lvl=0) {
         printf("\n"); for(int i=0; i<lvl; i++) printf("  ");
         switch(tp) {
            case PLUS: printf("+"); break;
            case MINUS: printf("-"); break;
            case TIMES: printf("*"); break;
            case DIVIDE: printf("/"); break;
            case POWER: printf("^"); break;
            case VALUE: printf("%g",val); break;
            case PRIOR: printf("[%d]",(int)val); break;
            case FUNCTION:
               switch((int)val) {
                  case 1: printf("ln"); break;
                  case 2: printf("exp"); break;
               } break;
         }
         printf(" (%d)",parens);
         if(l) l->print(lvl+1);
         if(r) r->print(lvl+1);
      }
      symb* CopyTree() {
         symb *y=new symb;
         *y=*this;
         if(l) y->l=l->CopyTree();
         if(r) y->r=r->CopyTree();
         return y;
      }
   double solve(double *prs) {
	   switch (tp) {
	      case PLUS: return l->solve(prs)+r->solve(prs); break;
	      case MINUS: return l->solve(prs)-r->solve(prs); break;
	      case TIMES: return l->solve(prs)*r->solve(prs); break;
	      case DIVIDE: return l->solve(prs)/r->solve(prs); break;
	      case POWER: return pow(l->solve(prs),r->solve(prs)); break;
	      case VALUE: return val; break;
	      case PRIOR: return (int(val)==-1)? *prs :prs[int(val)-1]; break;
	      case FUNCTION: //for functions, the left branch is a dummy
	         switch(int(val)) {
	            case 1: return log(r->solve(prs)); break;
	            case 2: return exp(r->solve(prs)); break;
	            default: printf("\nERROR: Unknown function!");
	         } break;
	      default:
	         printf("\nERROR: Trying to solve invalid type!");
	         return 0; break;
	   }
		return 0;
	}

};



#endif
