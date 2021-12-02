#ifndef _COALMAIN_H_
#define _COALMAIN_H_

#include "migrmat.h"

//------------------------------------------------------------------------------

template <class T> class myIntelligentVector {
   T *_data;
   int _size;
   double _mean, _sd;
public:
   myIntelligentVector(int);
   ~myIntelligentVector( ) { delete[] _data; }
   T& operator[ ] (int i) { return _data[i]; }
   void compute_moments();
   double mean() {return _mean;}
   double sd() {return _sd;}
};

// Note the syntax for out-of-line definitions.
template <class T>
myIntelligentVector<T>::myIntelligentVector(int n) {
   try {
   	_data = new T[n];
   	_size = n;
   }
   catch (...) {
   	if (_data) delete[] _data;
      _data=NULL;
      _size=0;
      cout 	<< "Unable to allocate memory for myIntelligentVector"
            << endl;
   }
};

template <class T> void
myIntelligentVector<T>::compute_moments() {
   double s=0, s2=0;
   _mean=_sd=0;
   for (int i=0; i<_size; ++i) {
    	s+=_data[i];
      s2+=_data[i]*_data[i];
   }
   _mean=s/_size;
   if (_size>1)
   	_sd=sqrt((s2-s*s/_size)/(_size-1));
};

//------------------------------------------------------------------------------
//A class derived from TMigrationMatrix to store means
class TMeanMat : public TMigrationMatrix {

   long _num_updates;

	public :
   	TMeanMat() : TMigrationMatrix() {_num_updates=0;};
      TMeanMat(const int& size, const float& m) :
      				TMigrationMatrix(size, ISLAND, m) {
      	_num_updates=0;
      };

      ~TMeanMat() {};
      int add_elem(const int& i, const int& j, MatElemType& e ) {
      	if (!_data) return 0;
         _data[i][j]=e;
         return 1;
      }
      int compute_mean(const int& n);
      int compute_sd(const int& n,  const TMeanMat& mean);
      int update_with(const TMigrationMatrix& MM);
		int update_with_square(const TMigrationMatrix& MM);

};

#endif
