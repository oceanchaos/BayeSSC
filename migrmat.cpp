//------------------------------------------------------------------------------
// (c) L. Excoffier , LGB, University of Geneva, August 1995
//
//------------------------------------------------------------------------------

#include "migrmat.h"
#include "public.h"
//------------------------------------------------------------------------------
//Initialization is done by reading from a file
TMigrationMatrix::TMigrationMatrix(const char* name) {
	_data=NULL;
   _Pmigr=0.0;
   read_from_file(name);
};
//------------------------------------------------------------------------------
TMigrationMatrix::TMigrationMatrix(	const int& size, const Matrix_Type& MT,
												const float& m) {
	_data=NULL;
   _Pmigr=0.0;
   EmigrFile[0]='\0';

   if (MT==STEPPING_STONE_3D) {_size=size*size; }
   else _size=size;
   
	allocate_matrix();

   	int i,j;
      int grid_size=size;
      int row, col;

	switch (MT) {
   //Island migration model
   	case ISLAND :  //Test ok
      	for (i=0; i<size; ++i) {
         	for (j=0; j<size; ++j) {
            	if (i!=j) _data[i][j]=m;
            }
         }
      	break;

   	//Circular stepping stone model
   	case CIRCULAR_STEPPING_STONE :  //test ok
      	for (i=0; i<size-1; ++i) {
         	_data[i][i+1]=m;
         	_data[i+1][i]=m;
         };
         _data[size-1][size-2]=m;
         _data[size-1][0]=m;
         _data[0][size-1]=m;
         break;

   	//2D stepping stone model, not circular
   	case STEPPING_STONE_2D :    //Test ok
      	for (i=0; i<size-1; ++i) {
         	_data[i][i+1]=m;
         	_data[i+1][i]=m;
         };
         _data[size-1][size-2]=m;
         break;

      //3D stepping-stone model: the number of population should be arranged on
      //a size by size square grid
      //This is analogous to a torous
      //The following procedure is not absolutely optimal, but I a sure all migration
      //rates are properly set
   	case STEPPING_STONE_3D :
      	for (i=0; i<grid_size; ++i) {
      		for (j=0; j<grid_size; ++j) {
            	row=(i*grid_size)+j;
               //There are four migration rates to change
               //West
               if (j==0) col=row+grid_size-1; //Go to the end of line
               else col=row-1;
               _data[row][col]=m;
               //East
               if (j==(grid_size-1)) col=row-grid_size+1; //Go to the beginning of line
               else col=row+1;
               _data[row][col]=m;
               //North
               if (i==0) col=(grid_size-1)*grid_size+j;  //Go to the last line
               else col=row-grid_size;
               _data[row][col]=m;
               //South
               if (i==grid_size-1) col=j;  //Go to the first line
               else col=row+grid_size; 
               _data[row][col]=m;
         	}
      	}
      	break;

      //Random island migration rate model:
      //The migration rates are assumed of the order of m but are alowed to vary.
      //The expectation of the migration rate is m, and can vary between 0 an 2m
   	case RANDOM :   //test ok
			for (i=0; i<size; ++i) {
         	for (j=0; j<size; ++j) {
            	if (i!=j) _data[i][j]=2.0*m*mtr();
            }
         }
      	break;
   	default : //Island migration model
			for (i=0; i<size; ++i) {
         	for (j=0; j<size; ++j) {
            	if (i!=j) _data[i][j]=m;
            }
         }
      break;
   };

};
//------------------------------------------------------------------------------
TMigrationMatrix&
TMigrationMatrix::operator=(const TMigrationMatrix& MM) {

	de_allocate_matrix();
 	_size=MM._size;
   _Pmigr=MM._Pmigr;
   strcpy(EmigrFile,MM.EmigrFile);

   //Allocating memory to the matrix
   try {
   	// TMigrationMatrix matrix allocation
      _data = new MatElemType*[_size];      // STEP 1: SET UP THE ROWS.
      for (int j=0; j<_size; ++j)  {
      	_data[j] = new MatElemType[_size];  // STEP 2: SET UP THE COLUMNS
      }
   }
   catch (...) {
   	de_allocate_matrix();
      cout << "TMigrationMatrix::operator=(): memory error !" << endl;
      return *this;
   }

   //Copying elements
   for (int i=0; i<_size; ++i)
      for ( int j=0; j<_size; ++j)
         _data[i][j]=MM._data[i][j];
   return *this;
};
//------------------------------------------------------------------------------
//Loro_16_9_99
int TMigrationMatrix::operator==(const TMigrationMatrix& MM){
	if (_size!=MM._size) return 0;
   if (!_data || !MM._data) return 0;
   for (int i=0; i<_size; ++i) {
   	for (int j=0; j<_size; ++j) {
      	if (_data[i][j]!=MM._data[i][j]) return 0;
      }
   }
   return 1;
}
//------------------------------------------------------------------------------
//Destructor
TMigrationMatrix::~TMigrationMatrix() {
   de_allocate_matrix();
}
//------------------------------------------------------------------------------
int
TMigrationMatrix::read_from_file()
{
	ifstream is(EmigrFile);
	if (!is)
	{
		cout << "TMigrationMatrix::read_from_file(): bad stream" << endl;
		return 0;
	}
	else
	{
		is >> *this;
		return 1;
	}
}
//------------------------------------------------------------------------------
int
TMigrationMatrix::read_from_file(const char* filename) {
   strcpy(EmigrFile,filename);
   read_from_file();
   return 1;
}
//------------------------------------------------------------------------------
int
TMigrationMatrix::allocate_matrix() {

   if (_data) de_allocate_matrix();

  		try {
			// TMigrationMatrix matrix allocation
         _data = new MatElemType*[_size];      // STEP 1: SET UP THE ROWS.
         for (int j=0; j<_size; ++j)
            _data[j] = new MatElemType[_size];  // STEP 2: SET UP THE COLUMNS
		}
		catch (...) {
			de_allocate_matrix();
			cout << "TMigrationMatrix::allocate_matrix(): memory error !" << endl;
         return 0;
		}
      // Initialisation
      for ( int i=0; i<_size; ++i)
      for ( int j=0; j<_size; ++j)
         _data[i][j]=0.0;
   
   return 1;
}
//------------------------------------------------------------------------------
int
TMigrationMatrix::de_allocate_matrix() {

   if (_data) {
      for (int i=_size; i>0;  --i)  {
         if (_data[i-1]) {
            delete[] _data[i-1];
         }
      }
      delete[] _data;
      _data=NULL;
      _size=0;
   }
   return 1;
}
//------------------------------------------------------------------------------
int
TMigrationMatrix::display_matrix(ostream& os) const {
   if (_data)
   for (int i = 0; i < _size; ++i) {
      for (int j = 0; j < _size; ++j)
             os 	<< setiosflags (ios::showpoint | ios::fixed | ios::right)
						<< setw(10)
						<< setprecision(4)
             		<< _data[i][j] << " ";
      os << "\n";
   }
   else os << "Matrix contains no data" << endl;
   return 1;
}
//------------------------------------------------------------------------------
int TMigrationMatrix::compute_total_migration_prob() {

   //Compute the total sum of the migration probabilities
  	_Pmigr=0.0;
   for (int i=0; i<_size; ++i) {
   	for (int j=0; j<_size; ++j) {
      	_Pmigr+=_data[i][j];
   	}
   }
   if (_Pmigr==0.0) return 0; //There is no possible migration then...

   //Set relative migration probabilities
   for (int i=0; i<_size; ++i) {
   	for (int j=0; j<_size; ++j) {
      	_data[i][j]/=_Pmigr;
   	}
   }
   return 1;
};
//------------------------------------------------------------------------------
//Returns a pair of population indices recording a migration event:
//The first index represents the source population
//The second index represents the sink population
//!!! Assumes that compute_total_migration_prob() has been called already
int TMigrationMatrix::get_migration_indices(indices& indx){
   //Second, get the index of the source  and sink populations
   double randnum=mtr();
   double tot=0.0;
   for (int i=0; i<_size ; ++i) {
   	for (int j=0; j<_size && tot<randnum; ++j) {
   		tot+=_data[i][j];
      	indx[0]=i;
         indx[1]=j;
      }
   }
 	return 1;
};
//------------------------------------------------------------------------------
istream& operator>>(istream& is, TMigrationMatrix& B) {
   float mat_elem;
   if (!is) {
		cout << "operator>>(ifstream& is, TMigrationMatrix& S): bad stream " << endl;
		return is;
	}
	//Begins by reading the header containing the structure of the matrix
   is >> B._size;

	if (is.bad() )
	{
		cout 	<< "operator>>(istream& is, TMigrationMatrix& S): "
      		<< "bad stream" << endl;
		return is;
	}
	//Create matrix in memory
	B.allocate_matrix();

 	register int cy;
	for (register int cx=0; cx < B._size; ++cx)
	{
		for (cy = 0; cy<B._size; ++cy)
		{
         is >> mat_elem;
         B(cx,cy)=mat_elem;
   
			if (!is)	break;
		};
	}
	if (!is)
		cout << "operator>>(istream& is, TMigrationMatrix& M): bad stream" << endl;
	return is;
}
//------------------------------------------------------------------------------
int TMigrationMatrix::set_elems(const MatElemType& m) {
	if (!_data) return 0;
   for (int i=0; i<_size; ++i)
   	for (int j=0; j<_size; ++j) _data[i][j]=m;
   return 1;
}

ostream& operator<<(ostream& os, const TMigrationMatrix& S) {
	S.display_matrix(os);
   os << endl;
   return os;
};
