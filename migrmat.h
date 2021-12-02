#ifndef _MIGRATION_H
#define _MIGRATION_H

#include "cstring.h"
#include <fstream>
#include <iomanip>
#include <float.h>
#include <cmath>

using namespace std;
//typedef my_string string;

//A two integers array for storing mîgration events coordinates
typedef int indices[2];
typedef float MatElemType;

enum Matrix_Type {CIRCULAR_STEPPING_STONE, STEPPING_STONE_2D, STEPPING_STONE_3D, ISLAND, RANDOM};

	// A class to store and manage the migration flux between a set of populations
	// By default it is a square matrix
	// It also contains two arrays containing the net emigration and immigration
	// flux of each population
	//The migration matrix must be build such the the total emigration rate should
	//be equal to the total immigration rate
class TMigrationMatrix {
	public:
		TMigrationMatrix() {
         _size=0;
   		_Pmigr=0;
         _data=NULL;
      };
      TMigrationMatrix(const char* name);
      TMigrationMatrix(const int& n) {
         _size=n;  
   		_Pmigr=0;
         _data=NULL;
      	allocate_matrix(n);
         EmigrFile[0]='\0';
      };
      TMigrationMatrix(const int& size, const Matrix_Type& MT, const float& m);
      TMigrationMatrix(const TMigrationMatrix& MM) {
      	if (this!=&MM) *this=MM;
      }
      ~TMigrationMatrix();
		int read_from_file();
      int read_from_file(const char* filename);
      int de_allocate_matrix();
      int allocate_matrix(const int n) {
      	de_allocate_matrix();
         _size=n;
         allocate_matrix();
         return _size;
      };
      int display_matrix(ostream&) const ;
		MatElemType& operator()(const int& coord_X,const int& coord_Y) const {
         if ((coord_X < _size) && (coord_Y < _size ) && (coord_X>=0)  && (coord_Y>=0))
         	return _data[coord_X][coord_Y];
         cout 	<< "TMigrationMatrix : Array index out of range ("
         		<< coord_X << ", " << coord_Y << ")" << endl;
         return _data[0][0]; //Should be zero  
      };
                                                  
      int compute_total_migration_prob();
      const double& migr_prob() {return _Pmigr;}
      
      //Returns a pair of population indices recording a migration event:
      //The first index represents the source population
      //The second index represents the sink population
      int get_migration_indices(indices& indx);
      int size() const {return _size;}

      TMigrationMatrix& operator=(const TMigrationMatrix& MM);
      //Loro_16_9_99
      int operator==(const TMigrationMatrix& MM);

      //set all matrix elements equal to m
      int set_elems(const MatElemType& m);


	protected :
		MatElemType  **_data;   //Stores the migration probabilties
		int    		_size;      //The size of the matrix (number of populations)
      double 		_Pmigr;     //The total probability of a migration

		char EmigrFile[200]; //A file from which to read the migration matrix
		friend istream& operator>>(istream& is, TMigrationMatrix& S);
      
		friend ostream& operator<<(ostream& os, const TMigrationMatrix& S);

		int allocate_matrix();
};

extern ostream& operator<<(ostream& os, const TMigrationMatrix& S);
extern istream& operator>>(istream& is, TMigrationMatrix& S);

#endif
