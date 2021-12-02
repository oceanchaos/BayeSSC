#ifndef _MY_Cstring_HPP_
#define _MY_Cstring_HPP_


#include <string.h>   //for functions strlen strcpy strcat
#include <iostream>  //for << operator
//#include <fstream.h>   //for read_token
#include <stdio.h>    //for sprint
//#include <strstream.h>

using namespace std;
const int NPOS=-1;

class my_string{
	private:
   	char* _data;
      int _del_inc;
      int _nb_char;
      int _tot_space;
      char* _reserve;
   protected:
   	int tot_space(){return _tot_space;};
   public:

   	my_string(){
      	 _data=new char[10];
         _data[0]='\0';
         _del_inc=10;
         _tot_space=10;
      	 _nb_char=0;
         _reserve=NULL;
      };

   	my_string(const char& c){
      	 _data=new char[10];
         _data[0]=c;
         _data[1]='\0';
         _del_inc=10;
         _tot_space=10;
      	 _nb_char=1;
         _reserve=NULL;
      };


      my_string(int i){
      	 _data=new char[i+2];
         _data[0]='\0';
         _del_inc=10;
         _tot_space=i+2;
      	 _nb_char=0;
         _reserve=NULL;
      };

      my_string(const my_string& x) {
      	 _data=new char[10];
         _data[0]='\0';
         _del_inc=10;
         _tot_space=10;
      	 _nb_char=0;
      	 *this=x;
         _reserve=NULL;
      };
      my_string(char* x) {
         _nb_char=strlen(x);
         _del_inc=10;
         if(_nb_char<10000){
            _data=new char[_nb_char+2];
            strcpy(_data,x);
            _tot_space=_nb_char+2;
         }
         else {
      		_data=new char[10];
         	_data[0]='\0';
         	_del_inc=10;
         	_tot_space=10;
      		_nb_char=0;
         }
      	_reserve=NULL;
      };

      ~my_string(){if(_data) delete[] _data; if(_reserve) delete[] _reserve;};

      char* c_str() const {return _data;};
      int length() const {return _nb_char;};
      void Lower();
      void Upper();
      void read_token(istream& is);
      //void read_token(strstream& is);
      my_string& operator=(const my_string& x);
      my_string& operator=(char* x);
      my_string& operator=(const char& x);
      my_string& operator=(int x);
      my_string& operator=(float x);
      my_string& operator+=(const my_string& x);
      my_string& operator+=(char* x);
      my_string& operator+=(char x); // new
      char operator[](int pos) const {return _data[pos];};  //new
      char&  operator[](int pos) {return _data[pos];};  //new
      friend my_string operator+(const my_string& x, const my_string& y);
      friend my_string operator+(const my_string& x ,char* y);
      friend ostream& operator<<(ostream& os, const my_string& x);
      friend istream& operator>>(istream& is, my_string& s);
      int operator==(const my_string& x) const{return !(strcmp(_data, x.c_str()));};
      int operator==(char* x) const {return !(strcmp(_data, x));};
      int operator!=(const my_string& x) const {return (strcmp(_data, x.c_str()));};
      int operator!=(char* x)  const {return (strcmp(_data, x));};
      int operator<(const my_string& x)const {return ((strcmp(_data, x.c_str()))< 0) ;};  //new
      int operator>(const my_string& x)const {return ((strcmp(_data, x.c_str()))> 0) ;};  //new
      //added since 21_3_97:
      void to_lower();      //converts my_string to lower case
      void to_upper();      //converts my_string to upper case
      bool read_to_delim(istream& is, char delim); //read to delimiter, return "true" if found before eof or "\n"
      void read_to_delim(istream& is); //read to delimiter
      my_string extract_sub_str(char delim); //extracts until the next occurence of delim or end of string //new 11_11_97
      my_string& remove( int pos, int n );
      int find_first_of( const char& s, int pos ) const;
      int find_first_of( const char& s) const;
      int read_line(istream& is);
      int is_null();
      int contains(const char* pat) const;
      int contains(char pat) const;
		int contains(const my_string& s) const;
      void assign( const my_string& s ) {*this=s;};
      int rfind( const my_string& s ); //return the position of patter s if not found return -1
      void rm_path(); //stef_28_3_98 removes path info from a file name

};

//Christian: from coalmain.cpp
my_string extract_path(const my_string& s);
my_string remove_extension(const my_string& s);
void my_strrev(char * str);


#endif







