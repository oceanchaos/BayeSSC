#include "cstring.h"
//for debug
//#include "a_glbdef.hpp"
//#include "a_glbobj.hpp"
//end debug

//------------------------------------------------------------------------------
//diverse functions
void
my_string::Lower(){
	int len=strlen(_data);
  	for(int i=0; i<len; i++){
   	if(_data[i]=='A') _data[i]='a';
      if(_data[i]=='B') _data[i]='b';
      if(_data[i]=='C') _data[i]='c';
      if(_data[i]=='D') _data[i]='d';
      if(_data[i]=='E') _data[i]='e';
      if(_data[i]=='F') _data[i]='f';
      if(_data[i]=='G') _data[i]='g';
      if(_data[i]=='H') _data[i]='h';
      if(_data[i]=='I') _data[i]='i';
      if(_data[i]=='J') _data[i]='j';
      if(_data[i]=='K') _data[i]='k';
      if(_data[i]=='L') _data[i]='l';
      if(_data[i]=='M') _data[i]='m';
      if(_data[i]=='N') _data[i]='n';
      if(_data[i]=='O') _data[i]='o';
      if(_data[i]=='P') _data[i]='p';
      if(_data[i]=='Q') _data[i]='q';
      if(_data[i]=='R') _data[i]='r';
      if(_data[i]=='S') _data[i]='s';
      if(_data[i]=='T') _data[i]='t';
      if(_data[i]=='U') _data[i]='u';
      if(_data[i]=='V') _data[i]='v';
      if(_data[i]=='W') _data[i]='w';
      if(_data[i]=='X') _data[i]='x';
      if(_data[i]=='Y') _data[i]='y';
      if(_data[i]=='Z') _data[i]='z';
   }
};

void
my_string::Upper(){ //convert to upper case
	int len=strlen(_data);
   for(int i=0; i<len; i++){
   	if(_data[i]=='a') _data[i]='A';
      if(_data[i]=='b') _data[i]='B';
      if(_data[i]=='c') _data[i]='C';
      if(_data[i]=='d') _data[i]='D';
      if(_data[i]=='e') _data[i]='E';
      if(_data[i]=='f') _data[i]='F';
      if(_data[i]=='g') _data[i]='G';
      if(_data[i]=='h') _data[i]='H';
      if(_data[i]=='i') _data[i]='I';
      if(_data[i]=='j') _data[i]='J';
      if(_data[i]=='k') _data[i]='K';
      if(_data[i]=='l') _data[i]='L';
      if(_data[i]=='m') _data[i]='M';
      if(_data[i]=='n') _data[i]='N';
      if(_data[i]=='o') _data[i]='O';
      if(_data[i]=='p') _data[i]='P';
      if(_data[i]=='q') _data[i]='Q';
      if(_data[i]=='r') _data[i]='R';
      if(_data[i]=='s') _data[i]='S';
      if(_data[i]=='t') _data[i]='T';
      if(_data[i]=='u') _data[i]='U';
      if(_data[i]=='v') _data[i]='V';
      if(_data[i]=='w') _data[i]='W';
      if(_data[i]=='x') _data[i]='X';
      if(_data[i]=='y') _data[i]='Y';
      if(_data[i]=='z') _data[i]='Z';
   }  
};

void
my_string::to_lower(){
  Lower();
};

void
my_string::to_upper(){ //convert to upper case
	Upper();
};
//Loro_29_8_98 There are bugs in this routine: do not use or debug it first
void
my_string::read_token(istream& is){ //not tested
//1) clean:
	if(_data) delete[] _data;
   _data=new char[10];
   _data[0]='\0';
   _del_inc=10;
   _tot_space=10;
   _nb_char=0;
   _reserve=NULL;
//2) skip blanks
	char current=' ';
   for(;current==' ' || current=='\t' || current=='\n' && is;){
   	 is.get(current);  //stop at first non blank char
   }
//3)take as long no blanks
   for(;current!=' ' && current!='\t' && current!='\n';){
   	if(_nb_char>=_tot_space-1){
      	int len=_tot_space+_del_inc;
   		if(_reserve) delete[] _reserve;
      	_reserve=new char[(strlen(_data)+2)];     //make spcace for buffer
      	if(!_reserve) return ;               //memory out
      	strcpy(_reserve, _data);                  //copy data to buffer
   		delete[] _data;
      	_data=new char[len+2];                    //make enought space for catenate
      	_tot_space=len+2;
      	strcpy(_data, _reserve);                  //copy bac
      	delete[] _reserve; _reserve=NULL;         //clean buffer;
         _tot_space=len;    								//store new len
      }
      _data[_nb_char]=current;
      _nb_char++;
      _data[_nb_char]='\0';
      is.get(current);
   }
   _data[_nb_char]='\0';
   if(current=='\n') is.putback(current);     //stef_8_5_98   read token is
                                              //supposed not to 'eat' newline

};

/*void
my_string::read_token(strstream& is){ //not tested
//1) clean:
	if(_data) delete[] _data;
   _data=new char[10];
   _data[0]='\0';
   _del_inc=10;
   _tot_space=10;
   _nb_char=0;
   _reserve=NULL;
//2) skipp blanks
	char current=' ';
   for(;current==' ' || current=='\t' || current=='\n';){
   	 is.get(current);  //stop at first not blank char
   }
//3)take as long no blanks
   for(;current!=' ' && current!='\t' && current!='\n';){
   	if(_nb_char>=_tot_space-1){
      	int len=_tot_space+_del_inc;
   		if(_reserve) delete[] _reserve;
      	_reserve=new char[(strlen(_data)+2)];     //make spcace for buffer
      	if(!_reserve) return ;               //memory out
      	strcpy(_reserve, _data);                  //copy data to buffer
   		delete[] _data;
      	_data=new char[len+2];                    //make enought space for catenate
      	_tot_space=len+2;
      	strcpy(_data, _reserve);                  //copy bac
      	delete[] _reserve; _reserve=NULL;         //clean buffer;
         _tot_space=len;    								//store new len
      }
      _data[_nb_char]=current;
      _nb_char++;
      _data[_nb_char]='\0';
      is.get(current);
   }
   _data[_nb_char]='\0';



};*/


//------------------------------------------------------------------------------
//copy operators
my_string&
my_string::operator=(const char& x){

   if(_tot_space < 4){
   	delete[] _data;
      _data=new char[12];
      _tot_space=12;
   }
   _nb_char=1;
   _data[0]=x;
   _data[1]='\0';
   return *this;
};

my_string&
my_string::operator=(const my_string& x){
	int len=strlen(x.c_str());
   if(len > 10000){
		return *this;
   }
   if(_tot_space < len+2){
   	delete[] _data;
      _data=new char[len+2];
      _tot_space=len+2;
   }
   _nb_char=len;
   strcpy(_data, x.c_str());
   return *this;
};

my_string&
my_string::operator=(char* x){
	int len=strlen(x);
   if(len > 10000){
		return *this;
   }
   if(_tot_space < len+2){
   	delete[] _data;
      _data=new char[len+2];
      _tot_space=len+2;
   }
   _nb_char=len;
   strcpy(_data, x);
   return *this;
};

my_string&
my_string::operator=(int x){

   if(_tot_space < 20){
   	delete[] _data;
      _data=new char[22];
      _tot_space=22;
   }
   sprintf(_data,"%d",x);
    _nb_char=strlen(_data);
   return *this;
};

my_string&
my_string::operator=(float x){
   if(_tot_space < 40){
   	delete[] _data;
      _data=new char[42];
      _tot_space=42;
   }
   sprintf(_data,"%f",x);
    _nb_char=strlen(_data);
   return *this;
};
//------------------------------------------------------------------------------
//catenation operators
my_string&
my_string::operator+=(const my_string& x){
	int len=strlen(x.c_str())+strlen(_data);
   if(len > 10000){
		return *this;
   }
   if(_tot_space < len+2){
   	if(_reserve) delete[] _reserve;
      _reserve=new char[(strlen(_data)+2)];     //make spcace for buffer
      if(!_reserve) return *this;               //memory out
      strcpy(_reserve, _data);                  //copy data to buffer
   	delete[] _data;
      _data=new char[len+2];                    //make enought space for catenate
      _tot_space=len+2;
      strcpy(_data, _reserve);                  //copy bac
      delete[] _reserve; _reserve=NULL;         //clean buffer;
   }
   _nb_char=len;
   if (_data) strcat(_data, x.c_str());
   return *this;
};

my_string&
my_string::operator+=(char* x){
	int len=strlen(x)+strlen(_data);
   if(len > 10000){
		return *this;
   }
   if(_tot_space < len+2){
   	if(_reserve) delete[] _reserve;
      _reserve=new char[(strlen(_data)+2)];     //make spcace for buffer
      if(!_reserve) return *this;               //memory out
      strcpy(_reserve, _data);                  //copy data to buffer
   	delete[] _data;
      _data=new char[len+2];                    //make enought space for catenate
      _tot_space=len+2;
      strcpy(_data, _reserve);                  //copy bac
      delete[] _reserve; _reserve=NULL;         //clean buffer;
   }
   _nb_char=len;
   if (_data) strcat(_data, x);
   return *this;
};


my_string&
my_string::operator+=(char x){
	int len=strlen(_data)+1;
   if(len > 10000){
		return *this;
   }
   if(_tot_space < len+2){
   	if(_reserve) delete[] _reserve;
      _reserve=new char[(strlen(_data)+2)];     //make spcace for buffer
      if(!_reserve) return *this;               //memory out
      strcpy(_reserve, _data);                  //copy data to buffer
   	delete[] _data;
      _data=new char[len+_del_inc+1];                    //make enought space for catenate
      _tot_space=len+_del_inc;
      strcpy(_data, _reserve);                  //copy bac
      delete[] _reserve; _reserve=NULL;         //clean buffer;
   }

   if (_data) _data[_nb_char]=x;               //catenate character
   _data[len]='\0';
   _nb_char=len;
   return *this;

}; // new

my_string
operator+(const my_string& x, const my_string& y){
	int len=strlen(x.c_str())+strlen(y.c_str());
   if(len > 10000){
      my_string z;
      return z;
   }
   my_string z(len+1);
   if (z._data) strcpy(z._data, x.c_str());        //copy x
   if (z._data) strcat(z._data, y.c_str());           //catenate y
   z._nb_char=len;
   return z;
};

my_string
operator+(const my_string& x, char* y){
	int len=strlen(x.c_str())+strlen(y);
   if(len > 10000){
      my_string z;
      return z;
   }
   my_string z(len+1);
   if (z._data) strcpy(z._data, x.c_str());        //copy x
   if (z._data) strcat(z._data, y);           //concatenate y
   z._nb_char=len;
   return z;
};

//------------------------------------------------------------------------------
// new since 21_3_97 //does not read the delimiter
//Loro_28_8_98 Also stops if end of line
bool
my_string::read_to_delim(istream & is, char delim = '\n'){
	(*this)="";
	char curr_char;
   is.get(curr_char);
   //Loro_29_8_98
   for(;curr_char!=EOF && curr_char!=delim && curr_char!='\n' && curr_char!='\0';){
   	(*this)+=curr_char;             //add curr_char to array
      is.get(curr_char);
   }
   if (curr_char==delim) return 1;
   return 0;
};

void
my_string::read_to_delim(istream & is){
	char delim = '\n';
	(*this)="";
	char curr_char;
   is.get(curr_char);
   //Loro_29_8_98
   for(;curr_char!=EOF && curr_char!=delim && curr_char!='\0';){
		(*this)+=curr_char;             //add curr_char to array
      is.get(curr_char);
   }
};

my_string&
my_string::remove( int pos, int n ){  //premiere position =0 !!!

	if(n<=0) return *this;
   if(pos<=0) return *this;

	if(n+pos >= _nb_char){  //all is removed
   	_data[pos]='\0';
   }
   else{
      for(int i=pos; i<_nb_char-n; i++){
      	_data[i]=_data[i+n];
      }
      _data[_nb_char-n]='\0';
   }
   _nb_char=strlen(_data);
   return *this;
};

my_string
my_string::extract_sub_str(char delim){ //extracts until the next occurence of delim or end of string //new 11_11_97

	my_string sub_str;
   int nb_removed=0;
   if(!_nb_char) sub_str;   //return empty string.
   char* buffer=new char[_nb_char+1];   //stef_3_6_98 one char was missing
   //first remove some delimiters at the beginning
   int i=0;
   for(;i<_nb_char && _data[i]==delim; ){
   	i++;
      nb_removed++;
   }
   //then copy the next series of charaters into the buffer until delim occures.
   int pos=0;
   for(;i<_nb_char && _data[i]!=delim; ){
      buffer[pos]=_data[i];
      i++;
      nb_removed++;
      pos++;
   }
   buffer[pos]='\0';
   if(!pos){
      //now remove the first chars.
   	for(i=0;i<_nb_char-nb_removed; i++){
   		_data[i]=_data[i+nb_removed];
   	}
   	_data[_nb_char-nb_removed]='\0';
   	_nb_char-= nb_removed;
		delete[] buffer;
   	return sub_str;   //return empty string.
   }
   sub_str=buffer;   //copy buffer into sub_string;
   //now remove the first chars.
   for(i=0;i<_nb_char-nb_removed; i++){
   	_data[i]=_data[i+nb_removed];
   }
   _data[_nb_char-nb_removed]='\0';
   _nb_char-= nb_removed;
	delete[] buffer;
   return sub_str;
};

int
my_string::find_first_of( const char& s, int pos ) const{

   if(pos > _nb_char) return NPOS;
	for(int i=0;i<_nb_char;i++){
   	if(_data[i]==s) return i;
   }
   return NPOS;
};

int
my_string::find_first_of( const char& s) const{

     //for debug
      //AOS << " ** Entering [find_first_of] with this=" <<  _data << " and s=" << s << endl;
      //end debug

	for(int i=0;i<_nb_char;i++){

      //for debug
      //AOS << " ** comparing s and:" << _data[i] << " at pos:" << i << endl;
      //end debug

      if(_data[i]==s) return i;

   }
   return NPOS;

};

int
my_string::read_line(istream& is){
	(*this)="";
	char curr_char;
   is.get(curr_char);              //skip blanks and null characters....
   //Loro_29_8_98
   for(;curr_char!=EOF && (curr_char==' ' ||  curr_char=='\t') && curr_char!='\0' && is; )  is.get(curr_char);
   for(;curr_char!=EOF && curr_char!='\n';){
   	  if(is){
   		(*this)+=curr_char;             //stef_22_6_98 if is is bad then stop
         //Loro_29_8_98
         if (curr_char=='\0') break;
   	  }
   	  else{
        	//Loro_29_8_98
   	    return 0;
   	  }
   	  is.get(curr_char);
   }
   return 1;
};

int
my_string::is_null(){
	if(_nb_char) return 0;
   else return 1;
};

int
my_string::contains(const char* pat) const{
   int pat_len=strlen(pat);
   if(pat_len==0 || pat_len> _nb_char) return 0;  //pat to long or empty
   for(int i=0; i<_nb_char-pat_len+1; i++){
   	int equal=1;
      for(int j=0;j<pat_len && equal;j++){
      	equal=(_data[j+i]==pat[j]);
      }
      if(equal) return 1;         //if here equal is true we found pat and we return true
   }
	return 0; //not found then return 0;
};

int
my_string::contains(char pat) const{

   for(int i=0; i<_nb_char; i++){
      if(_data[i]==pat) return 1;
   }
	return 0; //not found then return 0;
};

int
my_string::contains(const my_string& s) const{

  return contains(s.c_str());
};

int
my_string::rfind( const my_string& s ){ //return the position of pattern s if not found return -1
	int pat_len=strlen(s.c_str());
   if(pat_len==0 || pat_len> _nb_char) return NPOS;  //pat to long or empty //stef_27_3
   for(int i=0; i<_nb_char-pat_len+1; i++){
   	int equal=1;
      for(int j=0;j<pat_len && equal;j++){
      	equal=(_data[j+i]==(s.c_str())[j]);
      }
      if(equal) return i;         //if here equal is true we found pat and we return position
   }
	return NPOS; //not found then return -1;


};


//------------------------------------------------------------------------------
ostream&
operator<<(ostream& os, const my_string& x) {
	os	<< x.c_str();
   return os;
}

istream&
operator>>(istream& is, my_string& s){
	char curr_char;
   s="";
   is.get(curr_char);
   for(;curr_char!=EOF && curr_char!='\n' && (curr_char==' ' || curr_char=='\t');){ //skipp blanks before my_string
   	is.get(curr_char);
   }

   for(;curr_char!=EOF && curr_char!='\n' && curr_char!=' ';){
   	s+=curr_char;             //add curr_char to array
      is.get(curr_char);
   }
   //if(curr_char=='\n') is.putback(curr_char); //it seams that this operator leaves '\n' in stream???   //stef_10_6_98
   return is;
};

//------------------------------------------------------------------------------
void
my_string::rm_path(){ //stef_28_3_98 removes path info from a file name

	int pos=0;
	for(int i=_nb_char-1;i>=0;i--){
   	if(_data[i]=='/' || _data[i]=='\\'){
      	pos=i;
         i=-1;
      }
   }
   pos++; //goto next character after '/'
   int j=0;   //where we write
   for(;pos<_nb_char;){
   	_data[j]=_data[pos];
      j++;pos++;
   }
   _data[j]='\0';  //terminate
   _nb_char=strlen(_data);
    return;

};
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void my_strrev(char * str) {
   int size=strlen(str);
   char *buf = new char[size+1];
   for (int i=0; i<size; ++i) {
      buf[i]=str[size-1-i];
   }
   for (int i=0; i<size; ++i) {
      str[i]=buf[i];
   }
   delete[] buf;
}
//------------------------------------------------------------------------------
my_string remove_extension(const my_string& s) {
  	my_string sout;
  	char sin[400];
   char *a=strrchr(s.c_str(),'.');
   if(!a) return s;
   char *b=s.c_str();
   int i=0;
   for( ; a!=b; i++, b++)
      sin[i]=*b;
   sin[i]='\0';
   sout=sin;
   return sout;
}
//------------------------------------------------------------------------------
my_string extract_path(const my_string& s) {
	my_string sout;
   char sin1[400];
	char *a=strrchr(s.c_str(),'\\');
	if(!a) return ""; //no backslash
	char *b=s.c_str();
	int i=0;
	for( ; b!=a; i++, b++)
      sin1[i]=*b;
   sin1[i]='\0';
   sout=sin1;
   return sout;
}
//------------------------------------------------------------------------------

