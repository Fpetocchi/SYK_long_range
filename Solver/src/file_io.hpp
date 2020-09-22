#ifndef ___FILE_IO___
#define ___FILE_IO___
//
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <random>
#include <complex>
#include <cmath>
#include <math.h>
#include <time.h>
#include <string>
//
#include "string.h"
#include <cassert>
#include <sys/stat.h>
//
using namespace std;
//
#ifndef CPLX
#define CPLX std::complex<double>
#define INPUT_LINELEN 1024
#endif


//============================================================================//


char* delimiter=" =";


//------------------------------------------------------------------------------


void uncomment_line( char *line, int len)
{
   int i;
   int len1=strlen(line);
   //
   for(i=0;i<len1;i++)
   {
       if(line[i]=='#')
       {
          line[i]='\0';
          break;
       }
   }
}


//------------------------------------------------------------------------------


// read first line in file containing the string flag
char* find_line_with_flag( char *file, const char *flag, char *line )
{
   std::ifstream in;
   char *pch=0;
   try
   {
      in.exceptions(std::ifstream::eofbit |std::ifstream::badbit | std::ifstream::failbit);
      in.open(file,std::ios::in);
      do
      {
         in.getline(line,INPUT_LINELEN);
         uncomment_line(line,INPUT_LINELEN);
         if((pch=strstr(line,flag))) break;
      }while(pch==NULL);
      //
      in.close();
   }
   catch(...)
   {
      std::cerr << "some error occurred in file " << file << std::endl;
      std::cerr << "while looking for line with " << flag << std::endl;
      throw;
   }
   return pch;
}


//------------------------------------------------------------------------------


// split line into n words, separated bt characters in delim
// delim-characters in line are replaced by \0
int strtok(char *line, std::vector<char*> &words, const char *delim )
{
   char *pch;
   int n=0;
   words.resize(0);
   pch=strtok(line,delim);
   while(pch!=NULL)
   {
      n++;
      words.resize(n);
      words[n-1]=pch;
      pch=strtok(NULL,delim);
   }
   return n;
}


//------------------------------------------------------------------------------


//double datatype
void read_param( std::vector<char*> &words, int m, double &data)
{
   try
   {
      if((int)words.size()<m+1) throw("read_param: too few words in line");
      if(!(std::stringstream(words[m]) >> data)) throw("read_param_tvector: cannot interpret arg");
   }
   catch(...)
   {
      std::cerr << "error in read_param(double) " << std::endl;
      throw;
   }
}

//int datatype
void read_param( std::vector<char*> &words, int m, int &data )
{
   try
   {
      if((int)words.size()<m+1) throw("read_param: too few words in line");
      if(!(std::stringstream(words[m]) >> data)) throw("read_param_tvector: cannot interpret arg");
   }
   catch(...)
   {
      std::cerr << "error in read_param(int) " << std::endl;
      throw;
   }
}

//bool datatype
void read_param(std::vector<char*> &words, int m, bool &data)
{
   try
   {
     if(!(std::stringstream(words[m]) >> std::boolalpha >> data)) throw("read_param_tvector: cannot interpret arg");
   }
   catch(...)
   {
      std::cerr << "error in read_param(bool) " << std::endl;
      throw;
   }
}

//cmplx datatype
void read_param(std::vector<char*> &words,int m,CPLX &data)
{
   double real, imag;
   try
   {
      if((int)words.size()<m+2) throw("read_param: too few words in line");
      //if(!(std::stringstream(words[m]) >> data.real()))
      if(!(std::stringstream(words[m]) >> real))
         throw("read_param_tvector: cannot interpret arg");
      //if(!(std::stringstream(words[m+1]) >> data.imag()))
      if(!(std::stringstream(words[m+1]) >> imag))
         throw("read_param_tvector: cannot interpret arg");
   }
   catch(...)
   {
      std::cerr << "error in read_param(cmplx) " << std::endl;
      throw;
   }
      data.real(real);
      data.imag(imag);
}


//------------------------------------------------------------------------------


//generic template
template<typename T> void find_param( char *file, const char *flag, T &x )
{
   try
   {
      char line[INPUT_LINELEN],*pch;
      std::vector<char*> words;
      find_line_with_flag(file,flag,line);
      pch=strstr(line,flag);
      strtok(pch+strlen(flag),words,delimiter);
      read_param(words,0,x);
   }
   catch(...)
   {
      std::cerr << "error in find_param(x) template" << std::endl;
      throw;
   }
}


void find_param(char *file, const char *flag, char data[2])
{
   try
   {
      char line[INPUT_LINELEN],*pch;
      std::vector<char*> words;
      find_line_with_flag(file,flag,line);
      pch=strstr(line,flag);
      strtok(pch+strlen(flag),words,delimiter);
      //for(int len=0; len < words.size(); len++) printf ("%s\n",words[len]);
      //printf ("end\n\n");
      strcpy(data,words[0]);
      //read_param(words,0,x);
   }
   catch(...){
      std::cerr << "error in read_param(char)" << std::endl;
      throw;
   }
}


//-------------------------------------------------------------------------


bool PathExist(char* s)
{
  struct stat buffer;
  return (stat (s, &buffer) == 0);
}


//-------------------------------------------------------------------------


void print_line_minus(int width) {
  std::cout << std::string(width, '-') << std::endl;
}

void print_line_plus(int width) {
  std::cout << std::string(width, '+') << std::endl;
}

void print_line_equal(int width) {
  std::cout << std::string(width, '=') << std::endl;
}

void print_line_star(int width) {
  std::cout << std::string(width, '*') << std::endl;
}

void print_line_dot(int width) {
  std::cout << std::string(width, '.') << std::endl;
}




#undef CPLX
#endif
