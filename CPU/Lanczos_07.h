//Lanczos_07.h
//c++ class for performing a Lanczos diagonalization
//Roger Melko, November 2007

#ifndef LANCZOS_07
#define LANCZOS_07


#include <iostream>
#include <limits.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <iomanip>
#include <vector>
using namespace std;

#include "GenHam.h"

class LANCZOS{

  public:
    //Data
    int Dim; //dimension
    //    Array<long double,1> Psi;  //eigenvector

   //Methods
   LANCZOS(const int);
   double Diag(const GENHAM &, const int, const int, vector< long double > &);
   void tred3(vector< vector<double> >& , vector<double>& , vector<double>& e, const int );

  private:
   int STARTIT;
   long double CONV_PREC; //convergence precision

   vector<long double> V0;  
   vector<long double> V1;    //Ground state vector
   vector<long double> V2;

   void apply( vector<long double> &, const GENHAM &, const vector<long double>&);  //apply H to |V>
   void Normalize(vector<long double>& );
   int tqli2(vector<long double>& , vector<long double>& , int , vector< vector<long double > > & , const int );

}; //LANCZOS class


#endif
