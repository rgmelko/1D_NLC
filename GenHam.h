//c++ class for creating a general Hamiltonian in c++ vectors
//Roger Melko, November 2007
#ifndef GenHam_H
#define GenHam_H

#include <iostream>
#include <vector>
using namespace std;

#include <blitz/array.h>
BZ_USING_NAMESPACE(blitz)

typedef long double h_float;  //precision for Hamiltonian storage

class GENHAM{

  public:
    int Vdim; //dimenson of reduced Hilbert space

    vector<vector<long> > PosHam;
    vector<vector<h_float> > ValHam;
    //vector<double> DiagHam;

    vector<long> Basis;
    vector<long> BasPos;

    Array<int,2> Bond;
    
    Array<double,2> Ham;  //full hamiltonian

    GENHAM(const int,const h_float J_, const h_float h_,const h_float Q_,const int Sz); 
    void printg();
    //double at(const int , const int );
    Array<double,1> apply(const Array<double,1>&);

    void Bonds_16B();
    void SparseHamJQ();

  private:
    int Nsite; //number sites


    h_float JJ; //heisenberg exchange value
    h_float hh; //next-nearest neighbor exchange value
    h_float QQ; //ring-exchange value

    double HdiagPart(const long);
    double HOFFdBondX(const int, const long);
    double HOFFdBondY(const int, const long);

};

#endif
