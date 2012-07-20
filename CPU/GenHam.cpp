#include "GenHam.h"

//----------------------------------------------------------
GENHAM::GENHAM(const int Ns, const long double J_, const long double h_, vector <pair<int,int> > BBond_, bool Field)  
  : JJ(J_), hh(h_), Bond(BBond_)
//create bases and determine dim of full Hilbert space
{
  int Dim;
  Nsite = Ns;
  LowField = Field;

  Dim = 2;  //  S=1/2 models : two states
  for (int ch=1; ch<Nsite; ch++) Dim *=2;

  BasPos.resize(Dim,-1); //initialization 

  Vdim=0;
  unsigned long temp;    //create basis (16 site cluster)

  for (unsigned long i1=0; i1<Dim; i1++) 
  {
      temp = 0;
      for (int sp =0; sp<Nsite; sp++)
          temp += (i1>>sp)&1;  //unpack bra
          Basis.push_back(i1);
          BasPos.at(i1)=Basis.size()-1;
          Vdim++;
  }//Dim

  //cout<<"Vdim "<<Vdim<<" "<<Dim<<endl;

}//constructor


//----------------------------------------------------------
void GENHAM::printg()
{
  int i,j;
  vector<int> tempP;
  vector<long double> tempV;

  for (i=0; i<PosHam.size(); i++){
    //cout<<PosHam[i][0]<<" * ";
    cout<<i+1<<" * ";
    for (j=0; j<=PosHam[i][0]; j++){
      cout<<"("<<PosHam[i][j]+1<<","<<ValHam[i][j]<<") ";
    }
    cout<<endl;
  }

}//print


//----------------------------------------------------------
void GENHAM::SparseHamJQ()
{
  int ii, jj;

  vector<long> tempBas;
  vector<long double> tempH;
  unsigned long tempi, tempj, tempod;
  int si, sj,sk,sl;
  double tempD;

  for (ii=0; ii<Basis.size(); ii++){
    tempH.clear(); 
    tempBas.clear();

    tempi = Basis.at(ii);
    tempBas.push_back(0); //first element (Row size)
    tempH.push_back(0); //make base 0

    //-----1:   diagonal 
    tempBas.push_back(BasPos.at(tempi));  
    tempD = (*this).HdiagPart(tempi,Nsite);
    tempH.push_back(tempD); 

    for (int T0=0; T0<Nsite; T0++){ //T0 is your square index

      si = T0; //si = Bond(T0,0); //the lower left bond spin is not always T0
      //if (si != T0) cout<<"Square error 2\n";
      //-----2:   first bond (Horizontal)
      tempod = tempi;
      // sj = Bond(T0,1); 
      tempod = tempod^(1<<si);   //toggle bit 

      if (BasPos.at(tempod) > ii){ //build only upper half of matrix
        tempBas.push_back(BasPos.at(tempod));
        tempD = (*this).HOFFdBondX(T0,tempi);
        tempH.push_back(tempD); 
      }
  
    }//si

    tempBas.at(0) = tempBas.size()-1;
    //cout<<tempBas.at(0)<<" "<<tempBas.size()<<" "<<tempH.size()<<endl;

    //bubble sort (slow) 
    long stemp;
    bool noswap = false;
    while (noswap == false){
      noswap = true; 
      for (int i2=1; i2<tempBas.size()-1; i2++){ //ignore 0 element
        if (tempBas.at(i2) > tempBas.at(i2+1) ) {
          stemp = tempBas.at(i2);
          tempBas.at(i2) = tempBas.at(i2+1);
          tempBas.at(i2+1) = stemp;
          tempD = tempH.at(i2);
          tempH.at(i2) = tempH.at(i2+1);
          tempH.at(i2+1) = tempD;
          noswap = false;
        }
      }//i2
    }//while

    PosHam.push_back(tempBas);
    ValHam.push_back(tempH);

  }//ii       

}//Heisenberg

//----------------------------------------------------------
double GENHAM::HdiagPart(const long bra, int Sites){

  int S0b,S1b ;  //spins (bra 
  int T0,T1;  //site
  double valH = 0;

  for (int Ti=0; Ti<Bond.size(); Ti++){
    //***HEISENBERG PART

    T0 = Bond[Ti].first; //T0 = Bond(Ti,0); //lower left spin
    S0b = (bra>>T0)&1;  
    //if (T0 != Ti) cout<<"Square error 3\n";
    T1 = Bond[Ti].second; //T1 = Bond(Ti,1); //first bond
    S1b = (bra>>T1)&1;  //unpack bra

    valH += -JJ*2*(S0b-0.5)*2*(S1b-0.5);

  }//T0

  T0=0;
  T1=Sites-1;
  S0b = (bra>>T0)&1;
  S1b = (bra>>T1)&1;
  if(LowField){ valH += -JJ*2*((S0b-0.5) + (S1b-0.5)); }

  //cout<<bra<<" "<<valH<<endl;

  return valH;

}//HdiagPart

//----------------------------------------------------------
double GENHAM::HOFFdBondX(const int si, const long bra){

  double valH;

  valH = -hh; //contribution from the J part of the Hamiltonian

  return valH;

}//HOFFdPart



