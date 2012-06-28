#include "header.h"
#include "lapack.h"
#include <vector>
#include <fstream>
#include <blitz/array.h>
BZ_USING_NAMESPACE(blitz)
using namespace std;

int main()
{
  int Dim=65536;
  int Vdim=0;
  int Nsite=16;
  vector <int> basis, basisStart(Nsite+2,0);

  int temp;
  for(int ups=0; ups<=16; ups++){
    for (unsigned int i=0; i<Dim; i++){
      temp = 0;
      for (int sp =0; sp<Nsite; sp++){ temp += (i>>sp)&1; }
      if (temp==ups){ basis.push_back(i); basisStart[ups+1]++; Vdim++;}
    }
    basisStart[ups+2]+=basisStart[ups+1];
  }

  cout<<"Dim "<<Dim<<"   Vdim "<<Vdim<<"   BasisStates "<<basisStart[17]<<endl;

  int Bsize = 9;
  int Bdim = 512;
  int Adim = 128;
  double Z = 0;
  Array<double,2> supermat(Adim,Bdim);
  Array<double,2> dmA(Adim,Adim), dmB(Bdim,Bdim);
  vector<double> enteigs;
  double renyi[3], vN[3];
  int a(0),b(0),c(0);
  double temp1, temp2, temp3;
  vector <double> eigs;
  double initial(0),final(0),increment(0);
  ifstream paramin("param.txt");
  paramin >> initial >> final >> increment;
  paramin.close();
  int sectorDim[17]={1,16,120,560,1820,4368,8008,11440,12870,11440,8008,4368,1820,560,120,16,1};

  for(double temperature=initial; temperature<final; temperature+=increment){

    Z=0;
    dmA=0; dmB=0;

    ifstream getZ("//Volumes/sask_raid/region1/separate/all_values.dat");
    for(int i=0;i<65536;i++){getZ>>temp1;Z+=exp(-temp1/temperature);}
    getZ.close();
    cout<<"Done getting Z="<<Z<<endl ;

    ifstream eig_in("//Volumes/sask_raid/region1/separate/all_values.dat");
    ifstream vec_in("//Volumes/sask_raid/region1/separate/all_vectors.dat");
	vN[2]=0; renyi[2]=0;
    for(int sector=0; sector<=16; sector++){
      
      eigs.resize(sectorDim[sector],-999);
      for(int i=0; i<sectorDim[sector]; i++){
	eig_in>>eigs[i];
	eigs[i]=exp(-eigs[i]/temperature)/Z; 
	renyi[2]+=eigs[i]*eigs[i];
	temp2=log(eigs[i]);
	if(!(temp2>-1000000000)){temp2=0;}
	vN[2]+=-eigs[i]*temp2;
      }
      a=0;b=0;c=0;temp1=0;    
      
      for(int w=0; w<sectorDim[sector]; w++){ 
	supermat=0;
	for(int i=0; i<sectorDim[sector]; i++){ 
	  
	  a = basis[i+basisStart[sector]];
	  b = (((a>>11)&31)<<2)+(((a>>7)&1)<<1)+((a>>3)&1);
	  c = (((a>>8)&7)<<6)+(((a>>4)&7)<<3)+(a&7);
	  vec_in>>supermat(b,c); 
	}
	for(int l=0; l<Adim; l++){
	  for(int j=0; j<Adim; j++){
	    temp1=0;
	    for(int k=0; k<Bdim; k++){
	      temp1 += supermat(l,k)*supermat(j,k);
	    }
	    dmA(l,j) += temp1*eigs[w];
	  }
	}	
	for(int l=0; l<Bdim; l++){
	  for(int j=0; j<Bdim; j++){
	    temp1=0;
	    for(int k=0; k<Adim; k++){
	      temp1 += supermat(k,l)*supermat(k,j);
	    }
	    dmB(l,j) += temp1*eigs[w];
	  }
	}
      }
      
      cout << "sector " << sector << " done\n";
	cout.flush();
    }
    eig_in.close();
    vec_in.close();


    // ---- Region B --------  
    enteigs.resize(0);
    diagWithLapack_R(dmA,enteigs);
    renyi[0]=0; vN[0]=0;
    temp2=0; temp3=0;
    
    for(int s=0; s<enteigs.size(); s++){
      temp2+=enteigs[s];
      renyi[0]+=enteigs[s]*enteigs[s];
      temp3=log(enteigs[s]);
      if(!(temp3>-1000000000)){temp3=0;}
      vN[0]+=-enteigs[s]*temp3;
    }
    cout << temperature << "     " << -log(renyi[0]) << endl;
    cout << "eigenvalue sum = " << temp2 << endl;
    

    ofstream renout("Renyi_region3c.dat",ios::app);
    ofstream vnout("vonNeumann_region3c.dat",ios::app);
 

    renout<<temperature <<"      " <<setprecision(12) <<-log(renyi[0])<<endl;
    vnout <<temperature <<"      " <<setprecision(12) <<vN[0] << endl;
    
    renout.close();
    vnout.close();

    // ---- Region A ----------
    enteigs.resize(0);
    diagWithLapack_R(dmB,enteigs);
    renyi[1]=0; vN[1]=0;
    temp2=0; temp3=0;
    
    for(int s=0; s<enteigs.size(); s++){
      temp2+=enteigs[s];
      renyi[1]+=enteigs[s]*enteigs[s];
      temp3=log(enteigs[s]);
      if(!(temp3>-1000000000)){temp3=0;}
      vN[1]+=-enteigs[s]*temp3;
    }

    ofstream renout2("Renyi_region3.dat",ios::app);
    ofstream vnout2("vonNeumann_region3.dat",ios::app);
 

    renout2<<temperature <<"      " <<setprecision(12) <<-log(renyi[1])<<endl;
    vnout2 <<temperature <<"      " <<setprecision(12) <<vN[1] << endl;
    
    renout2.close();
    vnout2.close();

    // ---- mutual info ------
    ofstream renout3("Renyi_region3m.dat",ios::app);
    ofstream vnout3("vonNeumann_region3m.dat",ios::app);
 

    renout3<<temperature <<"      " <<setprecision(12) <<-log(renyi[1])-log(renyi[0])+log(renyi[2])<<endl;
    vnout3 <<temperature <<"      " <<setprecision(12) <<vN[1]+vN[0]-vN[2] << endl;
    
    renout3.close();
    vnout3.close();
  }
  
  return 0;
}
