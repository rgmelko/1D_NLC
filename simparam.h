#ifndef PARAMSIM_H
#define PARAMSIM_H


//Class to read in the simulation parameters from a file
// Adapted to the JQ ED code

class PARAMS
{
  public:
    double NN_; //the number of lattice sites
    double JJ_; //the heisenberg exchange
    double hh_; //the next-nearest neighbor heisenberg exchange

    int valvec_; //  1 for -values only, 2 for vals AND vectors
    // FULL_DIAG?

    PARAMS(){
      //initializes commonly used parameters from a file
      ifstream pfin;
      pfin.open("param.dat");
    
      pfin >> NN_;
      pfin >> JJ_;
      pfin >> hh_;
      pfin >> valvec_;
    
      pfin.close();
    }//constructor

}; //PARAMS

#endif
