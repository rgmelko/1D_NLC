/*******************************************************************
Linked Cluster Expansion program for the 1D TFIM model
Roger Melko, Ann Kallin, Katie Hyatt June 2012
based on a Lanczos code from November 2007
********************************************************************/

using namespace std;

#include <fstream>
#include <limits.h>
#include <cstdio>
#include <time.h>
#include "CPU/Lanczos_07.h"
#include "CPU/GenHam.h"
#include "CPU/simparam.h"
#include "../CUDA/Lanczos/lanczos.h"
#include "graphs.h"

int main(int argc, char **argv) 
{

    double energy;

    PARAMS prm;  //Read parameters from param.dat  : see simparam.h
    double J;
    double h;

    vector< long double > eVec;

    J=prm.JJ_;
    h=prm.hh_;

    //ifstream fin("rectanglegraphs.dat");
    //bool TypeFlag = false;
    /*fin >> TypeFlag;
    fin.close();*/
    vector< graph > fileGraphs;
    vector< double > WeightHigh;

    ReadGraphsFromFile(fileGraphs, "rectanglegraphs.dat");

    //int HowMany = 30;

    ofstream fout("output_2D.dat");
    fout.precision(10);
    cout.precision(10);

    J=1;
    
    for(int hh=1; hh<10; hh++) 
    {
        h = hh;

        WeightHigh.push_back(-h); //Weight for site zero
        double RunningSumHigh = WeightHigh[0];

        d_hamiltonian* HamilLancz = (d_hamiltonian*)malloc(sizeof(d_hamiltonian));
        parameters* data = (parameters*)malloc(sizeof(parameters));
        double** groundstates = (double**)malloc(sizeof(double*));
        double** eigenvalues = (double**)malloc(sizeof(double*));
        eigenvalues[ 0 ] = (double*)malloc(sizeof(double));
        int* NumElem = (int*)malloc(sizeof(int));
        int** Bonds = (int**)malloc(sizeof(int*));
        
        unsigned int i = 1;

        while ( i<fileGraphs.size() )//&& fileGraphs.at(i).NumberSites < 14) //skip the zeroth graph
        {
            cout<<fileGraphs[i].NumberSites<<endl; 
            if ( (argv[0] == "--gpu" || argv[0] == "-g") && fileGraphs[i].NumberSites > 14)
            {

                Bonds[ 0 ] = (int*)malloc(sizeof(int)*3*fileGraphs[i].NumberSites);
                for(unsigned int k = 0; k < fileGraphs[i].NumberSites; k++)
                {
                    Bonds[ 0 ][ k ] = k;
                    Bonds[ 0 ][ k + fileGraphs[i].NumberSites ] = fileGraphs[i].AdjacencyList.at(2*k).second;
                    Bonds[ 0 ][ k + 2*fileGraphs[i].NumberSites ] = fileGraphs[i].AdjacencyList.at(2*k + 1).second;
                }
                    
                data[ 0 ].Sz = 0;
                data[ 0 ].dimension = 2;
                data[ 0 ].J1 = J;
                data[ 0 ].J2 = h;
                data[ 0 ].modelType = 2;
                ConstructSparseMatrix( 1, Bonds, HamilLancz, data, NumElem, 1);
                lanczos( 1, NumElem, HamilLancz, groundstates, eigenvalues, 200, 1, 1e-12);
                
                energy = eigenvalues[ 0 ][0];
            }
            
            else
            {
                GENHAM HV(fileGraphs[i].NumberSites, J, h, fileGraphs[i].AdjacencyList, fileGraphs[i].LowField);

                LANCZOS lancz(HV.Vdim);  //dimension of reduced Hilbert space (Sz sector)
                HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos
            
                energy = lancz.Diag(HV, 1, 1, eVec); // Hamiltonian, # of eigenvalues to converge, 1 for -values only, 2 for vals AND vectors
            }

            WeightHigh.push_back(energy);
            for (unsigned int j = 0; j < fileGraphs[i].SubgraphList.size(); j++)
                WeightHigh.back() -= fileGraphs[i].SubgraphList[j].second * WeightHigh[fileGraphs[i].SubgraphList[j].first];

            cout<<"h="<<h<<" J="<<J<<" graph #"<<i<<"  ";
            //        cout<<" energy "<<setprecision(12)<<energy<<endl;
            //        cout<<"WeightHigh["<<i<<"] = "<<WeightHigh.back()<<endl;
            RunningSumHigh += WeightHigh.back();
            cout<<"RunningSumHigh = "<<RunningSumHigh;
            cout<<endl;
            i++;
        } 
        /*
        if( argv[0] == "--gpu" || argv[0] == "-g" )
        {
            while ( i < fileGraphs.size() )
            {
                i += 30;
                if (fileGraphs.at(i).NumberSites == 18 )
                {
                    HowMany = 2;
                }
                for( int j = 0; j < HowMany; j++)
                {
                
                    Bonds[ j ] = (int*)malloc(sizeof(int)*3*fileGraphs.at(i - j).NumberSites);
                    for(unsigned int k = 0; k < fileGraphs.at(i - j).NumberSites; k++)
                    {
                        Bonds[ j ][ k ] = k;
                        Bonds[ j ][ k + fileGraphs.at(i - j).NumberSites ] = fileGraphs.at(i - j).AdjacencyList.at(2*k).second;
                        Bonds[ j ][ k + 2*fileGraphs.at(i - j).NumberSites ] = fileGraphs.at(i - j).AdjacencyList.at(2*k + 1).second;
                    }
                    
                    data[ j ].Sz = 0;
                    data[ j ].dimension = 2;
                    data[ j ].J1 = J;
                    data[ j ].J2 = h;
                    data[ j ].modelType = 2;
                    eigenvalues[ j ] = (double*)malloc(sizeof(double));
                }
                
                ConstructSparseMatrix(HowMany, Bonds, HamilLancz, data, NumElem, 1);
                lanczos(HowMany, NumElem, HamilLancz, groundstates, eigenvalues, 200, 1, 1e-12);
                
                for( int j = 0; j < HowMany; j++)
                {
                    energy = eigenvalues[ HowMany - 1 - j ][0];
                    WeightHigh.push_back(energy);
                    for( unsigned int k = 0; k < fileGraphs.at(i - j).SubgraphList.size(); k++)
                    {
                        WeightHigh.back() -= fileGraphs.at(i - j).SubgraphList[k].second * WeightHigh[fileGraphs.at(i - j).SubgraphList[k].first];

                        cout<<"h="<<h<<" J="<<J<<" graph #"<<i - j<<"  ";
                        //cout<<" energy "<<setprecision(12)<<energy<<endl;
                        //cout<<"WeightHigh["<<i<<"] = "<<WeightHigh.back()<<endl;
                        RunningSumHigh += WeightHigh.back();
                        cout <<"RunningSumHigh = "<< RunningSumHigh;
                        cout<<endl;
                    }
                    free(Bonds[j]);
                    cudaFree(groundstates[j]);
                    cudaFree(eigenvalues[j]);
                    cudaFree(HamilLancz[j].rows);
                    cudaFree(HamilLancz[j].cols);
                    cudaFree(HamilLancz[j].vals);
                }
            }
        }    

        else 
        {
            while ( i < fileGraphs.size() )
            {
        //---High-Field---
                GENHAM HV(fileGraphs.at(i).NumberSites, J, h, fileGraphs.at(i).AdjacencyList, fileGraphs.at(i).LowField);

                LANCZOS lancz(HV.Vdim);  //dimension of reduced Hilbert space (Sz sector)
                HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos
                energy = lancz.Diag(HV, 1, prm.valvec_, eVec); // Hamiltonian, # of eigenvalues to converge, 1 for -values only, 2 for vals AND vectors
                WeightHigh.push_back(energy);
                for (unsigned int j = 0; j<fileGraphs.at(i).SubgraphList.size(); j++)
                    WeightHigh.back() -= fileGraphs.at(i).SubgraphList[j].second * WeightHigh[fileGraphs.at(i).SubgraphList[j].first];

                cout<<"h="<<h<<" J="<<J<<" graph #"<<i<<"  ";
                //cout<<" energy "<<setprecision(12)<<energy<<endl;
                //cout<<"WeightHigh["<<i<<"] = "<<WeightHigh.back()<<endl;
                RunningSumHigh += WeightHigh.back();
                cout <<"RunningSumHigh = "<< RunningSumHigh;
                cout<<endl;
                i++;
            }
        }*/

        fout<<"h= "<<h<<" J= "<<J;
        fout <<" Energy= "<< RunningSumHigh<< endl<<endl;

        WeightHigh.clear();
        RunningSumHigh=0;
    
    }
    
    fout.close();
    return 0;

}
