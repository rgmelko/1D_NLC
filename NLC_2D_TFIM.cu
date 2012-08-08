/*******************************************************************
Linked Cluster Expansion program for the 1D TFIM model
Roger Melko, Ann Kallin, Katie Hyatt June 2012
based on a Lanczos code from November 2007
********************************************************************/

using namespace std;

#include <fstream>
#include <limits.h>
#include <cstdio>
#include <algorithm>
#include <time.h>
#include "CPU/Lanczos_07.h"
#include "CPU/GenHam.h"
#include "CPU/simparam.h"
#include "GPU/lanczos.h"
#include "../Graphs/graphs.h"

#define HOW_MANY_16 30 
#define HOW_MANY_18 2 

bool order_16( Graph g ){ return g.Order == 16; };
bool order_18( Graph g ){ return g.Order == 18; };

int main(int argc, char **argv) 
{

    int CurrentArg = 1;
    bool gpuFlag = false;
    string InputFile;
    string OutputFile = "Output_2D.dat";
    while ( CurrentArg < argc )
    {
        if ( argv[ CurrentArg ] == string("-g") || argv[ CurrentArg ] == string("--gpu") )
        {
            gpuFlag = true;
        }
        if ( argv[ CurrentArg ] == string("-i") || argv[ CurrentArg ] == string("--input") )
        {
            InputFile = string(argv[ CurrentArg + 1 ]);
        }
        if ( argv[ CurrentArg ] == string("-o") || argv[ CurrentArg ] == string("--output") )
        {
            OutputFile = string(argv[ CurrentArg + 1 ]);
        }
        CurrentArg++;
    }

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
    vector< Graph > fileGraphs;
    vector< double > WeightHigh;

    ReadGraphsFromFile( fileGraphs, InputFile);

    //int HowMany = 30;

    ofstream fout( OutputFile.c_str() );
    fout.precision( 10 );
    cout.precision( 10 );

    int FirstCount; 
    FirstCount = (int) std::count_if(fileGraphs.begin(), fileGraphs.end(), order_16);
    int SecondCount; 
    SecondCount = (int) std::count_if(fileGraphs.begin(), fileGraphs.end(), order_18);

    J = 1;
    
    for( double hh = 4; hh <= 5; hh += 1 ) 
    {
        h = hh;

        WeightHigh.push_back( -h ); //Weight for site zero
        double RunningSumHigh = WeightHigh[ 0 ];

        d_hamiltonian* HamilLancz = (d_hamiltonian*) malloc(HOW_MANY_16 * sizeof( d_hamiltonian ) );
        parameters* data = (parameters*) malloc( HOW_MANY_16 * sizeof( parameters ) );
        double** groundstates = (double**) malloc( HOW_MANY_16 * sizeof( double* ) );
        double** eigenvalues = (double**) malloc( HOW_MANY_16 * sizeof( double* ) );
        eigenvalues[ 0 ] = (double*) malloc( HOW_MANY_16 * sizeof( double ) );
        int* NumElem = (int*) malloc( HOW_MANY_16 * sizeof( int ) );
        int** Bonds = (int**) malloc( HOW_MANY_16 * sizeof( int* ) );
        unsigned int GPUqueue[HOW_MANY_16];


        int GPUprocessed = 0;
        int GPUmax = HOW_MANY_16;
        int remaining = FirstCount;

        unsigned int i = 1;
        while ( i < fileGraphs.size() )//&& fileGraphs.at(i).Order < 14) //skip the zeroth graph
        {
            if ( gpuFlag && 
                 GPUprocessed < GPUmax &&
                 remaining >= GPUmax &&
                 (fileGraphs[ i ].Order == 16 || fileGraphs[ i ].Order == 18) )
            {
                GPUprocessed++;
                remaining--;
                GPUqueue[GPUprocessed] = i; //store the locations of the graphs we're going to process in parallel

                //energy = 0;
                Bonds[ GPUprocessed ] = ( int* ) malloc( sizeof( int ) * 3 * fileGraphs[ i ].Order );
                for ( unsigned int k = 0; k < fileGraphs[ i ].Order; k++ )
                {
                    Bonds[ GPUprocessed ][ k ] = k;
                    Bonds[ GPUprocessed ][ k + fileGraphs[ i ].Order ] = fileGraphs[ i ].AdjacencyList[ k ].second;
                    Bonds[ GPUprocessed ][ k + 2 * fileGraphs[ i ].Order ] = fileGraphs[ i ].AdjacencyList[ 2 * k + 1 ].second;
                }
                    
                data[ GPUprocessed ].nsite = fileGraphs[ i ].Order;
                data[ GPUprocessed ].Sz = 0;
                data[ GPUprocessed ].dimension = (fileGraphs[i].AdjacencyList.size() <= fileGraphs[i].Order) ? 1 : 2;
                data[ GPUprocessed ].J1 = 4*J;
                data[ GPUprocessed ].J2 = h;
                data[ GPUprocessed ].modelType = 2;
                i++;
                /*cudaEvent_t start, stop;
                cudaEventCreate(&start);
                cudaEventCreate(&stop);
                float time;
                cudaEventRecord(start, 0);
                */
            }
            
            if( GPUprocessed == GPUmax - 1 || remaining == 0)
            {
                ConstructSparseMatrix( GPUprocessed, Bonds, HamilLancz, data, NumElem, 0);
                lanczos( GPUprocessed, NumElem, HamilLancz, groundstates, eigenvalues, 200, 1, 1e-12);
                /*cudaEventRecord(stop, 0);
                cudaEventSynchronize(stop);
                cudaEventElapsedTime(&time, start, stop);
                cout<<"Time to do GPU work: "<<time<<endl;
                cudaEventDestroy(start);
                cudaEventDestroy(stop);
                */
                for( int j = 0; j < GPUprocessed; j++ )
                {
                    energy = eigenvalues[ j ][ 0 ];
                    WeightHigh.push_back( energy );
                    for ( unsigned int k = 0; k < fileGraphs[ GPUqueue[j] ].SubgraphList.size(); k++ )
                        WeightHigh.back() -= fileGraphs[ GPUqueue[j] ].SubgraphList[ k ].second * WeightHigh[ fileGraphs[ GPUqueue[i] ].SubgraphList[ k ].first ];

                    RunningSumHigh += fileGraphs[ GPUqueue[j] ].LatticeConstant * WeightHigh.back();
                }
                GPUprocessed = 0;
                //free(Bonds[0]);
                //cudaFree(HamilLancz[0].rows);
                //cudaFree(HamilLancz[0].cols);
                //cudaFree(HamilLancz[0].vals);
                //cudaFree(eigenvalues[0]);
                //cudaFree(groundstates[0]);
            }
            if ( remaining == 0 )
            {
                remaining = SecondCount;
                GPUmax = HOW_MANY_18;
            }
            
            else
            {
                GENHAM HV( fileGraphs[ i ].Order, J, h, fileGraphs[ i ].AdjacencyList, fileGraphs[ i ].LowField );

                LANCZOS lancz( HV.Vdim );  //dimension of reduced Hilbert space (Sz sector)
                HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos
            
                energy = lancz.Diag( HV, 1, 1, eVec ); // Hamiltonian, # of eigenvalues to converge, 1 for -values only, 2 for vals AND vectors
            
                WeightHigh.push_back( energy );
                for ( unsigned int j = 0; j < fileGraphs[ i ].SubgraphList.size(); j++ )
                    WeightHigh.back() -= fileGraphs[ i ].SubgraphList[ j ].second * WeightHigh[ fileGraphs[ i ].SubgraphList[ j ].first ];

                RunningSumHigh += fileGraphs[ i ].LatticeConstant * WeightHigh.back();
                i++;
            }
        } 
        fout<<"h= "<<h<<" J= "<<J;
        fout <<" Energy= "<< RunningSumHigh<<endl;

        WeightHigh.clear();
        RunningSumHigh=0;
    
    }
    
    fout.close();
    return 0;

}
