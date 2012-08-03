#include "magnetization.h"

double Magnetization( const vector< long double > & Eigenvector, const int NumberSpins)
{
    double chi = 0;
    double norm = 0;
    for ( unsigned int CurrentKet = 0; CurrentKet < Eigenvector.size(); CurrentKet++ )
    {
        int UpSpins = 0;
        for ( int CurrentSpin = 0; CurrentSpin < NumberSpins; CurrentSpin++ )
        {
            UpSpins += ( CurrentKet >> CurrentSpin) & 1;
        }
        double temp = ((2*UpSpins) - NumberSpins )* Eigenvector[ CurrentKet ] * Eigenvector[ CurrentKet ] * ((2*UpSpins) - NumberSpins );
        chi += temp;
        norm += Eigenvector[ CurrentKet ] * Eigenvector[ CurrentKet ];
        //cout<<"Current Ket: "<<CurrentKet<<" value in GS: "<<Eigenvector[ CurrentKet ]<<endl;
    }
    //chi = chi/(double)NumberSpins/norm;
    chi /= norm;
    chi /= NumberSpins;
    return abs(chi);
}
