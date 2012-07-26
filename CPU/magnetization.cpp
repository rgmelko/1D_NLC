#include "magnetization.h"

double Magnetization( vector< long double > & Eigenvector, int NumberSpins)
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

        chi += abs( (UpSpins * 2) - NumberSpins) * Eigenvector[ CurrentKet ] * Eigenvector[ CurrentKet ];
        norm += Eigenvector[ CurrentKet ] * Eigenvector[ CurrentKet ];
    }
    chi /= (double)NumberSpins;
    chi /= norm;
    return chi;
}
