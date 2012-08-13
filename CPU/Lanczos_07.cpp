# include "Lanczos_07.h"
#define EVECTS 0  

LANCZOS::LANCZOS(const int Dim_) : Dim (Dim_)
{
  //Dim = Dim_;

  STARTIT = 5;
  CONV_PREC = 1E-10;

  //  Psi.resize(Dim_);
  V0.resize(Dim_); 
  //Vorig.resize(Dim_);
  V1.resize(Dim_);  
  V2.resize(Dim_);

}//constructor


double LANCZOS::Diag(const GENHAM& SparseH, const int Neigen, const int Evects2, vector<long double>& Psi)
// Reduces the Hamiltonian Matrix to a tri-diagonal form 
{

    Psi.resize(Dim);
    int ii;
    int iter, MAXiter, EViter;
    int min;
    int Lexit;
    long double Norm;
    long double E0;
    long double Energy;
    vector<long double> Ord(Neigen); //a vector for the ordered eigenvalues
  
    int LIT;   //max number of Lanczos iterations
    LIT = 100;
 
    //Matrices
    vector<long double> alpha(LIT, 0.);
    vector<long double> beta(LIT, 0.);
    //For ED of tri-di Matrix routine (C)
    int nn, rtn;
    vector<long double> e(LIT, 0.);
    vector<long double> d(LIT, 0.); 
    vector< vector< long double> > Hmatrix;
    Hmatrix.resize(LIT);
    for( int ii = 0; ii < LIT; ii++)
    {
        Hmatrix[ii].resize(LIT);
    }

    //tensor indices
    //firstIndex i;    
    //secondIndex j;
    //thirdIndex k;   
  
    iter = 0;

    for (EViter = 0; EViter < Evects2; EViter++) {//0=get E0 converge, 1=get eigenvec
        iter = 0;
    //create a "random" starting vector
        V0.assign(Dim, 0.0);
        //V0[ 0 ] = 1.0;
        if( V0.size() == 4) V0[1] = 1.0;
        for (unsigned int vi=0; vi < V0.size(); vi++) 
        { 
            //V0[vi] = 1.0;
            if (vi == V0.size() - 1 && V0.size() != 4) V0[vi] = 1.0;
            //else if (vi%5 == 0) V0[vi] = -2.0;
            //else if (vi%7 == 0) V0[vi] = 3.0;
            //else if (vi%9 == 0) V0[vi] = -4.0;
        }
        //if(V0.size()>4){V0[V0.size()-1]=1.0;}
        Normalize(V0);  
   
        if (EViter == 1) for( ii = 0; ii < Dim; ii++) Psi[ii] = V0[ii] * (Hmatrix[0][min]);
    
        V1.assign(Dim, 0.);
        beta[0] = 0;  //beta_0 not defined
    
        //****** do V1=H|V0> below
        apply(V1, SparseH, V0);

        alpha[0] = 0;
        for (ii = 0; ii < Dim; ii++) alpha[0] += V0[ii] * V1[ii];
        for (ii = 0; ii < Dim; ii++) V1[ii] -= alpha[0] * V0[ii];
        Norm = 0;
        for ( ii=0; ii < Dim; ii++ ) Norm += V1[ii]*V1[ii];
        beta[1] = sqrt(Norm);

        for ( ii=0; ii < Dim; ii++ ) V1[ii] /= beta[1];

        if (EViter == 1)
        {
            for( ii = 0; ii < Dim; ii++) Psi[ii] += V1[ii] * (Hmatrix[1][min]);
        }
        // done 0th iteration
        Lexit = 0;   //exit flag
        E0 = 1.0;    //previous iteration GS eigenvalue
        while(Lexit != 1){

            iter++;

            //****** do V2=H|V1> below
            apply(V2, SparseH, V1);

            alpha[iter] = 0.;
            for (ii = 0; ii < Dim; ii++) alpha[iter] += V1[ii] * V2[ii];
            for (ii = 0; ii < Dim; ii++) V2[ii] -= (alpha[iter] * V1[ii] + beta[iter] * V0[ii]);
            Norm = 0;
            for (ii = 0; ii < Dim; ii++) Norm += V2[ii]*V2[ii];
            beta[iter+1] = sqrt(Norm);

            for (ii = 0; ii < Dim; ii++) V2[ii] /= beta[iter+1];

            if (EViter == 1) {
                for(ii = 0; ii < Dim; ii++) Psi[ii] += V2[ii] * (Hmatrix[iter + 1][min]);
            }

            V0 = V1;
            V1 = V2;

            if (iter > STARTIT && EViter == 0){

                //diagonalize tri-di matrix
                d[0] = alpha[0];
                for (ii=1; ii <= iter; ii++){
                    d[ii] = alpha[ii];
                    e[ii-1] = beta[ii];
                }
                e[iter] = 0;

                nn = iter+1;
                rtn = tqli2(d,e,nn,Hmatrix,0);

                //---determin vector (value) of minimal eigenvectors
                Ord.assign(Neigen,999.0); //initialize
                min = 0;
                for (unsigned int oo=0; oo<Ord.size(); oo++) 
                {
                    if ( d[0] < Ord[oo] ) {
                        Ord.insert(Ord.begin() + oo, d[0]);
                        Ord.pop_back();
                        break;
                    }//if
                }
                for (ii=1;ii<=iter;ii++){
                    if (d[ii] < d[min])  min = ii;
                    for (unsigned int o=0; o<Ord.size(); o++) {
                        if ( d[ii] < Ord[o] ) {
                            Ord.insert(Ord.begin() + o, d[ii]);
                            Ord.pop_back();
                            break;
                        }//if
                    }//o
                }//ii

            //********************* UNCOMMENT THIS TO PRINT YOUR ENERGIES
            //for (int o=0; o<Ord.size(); o++)  cout<<setprecision(12)<<Ord.at(o)<<"  ";     
            //cout<<endl;

            if ( (E0 - Ord.back() ) < CONV_PREC && iter > 12) {
                Lexit = 1;
                //cout<<"Lanc :"<<iter<<" ";
                //cout<<setprecision(12)<<d(min)<<"\n";     
                //for (int o=0; o<Ord.size(); o++)  cout<<setprecision(12)<<Ord.at(o)<<"  ";     
                //cout<<endl;
            }
            else {
                E0 = Ord.back(); //E0 = d(min);
            }

            if (iter == LIT-2) {
                LIT += 100;
                //cout<<LIT<<" Resize Lan. it \n";
                d.resize(LIT);
                e.resize(LIT);
                Hmatrix.resize(LIT);
                for(int jj = 0; jj < LIT; jj++) Hmatrix[jj].resize(LIT);
                alpha.resize(LIT);
                beta.resize(LIT);
            }//end resize


        }//end STARTIT

        if (EViter == 1 && iter == MAXiter) Lexit = 1;

    }//while
   
    if (EViter == 0){
        MAXiter = iter;
        //diagonalize tri-di matrix
        d[0] = alpha[0];
        for (ii=1;ii<=iter;ii++){
            d[ii] = alpha[ii];
            e[ii-1] = beta[ii];
        }
        e[iter] = 0;
        //calculate eigenvector
        for ( unsigned int kk = 0; kk < Hmatrix.size(); kk++) Hmatrix[kk].assign(LIT, 0.l);
        for ( ii = 0; ii <= iter; ii++) Hmatrix[ii][ii] = 1.l; //identity matrix
        nn = iter+1;
        rtn = tqli2(d,e,nn,Hmatrix,1);
        min = 0;
        for (ii=1; ii<=iter; ii++)
        {
            if (d[ii] < d[min])  min = ii;
        }
    }
    
    if (EViter == 0) Energy = E0;


  }//repeat (EViter) to transfrom eigenvalues H basis
  
  //Normalize(Psi,N);
  Normalize(Psi);

  //  cout<<Psi<<" Psi \n";

//  V2 = sum(Ham(i,j)*Psi(j),j);
//  for (ii=0;ii<Dim;ii++)
//    cout<<ii<<" "<<V2(ii)/Psi(ii)<<" EVdiv \n";

//  return E0;
  return Energy;

}//end Diag

/*********************************************************************/
void LANCZOS::apply(vector< long double > & U, const GENHAM& H, const vector< long double > & V)  //apply H to |V>
{
    int kk;
    U.assign(U.size(), 0.);

    for (int ii = 0; ii < Dim; ii++) 
    {
        for (int jj = 1; jj <= H.PosHam[ii][0]; jj++)
        {
            kk = H.PosHam[ii][jj]; //position index
            U[ii] += H.ValHam[ii][jj] * V[kk];
            if (ii != kk) U[kk] += H.ValHam[ii][jj] * V[ii]; //contribution to lower half
        }
    }

}//apply


/*********************************************************************/
void LANCZOS::Normalize(vector<long double>& V) //Normalize the input vector (length N)
{
  long double norm;
  unsigned int i;

  norm = 0.0;             
  for (i=0; i<V.size(); i++) norm += V[i]*V[i]; //  <V|V>
  norm = sqrt(norm);

  for (i=0; i<V.size(); i++) V[i] /= norm;
}

/*********************************************************************/
#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))
int LANCZOS::tqli2(vector<long double >& d, vector<long double>& e, int n, vector< vector< long double > >& z, const int Evects)
/***
 April 2005, Roger Melko, modified from Numerical Recipies in C v.2
 modified from www.df.unipi.it/~moruzzi/
 Diagonalizes a tridiagonal matrix: d[] is input as the diagonal elements,
 e[] as the off-diagonal.  If the eigenvalues of the tridiagonal matrix
 are wanted, input z as the identity matrix.  If the eigenvalues of the
 original matrix reduced by tred2 are desired, input z as the matrix
 output by tred2.  The kth column of z returns the normalized eigenvectors,
 corresponding to the eigenvalues output in d[k].
 Feb 23 2005: modified to use Blitz++ arrays
***/
{
    int m,l,iter,i,k;
    long double s,r,p,g,f,dd,c,b;

    for ( l=0; l < n; l++) 
    {
        iter = 0;
        do 
        { 
            for (m = l; m < n-1; m++) 
            { 
                dd = fabs( d[m] ) + fabs( d[m+1] );
                if (fabs( e[m] ) + dd == dd) break;
            }
            if (m != l) { 
                if (iter++ == 30) { 
                    cout <<"Too many iterations in tqli() \n";
                    return 0;
                }
                g = ( d[l+1] - d[l] )/( 2.0 * e[l] );
                r = sqrt( (g * g) + 1.0 );
                g = d[m] - d[l] + e[l] / ( g + SIGN(r,g) );
                s = c = 1.0;
                p = 0.0;
                for ( i = m - 1; i >= l; i--) 
                { 
                    f = s * e[i];
                    b = c * e[i];
                    if ( fabs( f ) >= fabs( g ) ) 
                    { 
                        c = g / f;
                        r = sqrt( (c*c) + 1.0 );
                        e[i + 1] = f * r;
                        c *= ( s = 1.0 / r );
                    }
                    else 
                    { 
                        s = f / g;
                        r = sqrt( (s*s) + 1.0 );
                        e[i + 1] = g * r;
                        s *= ( c = 1.0 / r );
                    }
                    g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    p = s * r;
                    d[i + 1] = g + p;
                    g = c * r - b;
                    /*EVECTS*/
                    if (Evects == 1) 
                    {
                        for ( k = 0; k < n; k++ ) 
                        { 
                            f = z[k][i + 1];
                            z[k][i + 1] = s * z[k][i] + c * f;
                            z[k][i] = c * z[k][i] - s * f;
                        }      
                    }//Evects
                }
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        } while (m != l);
    }
    return 1;
}


/*****************************************************************************/
///
/// @brief Householder reduces a real symmetric matrix a to tridiagonal form
///
/// @param a is the matrix
///
/// On output, a is replaced by the orthogonal matrix effecting the transformation.
/// The diagonal elements are stored in d[], and the offdiagonal elements are
/// stored in e[].
/// March 30 2006: modified to use Blitz++ arrays
///
void LANCZOS::tred3(vector< vector<double> > & a, vector<double> & d, vector<double>& e, const int n)
{
    int i,j,k,l;
    double f,g,h,hh,scale;
 
    for ( i = n-1; i > 0; i--)
    { 
        cout<<i<<endl;
        l = i - 1;
        h = scale = 0.0;
        if ( l > 0 )
        { 
            for ( k = 0; k <= l; k++) scale += fabs( a[i][k] );
            if ( scale == 0.0) 
            {
                e[i] = a[i][l];
                continue;
            }// .....skip transformation
          // .............................. use scaled a's for transformation
            for ( k = 0; k <= l; k++) 
            { 
                a[i][k] /= scale;
                h += a[i][k]*a[i][k];
            }
            f = a[i][l];
            g = sqrt(h);
            if ( f > 0.0) g = -g;
            e[i] = scale * g;
            h -= f * g;
            a[i][l] = f - g;
            f = 0.0;
            for ( j = 0; j <= l; j++ )
            { 
                a[j][i] = a[i][j] / h;
                g = 0.0;
                for ( k = 0; k <= j; k++ ) g += a[j][k] * a[i][k];
                for ( k = j + 1; k <= l; k++ ) g += a[k][j] * a[i][k];
                e[j] = g / h;
                f += e[j]*a[i][j];
            }
            hh = f / (h+h);
            for ( j = 0; j <= l; j++ )
            { 
                f = a[i][j];
                e[j] = g = e[j] - hh * f;
                for ( k = 0; k <= j; k++ ) a[j][k] -= ( f * e[k] + g * a[i][k]);
            }
        }
        else e[i] = a[i][l];
        d[i] = h;
    }
                          // eigenvectors
#if (EVECTS == 1)
    d[0] = 0.0;
    e[0] = 0.0;
    for ( i = 0; i < n; i++ )
    { 
        l = i - 1;
        if ( d[i] )
        { 
            for ( j = 0; j <= l; j++ )
            { 
                g=0.0;
                for ( k = 0; k <= l; k++ ) g += a[i][k] * a[k][j];
                for ( k = 0; k <= l; k++) a[k][j] -= g * a[k][i];
            }
        }
        d[i] = a[i][i];
        a[i][i] = 1.0;
        for ( j = 0; j <= l; j++ ) a[j][i] = a[i][j] = 0.0;
    }
#else 
    for ( i = 0; i < n; i++ ) d[i] = a[i][i];
#endif

    for ( i = 0; i < n - 1; i++ ) e[i] = e[i + 1];
    e[n - 1] = 0;
  // .......................................... complete transformation matrix
  // this line has to be commented!!
  //for (i=0;i<n;i++) for (j=i+1;j<n;j++) a(i,j)=a(j,i);
}
