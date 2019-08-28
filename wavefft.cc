/*  --------------------------------------------------------------------
 *
 *  Multivariate complex fourier transform, computed in place
 *  using mixed-radix fast fourier transform algorithm.
 *
 *  Translated from Fortran to C++ by
 *  A. Sazonov (sazonov@thsun1.jinr.ru)
 *  March 2000.
 *
 *  Original Fortran code by
 *  R. C. Singleton, Stanford Research Institute, Sept. 1968
 *  see source at http://www.netlib.org/go/wavefft.f
 *  or           http://www.numis.nwu.edu/ftp/pub/transforms/wavefft.f
 *
 *  -------------------------------------------------------------------
 */

#include <math.h>
#ifndef _WAVEFFT_H
  #include "wavefft.hh"
#endif

#include <iostream>
//#ifndef _STREAM_H
  //#include <stream.h>
//#endif

using namespace std;

void wavefft(double a[], double b[], int ntot, int n, int nspan, int isn)
{
/* ----------------------------------------------------------------
 *  Indexing of external arrays a[] and b[] is changed according to C++
 *  requirements, internal arrays ( nfac[], np[], at[], ck[], bt[], sk[] )
 *  are still using index starting from 1.
 *
 *  A.Sazonov (sazonov@thsun1.jinr.ru)
 * -----------------------------------------------------------------
 */

/*
c  multivariate complex fourier transform, computed in place
c    using mixed-radix fast fourier transform algorithm.
c  by r. c. singleton, stanford research institute, sept. 1968
c  arrays a and b originally hold the real and imaginary
c    components of the data, and return the real and
c    imaginary components of the resulting fourier coefficients.
c  multivariate data is indexed according to the fortran
c    array element successor function, without limit
c    on the number of implied multiple subscripts.
c    the subroutine is called once for each variate.
c    the calls for a multivariate transform may be in any order.
c  ntot is the total number of complex data values.
c  n is the dimension of the current variable.
c  nspan/n is the spacing of consecutive data values
c    while indexing the current variable.
c  the sign of isn determines the sign of the complex
c    exponential, and the magnitude of isn is normally one.
c  a tri-variate transform with a(n1,n2,n3), b(n1,n2,n3)
c    is computed by
c      call wavefft(a,b,n1*n2*n3,n1,n1,1)
c      call wavefft(a,b,n1*n2*n3,n2,n1*n2,1)
c      call wavefft(a,b,n1*n2*n3,n3,n1*n2*n3,1)
c  for a single-variate transform,
c    ntot = n = nspan = (number of complex data values), e.g.
c      call wavefft(a,b,n,n,n,1)
c  the data can alternatively be stored in a single complex array c
c    in standard fortran fashion, i.e. alternating real and imaginary
c    parts. then with most fortran compilers, the complex array c can
c    be equivalenced to a real array a, the magnitude of isn changed
c    to two to give correct indexing increment, and a(1) and a(2) used
c    to pass the initial addresses for the sequences of real and
c    imaginary values, e.g.
c       complex c(ntot)
c       real    a(2*ntot)
c       equivalence (c(1),a(1))
c       call wavefft(a(1),a(2),ntot,n,nspan,2)
c  arrays at(maxf), ck(maxf), bt(maxf), sk(maxf), and np(maxp)
c    are used for temporary storage.  if the available storage
c    is insufficient, the program is terminated by a stop.
c    maxf must be .ge. the maximum prime factor of n.
c    maxp must be .gt. the number of prime factors of n.
c    in addition, if the square-free portion k of n has two or
c    more prime factors, then maxp must be .ge. k-1.
c  ***************************************
*/

//---------------------------------------------------------------
// Notes below are added by A.Sazonov, May 2000:
//---------------------------------------------------------------
//  arrays at[maxf], ck[maxf], bt[maxf], sk[maxf], and np[maxp]
//  are used for temporary storage.
//     maxf must be >=   the maximum prime factor of n.
//     maxp must be > the number of prime factors of n.
//  in addition, if the square-free portion k of n has two or
//  more prime factors, then must be maxp >= k-1.
//---------------------------------------------------------------
// The advantage of C++ dynamic memory allocation is used
// to create arrays with sufficient lengths.
// The lengths of arrays at[maxf], ck[maxf], bt[maxf], sk[maxf]
// may be calculated at the begining of program. 
// The array np[maxp] initially created with length maxp+1
// will be resized as soon as it will be found insufficient
// for further calculations.
//---------------------------------------------------------------
  int maxp=209;		// initial value
  int maxf=1;           // maxf will be adjusted later
  int nf=32;
//  array storage in 'nfac' for a maximum of 32 prime factors of n.
//  it's surely enough if n <= 2^32, and minimal factor is 2
  int nfac[nf];

  if (n < 2) return;
  int inc=isn;
  double c72=0.30901699437494742;		// cos(72 deg)
  double s72=0.95105651629515357;		// sin(72 deg)
  double s120=0.86602540378443865;		// sin(120 deg)
  double pi=3.141592653589793;
  double rad=2.*pi;
  double c1, c2=0., c3=0., s1, s2=0., s3=0.;

  if (isn < 0) { s72=-s72; s120=-s120; rad=-rad; inc=-inc; }

  int nt=inc*ntot;
  int ks=inc*nspan;
  int kspnn, kspan=ks;
  int nn=nt-inc;
  int jc=ks/n;
  double radf=rad*double(jc)*0.5;
  int i=0, jf=0;

  for ( int i=1; i<nf; i++)  nfac[i]=0;

//  determine the factors of n
  int kt, m=0, k=n;

  while ( ((k/16)*16) == k )  { m++; nfac[m]=4; k/=16; }

  int j=3, jj=9;

  do {
       while ( (k%jj) == 0) { m++; nfac[m]=j; k/=jj ;}
       j+=2;
       jj=j*j;
     }                         while ( jj <= k ) ;

  if ( k <= 4 ) { kt=m; nfac[m+1]=k; if (k != 1) m++; }
  else
  {
    if ( ((k/4)*4) == k ) { m++; nfac[m]=2; k/=4; }

    kt=m; j=2;

    do {
         if( (k%j) == 0) { m++; nfac[m]=j; k/=j; }
         j=((j+1)/2)*2+1;
       }                        while ( j <= k ) ;
   }

   if ( kt != 0 )
   { j=kt; do { m++; nfac[m]=nfac[j]; j--; } while ( j != 0 ) ; }

//  find maximum prime factor
  for ( int i=1; i<nf; i++)
    if ( nfac[i] > maxf) { maxf=nfac[i]; }

//  compute Fourier transform
  int kk, k1, k2, k3=0, k4;
  double sd, cd;
  double aa, bb, ak, bk, aj, bj;
  double ajm, ajp, akm, akp, bjm, bjp, bkm, bkp;

// allocate temporary arrays, here np[] is not quaranteed to have enough size
// but it will be adjusted if data will not fit default size
  double *at, *ck, *bt, *sk;
  at=new double[maxf+1];
  bt=new double[maxf+1];
  ck=new double[maxf+1];
  sk=new double[maxf+1];

  int *np;
  np=new int[maxp+1];

L100:
    sd=radf/double(kspan);
    cd=sin(sd); cd*=cd; cd*=2.;
    sd=sin(sd+sd);
    kk=0;
    i++;

// transform for factor of 2  (including rotation factor)
    if ( nfac[i] != 2 ) goto L400;
         kspan=kspan/2;
         k1=kspan+2-1;
         do
         {  do
            {  k2=kk+kspan;
               ak=a[k2];           bk=b[k2];
               a[k2]=a[kk]-ak;     b[k2]=b[kk]-bk;
               a[kk]=a[kk]+ak;     b[kk]=b[kk]+bk;
               kk=k2+kspan;

             } while (kk <= nn-1);

             kk=kk-nn;

          } while ( kk <= jc-1);
          if ( kk > (kspan -1) ) goto L800 ;

          do
          {  c1=1.-cd;
             s1=sd;
             do
             { do
                { do
                  { k2=kk+kspan;
                    ak=a[kk]-a[k2];        bk=b[kk]-b[k2];
                    a[kk]=a[kk]+a[k2];     b[kk]=b[kk]+b[k2];
                    a[k2]=c1*ak-s1*bk;     b[k2]=s1*ak+c1*bk;
                    kk=k2+kspan;

                  } while (kk < nt-1);

                  k2=kk-nt;
                  c1=-c1;
                  kk=k1-k2-1;

                } while (kk > k2);

                ak=c1-(cd*c1+sd*s1);
                s1=(sd*c1-cd*s1)+s1;
                c1=2.-(ak*ak+s1*s1);
                s1=c1*s1;
                c1=c1*ak;
                kk+=jc;

              } while ( kk < k2 );

              k1+=inc; k1+=inc;
              kk=(k1-kspan+1)/2+jc-1;

            } while (kk <= jc+jc-1);

     goto L100;

//  transform for factor of 3 (optional code)
L320:
  do
  { do
    { k1=kk+kspan;            k2=k1+kspan;
      ak=a[kk];               bk=b[kk];
      aj=a[k1]+a[k2];         bj=b[k1]+b[k2];
      a[kk]=ak+aj;            b[kk]=bk+bj;
      ak=-0.5*aj+ak;          bk=-0.5*bj+bk;
      aj=(a[k1]-a[k2])*s120;  bj=(b[k1]-b[k2])*s120;
      a[k1]=ak-bj;            b[k1]=bk+aj;
      a[k2]=ak+bj;            b[k2]=bk-aj;
      kk=k2+kspan;

    } while ( kk < nn -1 );

      kk=kk-nn;

  }   while ( kk <= kspan -1 );

      goto L700;

//  transform for factor of 4
L400:
    if ( nfac[i] != 4 ) goto L600;
      kspnn=kspan;
      kspan=kspan/4;
L410:
      c1=1.;
      s1=0.;
L420:
      k1=kk+kspan;      k2=k1+kspan;      k3=k2+kspan;
      akp=a[kk]+a[k2];  akm=a[kk]-a[k2];
      ajp=a[k1]+a[k3];  ajm=a[k1]-a[k3];
      a[kk]=akp+ajp;
      ajp=akp-ajp;
      bkp=b[kk]+b[k2];  bkm=b[kk]-b[k2];
      bjp=b[k1]+b[k3];  bjm=b[k1]-b[k3];
      b[kk]=bkp+bjp;
      bjp=bkp-bjp;
      if ( isn < 0 ) goto L450;
      akp=akm-bjm;      akm=akm+bjm;
      bkp=bkm+ajm;      bkm=bkm-ajm;
      if ( s1 == 0 ) goto L460;
L430:
      a[k1]=akp*c1-bkp*s1;    b[k1]=akp*s1+bkp*c1;
      a[k2]=ajp*c2-bjp*s2;    b[k2]=ajp*s2+bjp*c2;
      a[k3]=akm*c3-bkm*s3;    b[k3]=akm*s3+bkm*c3;
      kk=k3+kspan;

      if ( kk <= nt-1 ) goto L420;
L440:
      c2=c1-(cd*c1+sd*s1);
      s1=(sd*c1-cd*s1)+s1;
      c1=2.0-(c2*c2+s1*s1);
      s1*=c1;
      c1*=c2;
      c2=c1*c1-s1*s1;
      s2=2.0*c1*s1;
      c3=c2*c1-s2*s1;
      s3=c2*s1+s2*c1;
      kk-=nt; kk+=jc;

      if ( kk <= kspan-1) goto L420;
      kk-=kspan; kk+=inc;

      if ( kk <= jc-1 ) goto L410;
      if ( kspan == jc ) goto L800;

      goto L100;

L450:
      akp=akm+bjm;   akm=akm-bjm;
      bkp=bkm-ajm;   bkm=bkm+ajm;

      if ( s1 != 0) goto L430;

L460:
      a[k1]=akp;     b[k1]=bkp;
      a[k2]=ajp;     b[k2]=bjp;
      a[k3]=akm;     b[k3]=bkm;
      kk=k3+kspan;

      if ( kk <= nt-1 ) goto L420;

      goto L440;

//  transform for factor of 5 (optional code)
L510:
      c2=c72*c72-s72*s72;
      s2=2.0*c72*s72;
  do
  { do
    { k1=kk+kspan; k2=k1+kspan; k3=k2+kspan;  k4=k3+kspan;
      akp=a[k1]+a[k4];          akm=a[k1]-a[k4];
      bkp=b[k1]+b[k4];          bkm=b[k1]-b[k4];
      ajp=a[k2]+a[k3];          ajm=a[k2]-a[k3];
      bjp=b[k2]+b[k3];          bjm=b[k2]-b[k3];
      aa=a[kk];                 bb=b[kk];
      a[kk]=aa+akp+ajp;         b[kk]=bb+bkp+bjp;
      ak=akp*c72+ajp*c2+aa;     bk=bkp*c72+bjp*c2+bb;
      aj=akm*s72+ajm*s2;        bj=bkm*s72+bjm*s2;
      a[k1]=ak-bj;              a[k4]=ak+bj;
      b[k1]=bk+aj;              b[k4]=bk-aj;
      ak=akp*c2+ajp*c72+aa;     bk=bkp*c2+bjp*c72+bb;
      aj=akm*s2-ajm*s72;        bj=bkm*s2-bjm*s72;
      a[k2]=ak-bj;              a[k3]=ak+bj;
      b[k2]=bk+aj;              b[k3]=bk-aj;
      kk=k4+kspan;
    }                       while ( kk < nn-1 );
      kk=kk-nn;
  }                         while ( kk <= kspan-1 ) ;
      goto L700;

//  transform for odd factors

L600:
      k=nfac[i];
      kspnn=kspan;
      kspan=kspan/k;
      if ( k == 3 ) goto L320;
      if ( k == 5 ) goto L510;
      if ( k == jf ) goto L640;
      jf=k;
      s1=rad/double(k);
      c1=cos(s1);       s1=sin(s1);
      if ( jf > maxf ) goto L998 ;
      ck[jf]=1.0;       sk[jf]=0.0;
      j=1;
   do
   {  ck[j]=ck[k]*c1+sk[k]*s1;
      sk[j]=ck[k]*s1-sk[k]*c1;
      k=k-1;
      ck[k]=ck[j];
      sk[k]=-sk[j];
      j=j+1;
    }  while ( j < k );
L640:
  do
  { do
    { k1=kk;
      k2=kk+kspnn;
      aa=a[kk];   bb=b[kk];
      ak=aa;      bk=bb;
      j=1;
      k1=k1+kspan;

      do
      { k2=k2-kspan;
        j++;
        at[j]=a[k1]+a[k2];
        ak=at[j]+ak;
        bt[j]=b[k1]+b[k2];
        bk=bt[j]+bk;
        j++;
        at[j]=a[k1]-a[k2];
        bt[j]=b[k1]-b[k2];
        k1=k1+kspan;
      }  while ( k1 < k2 );

      a[kk]=ak;
      b[kk]=bk;
      k1=kk;
      k2=kk+kspnn;
      j=1;

      do
      { k1+=kspan;   k2-=kspan;
        jj=j;
        ak=aa;       bk=bb;
        aj=0.0;      bj=0.0;
        k=1;

        do
        { k++;
          ak=at[k]*ck[jj]+ak;
          bk=bt[k]*ck[jj]+bk;
          k++;
          aj=at[k]*sk[jj]+aj;
          bj=bt[k]*sk[jj]+bj;
          jj+=j;
          if ( jj > jf ) jj-=jf;

        }                               while ( k < jf );

        k=jf-j;
        a[k1]=ak-bj;  b[k1]=bk+aj;
        a[k2]=ak+bj;  b[k2]=bk-aj;
        j++;

      }                                 while( j < k );

      kk=kk+kspnn;
    }                                   while ( kk <= nn-1 );
      kk-=nn;
  }                                     while ( kk <= kspan-1 );
//  multiply by rotation factor (except for factors of 2 and 4);
L700:
  if ( i == m ) goto L800;
  kk=jc;
  do
  { c2=1.0-cd;
    s1=sd;

    do
    { c1=c2;
      s2=s1;
      kk+=kspan;

      do
      { do
        { ak=a[kk];
          a[kk]=c2*ak-s2*b[kk];
          b[kk]=s2*ak+c2*b[kk];
          kk+=kspnn;
        }                           while ( kk <= nt-1 );

        ak=s1*s2;
        s2=s1*c2+c1*s2;
        c2=c1*c2-ak;
        kk=kk-nt+kspan;
      }                             while ( kk <= kspnn-1 );

      c2=c1-(cd*c1+sd*s1);
      s1=s1+(sd*c1-cd*s1);
      c1=2.0-(c2*c2+s1*s1);
      s1=c1*s1;
      c2=c1*c2;
      kk=kk-kspnn+jc;

    }                              while ( kk <= kspan -1 );

      kk=kk-kspan+jc+inc;

  }                                while ( kk <= jc+jc-1 );
      goto L100;

//  permute the results to normal order---done in two stages
//  permutation for square factors of n

L800:
      np[1]=ks;
      if ( kt == 0) goto L890;
      k=kt+kt+1;
      if ( m < k) k=k-1;
      j=1;
      np[k+1]=jc;
    do
    { np[j+1]=np[j]/nfac[j];
      np[k]=np[k+1]*nfac[j];
      j++;
      k--;
    } while ( j < k );

      k3=np[k+1]-1;
      kspan=np[2];
      kk=jc;
      k2=kspan;
      j=1;
      if ( n != ntot ) goto L850;
//  permutation for single-variate transform (optional code)
L820:
      do
      { ak=a[kk];
        a[kk]=a[k2];
        a[k2]=ak;
        bk=b[kk];
        b[kk]=b[k2];
        b[k2]=bk;
        kk+=inc;
        k2+=kspan;
      } while ( k2 < ks-1 );
L830:
    do
    { k2-=np[j];
      j++;
      k2+=np[j+1];
    } while ( k2 > np[j]-1 );
      j=1;
L840: if ( kk < k2 ) goto L820;
      kk+=inc;
      k2+=kspan;
      if ( k2 < ks-1 ) goto L840;
      if ( kk < ks-1 ) goto L830;
      jc=k3+1;
      goto L890;
//  permutation for multivariate transform;
L850:
       do
       { do
         { k=kk+jc+1;
           do
           { ak=a[kk];  a[kk]=a[k2];  a[k2]=ak;
             bk=b[kk];  b[kk]=b[k2];  b[k2]=bk;
             kk+=inc;   k2+=inc;
           }                         while ( kk < k-1 );

           kk=kk+ks-jc;
           k2=k2+ks-jc;
         }                           while ( kk < nt-1 );

         k2=k2-nt+kspan;
         kk=kk-nt+jc;
       }                             while ( k2 < ks-1 );
 L870:
       do { k2-=np[j]; j++;  k2+=np[j+1]; }  while ( k2 > np[j]-1 );

       j=1;
 L880:
       if ( kk < k2 ) goto L850;
       kk+=jc;
       k2+=kspan;
       if ( k2 < ks-1 ) goto L880;
       if ( kk < ks-1 ) goto L870;
       jc=k3+1;
 L890:
       if ( 2*kt+1 >= m ) {delete [] np; delete [] at; delete [] bt;
                           delete [] ck ; delete [] sk; return;}
       kspnn=np[kt+1];
//  permutation for square-free factors of n;
      j=m-kt;
      nfac[j+1]=1;
//L900:
      do { nfac[j]=nfac[j]*nfac[j+1];  j--; } while ( j != kt );

      kt++;
      nn=nfac[kt]-1;
      if ( nn > maxp )                          // was goto L998;
      { maxp=nn;
        delete [] np; np=new int[maxp+1];  }

      jj=0;
      j=0;
      goto L906;
L902:
      jj=jj-k2-1;
      k2=kk;
      k++;
      kk=nfac[k]-1;
L904:
      jj=kk+jj+1;
      if ( jj-1 >= k2 ) goto L902;
      np[j]=jj;
L906:
      k2=nfac[kt]-1;
      k=kt+1;
      kk=nfac[k]-1;
      j++;
      if ( j <= nn ) goto L904;
//  determine the permutation cycles of length greater than 1;
      j=0;
      goto L914;
L910:
      do { k=kk+1;  kk=np[k]-1;  np[k]=-kk-1; } while ( kk != j-1 );

      k3=kk;
L914:
      do { j++;  kk=np[j]-1; } while ( kk < -1 );

      if ( kk != j-1 ) goto L910;
      np[j]=-j;
      if ( j != nn ) goto L914;
      maxf=inc*maxf;
//  reorder a and b, following the permutation cycles;
      goto L950;
L924:
      do
      { do { j--; } while ( np[j] < 0 );

        jj=jc;
        do
        { kspan=jj;
          if ( jj > maxf ) kspan=maxf;
          jj-=kspan;
          k=np[j];
          kk=jc*k+i+jj-1;
          k1=kk+kspan;
          k2=0;
          do { at[k2]=a[k1]; bt[k2]=b[k1]; k2++; k1-=inc; } while ( k1 != kk );

          do
          { k1=kk+kspan;
            k2=k1-jc*(k+np[k]);
            k=-np[k];
            do { a[k1]=a[k2]; b[k1]=b[k2]; k1-=inc; k2-=inc; } while( k1 != kk );
            kk=k2;
          } while ( k != j );

          k1=kk+kspan;
          k2=0;
          do { a[k1]=at[k2]; b[k1]=bt[k2]; k2++; k1-=inc; } while ( k1 != kk );

        } while ( jj != 0 );
      } while ( j != 1 );
L950:
      j=k3+2;
      nt=nt-kspnn;
      i=nt-inc+1;
      if ( nt >= 0 ) goto L924;
      delete [] np; delete [] at; delete [] bt; delete [] ck ; delete [] sk;
      return;
//  error finish, insufficient array storage;
L998:
      delete [] np; delete [] at; delete [] bt; delete [] ck ; delete [] sk;
      isn=0;
      cout <<"Error: array bounds exceeded within subroutine wavefft.\n";
      return;
}

