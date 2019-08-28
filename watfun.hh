// Wavelet Analysis Tool
// S.Klimenko, University of Florida
// library of general functions

#ifndef WATFUN_HH
#define WATFUN_HH

//#include <stream>
#include <iostream>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592653589793
#define speedlight 299792458.0

// WAT functions

// Calculates polynomial interpolation coefficient using Lagrange formula.
inline double Lagrange(const int n, const int i, const double x)
{
    double c = 1.;
    double xn = x+n/2.-0.5;	// shift to the center of interpolation interval

    for(int j=0; j<n; j++) 
       if(j!=i) c *= (xn-j)/(i-j);

    return c;
}


// Nevill's polynomial interpolation.
// x - polynom argument (x stride is allways 1) 
// n - number of interpolation samples
// p - pointer to a sample with x=0.
// q - double array of length n
template<class DataType_t>
inline double Nevill(const double x0,
		    int n,   
		    DataType_t* p,
		    double* q)
{
   register int i;
   register double x = x0;
   register double xm = 0.5;

   n--;
   *q = *p;

   for(i=0; i<n; i++)
      q[i] = p[i] + (x--)*(p[i+1]-p[i]);

   while(--n >= 1){
      x = x0;

      q[0] += xm*(x--)*(q[1]-q[0]);
      if(n == 1) goto M0;
      q[1] += xm*(x--)*(q[2]-q[1]);
      if(n == 2) goto M0;
      q[2] += xm*(x--)*(q[3]-q[2]);
      if(n == 3) goto M0;
      q[3] += xm*(x--)*(q[4]-q[3]);
      if(n == 4) goto M0;
      q[4] += xm*(x--)*(q[5]-q[4]);
      if(n == 5) goto M0;
      q[5] += xm*(x--)*(q[6]-q[5]);
      if(n == 6) goto M0;

      for(i=6; i<n; i++)
	 q[i] += xm*(x--)*(q[i+1]-q[i]);

M0:   xm /= (1.+xm);
   }

   return *q;
}

// calculate significance for sign cross-correlation
inline double signPDF(const size_t m, const size_t k) 
{
   size_t i;
   double pdf = 0.;
   size_t n = m+k;
   double rho = (double(m)-double(k))/n;

   if(n < 1) return 0.;
   if(n < 100) {
      for(i=1; i<k+1; i++){ pdf -= log(double(m+i)/double(i)); }
      pdf -= log(double(n))-(n+1)*log(2.);
      pdf -= log(sqrt(2.*PI/n));  // normalization from Gaussian approximation
   }
   else pdf = n*rho*rho/2.;
   return pdf;
}


// survival probability for Gamma distribution
// calculates integral I=int(x*x^n*exp(-x)) from Y to infinity
// returns confidence = -log(I/Gamma(n))
// Y - low integration limit
// n - Gamma function parameter
inline double gammaCLa(double Y, int n){   
  double y = Y;
  double s = 1.;
  //  if(Y<0.) return 0.;
  for(int k=1; k<n; k++){ 
    s += y; 
    if((y*=Y/(k+1))>1.e290) break; 
  }
  return Y-log(s);
}

inline double gammaCL(double x, double n){   
  return x>=n ? x-n + (2+(n-1)/sqrt(2.))/(n+1)-(n-sqrt(n/4.)-1)*log(x/n) : 
    pow(x/n,sqrt(n))*(2+(n-1)/sqrt(2.))/(n+1);
}

// survival probability for LogNormal distribution
// calculates integral I=int(Pln(x)) from Y to infinity
// returns one side confidence = -log(I)
// Y - low integration limit
// p - x peak
// s - sigma
// a - asymmetry

inline double logNormCL(double Y, double p=0., double s=1., double a=0.0001){   
  if(a <= 0.) a = 0.0001;
  double z = a*sqrt(log(4.));
  z = sinh(z)/z; 
  z = 1+a*z*(Y-p)/s;
  z = z>0 ? fabs(log(z)/a-a) : 0.; 
  return erfc(z/sqrt(2.))/2.;
}

// calculate value of argument from LogNormal and Normal survival probability 
// calculates inverse of integral I=int(P(x) dx) from Y to infinity
// returns value of the argument Y
// CL - single side probability
// p - x peak
// s - sigma
// a - asymmetry

inline double logNormArg(double CL, double p=0., double s=1., double a=0.){   
  double z = fabs(sqrt(-2.14*log(2.*CL))-0.506);
  if(a <= 0.) return p+z*s;
  double f = a*sqrt(log(4.));
  z = (exp((z+a)*a)-1)/a;
  return p+z*s*f/sinh(f);
}


// Calculates polynomial interpolation coefficients using Lagrange formula.
// gap is allouded in the middle of interval, like
// n=13, m=5:  x x x x o o o o o x x x x 
inline void fLagrange(int n, int m, double* c)
{
    if(!(n&1)) n++;
    if(!(m&1)) m++;
    int i,j;

    for(i=0; i<n; i++) c[i]=0.; 

    for(i=0; i<n; i++) {
       if(abs(n/2-i)<=m/2) continue;
       c[i] = 1.;
       for(j=0; j<n; j++) {
	  if(abs(n/2-j)<=m/2) continue;
          if(j!=i) c[i] *= double(n/2-j)/(i-j);
       }
    }
    return;
}

// complete gamma function
inline double Gamma(double r)
{
   if(r<=0.) return 0.;
   double p0 = 1.000000000190015;
   double p1 = 76.18009172947146;
   double p2 = -86.50532032941677;
   double p3 = 24.01409824083091;
   double p4 =  -1.231739572450155;
   double p5 = 1.208650973866179e-3;
   double p6 = -5.395239384953e-6;
   double g = p0+p1/(r+1)+p2/(r+2)+p3/(r+3)+p4/(r+4)+p5/(r+5)+p6/(r+6);
   return g*pow(r+5.5,r+0.5)*exp(-(r+5.5))*sqrt(2*PI)/r;
}

// incomplete regularized gamma function
inline double Gamma(double r, double x)
{
   double a = r>0 ? 1/r : 0;
   double b = a;
   int n = int(10*x-r);
   if(n<100) n = 100;
   for(int i=1; i<n; i++) {
      a *= x/(r+i);
      b += a;
   }
   return b*exp(log(x)*r-x)/Gamma(r);
}

// inverse incomplete regularized gamma function (approximation)
inline double iGamma(double r, double p)
{
   double x = 700;
   double P = 1.;
   while(P>p) { x *= 0.999; P=Gamma(r,x); }
   return 2*x;
}


#endif // WATFUN_HH

















