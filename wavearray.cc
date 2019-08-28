/*-------------------------------------------------------
 * Package: 	Wavelet Analysis Tool
 * generic data container to use with DMT and ROOT
 * File name: 	wavearray.cc
 *-------------------------------------------------------
*/
// wavearray class is the base class for wavelet data amd methods . 

#include <time.h>
#include <iostream>
//#include <sstream>
//#include <strstream>

//#include "PConfig.h"
#ifdef __GNU_STDC_OLD
#include <gnusstream.h>
#else
#include <sstream>
#endif

#include "wavearray.hh"
#include "wavefft.hh"

extern "C" {
   typedef int (*qsort_func) (const void*, const void*);
}

using namespace std;

//: Default constructor
template<class DataType_t>
wavearray<DataType_t>::wavearray() : 
data(NULL), Size(0), Rate(1.), Start(0.) 
{   
   Slice = std::slice(0,0,0); 
}

//: allocates a data array with N=n elements.
template<class DataType_t>
wavearray<DataType_t>::wavearray(int n) : 
Rate(1.), Start(0.)
{ 
  if (n <= 0 ) n = 1;
  data = (DataType_t *)malloc(n*sizeof(DataType_t));
  Size = n;
  Slice = std::slice(0,n,1);
  *this = 0;
}

//: copy constructor 
template<class DataType_t>
wavearray<DataType_t>::wavearray(const wavearray<DataType_t>& a) :
data(NULL), Size(0), Rate(1.), Start(0.)
{ *this = a; }

// explicit construction from array
template<class DataType_t> template<class T> 
wavearray<DataType_t>::
wavearray(const T *p, unsigned int n, double r) : 
data(NULL), Size(0), Rate(1.), Start(0.)
{ 
   unsigned int i;
   if(n != 0 && p != NULL){
      data = (DataType_t *)malloc(n*sizeof(DataType_t));
      for (i=0; i < n; i++) data[i] = p[i];
      Size = n;
      Rate = r;
   }
   Slice = std::slice(0,n,1);
}

//: destructor
template<class DataType_t>
wavearray<DataType_t>::~wavearray()
{ 
   free(data);
}

//: operators =

template<class DataType_t>
wavearray<DataType_t>& wavearray<DataType_t>::
operator=(const wavearray<DataType_t>& a)
{
   if (this==&a) return *this;

   unsigned int i;
   unsigned int N = a.Slice.size();
   unsigned int m = a.Slice.stride();

   if (N>0 && a.data) {
      register DataType_t *pm  = a.data + a.Slice.start();
      wavearray<DataType_t>::resize(N);

      for (i=0; i < N; i++) { data[i] = *pm; pm+=m; }

      if(a.rate()>0.) 
	 start(a.start() + a.Slice.start()/a.rate());
      else          
	 start(a.start());

      rate(a.rate());
      Slice = std::slice(0,  size(),1);
      const_cast<wavearray<DataType_t>&>(a).Slice = std::slice(0,a.size(),1);
   }
   else if(data) {
      free(data);
      data = NULL;
      Size = 0;
      Start = a.start();
      Rate = a.rate();
      Slice = std::slice(0,0,0);
   }

   return *this;
}

template<class DataType_t>
wavearray<DataType_t>& wavearray<DataType_t>::
operator<<(wavearray<DataType_t>& a)
{
   unsigned int i;
   unsigned int N = limit(a);
   unsigned int n = Slice.stride();
   unsigned int m = a.Slice.stride();
   register DataType_t *p = a.data + a.Slice.start();

   if(size())
      for (i=Slice.start(); i<N; i+=n){ data[i]  = *p; p += m; }
   
     Slice = std::slice(0,  size(),1);
   a.Slice = std::slice(0,a.size(),1);
   return *this;
}

template<class DataType_t>
wavearray<DataType_t>& wavearray<DataType_t>::
operator=(const DataType_t c)
{
   unsigned int i;
   unsigned int n = Slice.stride();
   unsigned int N = limit();

   if(size())
      for (i=Slice.start(); i < N; i+=n) data[i]  = c;

   Slice = std::slice(0,  size(),1);
   return *this;
}

//: operators +=

template<class DataType_t>
wavearray<DataType_t>& wavearray<DataType_t>::
operator+=(wavearray<DataType_t> &a)
{
   unsigned int i;
   unsigned int N = limit(a);
   unsigned int n = Slice.stride();
   unsigned int m = a.Slice.stride();
   register DataType_t *p = a.data + a.Slice.start();

   if(size())
      for (i=Slice.start(); i<N; i+=n){ data[i]  += *p; p += m; }
   
     Slice = std::slice(0,  size(),1);
   a.Slice = std::slice(0,a.size(),1);
   return *this;
}

template<class DataType_t>
wavearray<DataType_t>& wavearray<DataType_t>::
operator+=(const DataType_t c)
{
   unsigned int i;
   unsigned int n = Slice.stride();
   unsigned int N = limit();

   if(size())
      for (i=Slice.start(); i < N; i+=n) data[i]  += c;

   Slice = std::slice(0,  size(),1);
   return *this;
}

//: operators -=

template<class DataType_t>
wavearray<DataType_t>& wavearray<DataType_t>::
operator-=(wavearray<DataType_t> &a)
{
   unsigned int i;
   unsigned int N = limit(a);
   unsigned int n = Slice.stride();
   unsigned int m = a.Slice.stride();
   register DataType_t *p = a.data + a.Slice.start();

   if(size())
      for (i=Slice.start(); i<N; i+=n){ data[i]  -= *p; p += m; }
   
     Slice = std::slice(0,  size(),1);
   a.Slice = std::slice(0,a.size(),1);
   return *this;
}

template<class DataType_t>
wavearray<DataType_t>& wavearray<DataType_t>::
operator-=(const DataType_t c)
{
   unsigned int i;
   unsigned int n = Slice.stride();
   unsigned int N = limit();

   if(size())
      for (i=Slice.start(); i < N; i+=n) data[i]  -= c;

   Slice = std::slice(0,  size(),1);
   return *this;
}

//: multiply all elements of data array by constant
template<class DataType_t>
wavearray<DataType_t>& wavearray<DataType_t>::
operator*=(const DataType_t c)
{
   unsigned int i;
   unsigned int n = Slice.stride();
   unsigned int N = limit();

   if(size())
      for (i=Slice.start(); i < N; i+=n) data[i]  *= c;

   Slice = std::slice(0,  size(),1);
   return *this;
}

// scalar production
template<class DataType_t>
wavearray<DataType_t>& wavearray<DataType_t>::
operator*=(wavearray<DataType_t>& a)
{
   unsigned int i;
   unsigned int N = limit(a);
   unsigned int n = Slice.stride();
   unsigned int m = a.Slice.stride();
   register DataType_t *p = a.data + a.Slice.start();

   if(size())
      for (i=Slice.start(); i<N; i+=n){ data[i]  *= *p; p += m; }
   
     Slice = std::slice(0,  size(),1);
   a.Slice = std::slice(0,a.size(),1);
   return *this;
}

// operator[](const std::slice &)
template<class DataType_t>
wavearray<DataType_t>& wavearray<DataType_t>::
operator[](const std::slice& s)
{
   Slice = s;
   if(limit() > size()){
      cout << "wavearray::operator[slice]: Illegal argument "<<limit()<<" "<<size()<<"\n";
      Slice = std::slice(0,size(),1);
   }
   return *this;
}

// operator[](const std::slice &)
template<class DataType_t>
DataType_t & wavearray<DataType_t>::
operator[](const unsigned int n)
{
   if(n >= size()){
      cout << "wavearray::operator[int]: Illegal argument\n";
      return data[0];
   }
   return data[n];
}

//: Dumps data array to file *fname in ASCII format.
template<class DataType_t>
void wavearray<DataType_t>::Dump(const char *fname, int app)
{
  int n=size();
  char mode[3] = "w";
  if (app == 1) strcpy(mode, "a");

  FILE *fp;

  if ( (fp = fopen(fname, mode)) == NULL ) {
     cout << " Dump() error: cannot open file " << fname <<". \n";
     return;
  };

  if(app == 0) {
    fprintf( fp,"# start time: -start %lf \n", this->Start );
    fprintf( fp,"# sampling rate: -rate %lf \n", this->Rate );
    fprintf( fp,"# number of samples: -size %d \n", (int)this->Size );
  }

  for (int i = 0; i < n; i++) fprintf( fp,"%e ", (float)data[i]);
  fclose(fp); 
}
  
//: Dumps data array to file *fname in binary format.
template<class DataType_t>
void wavearray<DataType_t>::DumpBinary(const char *fname, int app)
{
  int n = size() * sizeof(DataType_t);
  char mode[5];
  strcpy (mode, "w");

  if (app == 1) strcpy (mode, "a");

  FILE *fp;

  if ( (fp=fopen(fname, mode)) == NULL ) {
     cout << " DumpBinary() error : cannot open file " << fname <<". \n";
     return ;
  }

  fwrite(data, n, 1, fp);
  fclose(fp);
}

//: Dumps data array to file *fname in binary format as "short".
template<class DataType_t> 
void wavearray<DataType_t>::DumpShort(const char *fname, int app)
{
  int n = size();
  char mode[5] = "w";
  if (app == 1) strcpy (mode, "a");

  FILE *fp;
  if ( (fp = fopen(fname, mode)) == NULL ) {
    cout << " DumpShort() error : cannot open file " << fname <<". \n";
    return;
  }

  short *dtemp;
  dtemp=new short[n];

  for ( int i=0; i<n; i++ ) dtemp[i]=short(data[i]) ;

  n = n * sizeof(short);

  fwrite(dtemp, n, 1, fp);
  fclose(fp);
  delete [] dtemp;
}

//: Read data from file in binary format.
template<class DataType_t> 
void wavearray<DataType_t>::ReadBinary(const char *fname, int N)
{
   int step = sizeof(DataType_t);
   FILE *fp;
   double d;

   if ( (fp=fopen(fname,"r")) == NULL ) {
      cout << " ReadBinary() error : cannot open file " << fname <<". \n";
      exit(1);
   }


   if(N == 0){              // find the data length
      while(!feof(fp)){ 
	 fread(&d,step,1,fp);
	 N++;
      }
      N--;
      rewind(fp);
   }	

   if(int(size()) != N) resize(N);
   int n = size() * sizeof(DataType_t);

   fread(data, n, 1, fp);    // Reading binary record
   fclose(fp);
}

//: Read data from file as short.
template<class DataType_t> 
void wavearray<DataType_t>::ReadShort(const char *fname)
{
  short *dtmp;
  dtmp = new short[size()];
  int n = size() * sizeof(short);
  FILE *fp;

  if ( (fp=fopen(fname,"r")) == NULL ) {
     cout << " ReadShort() error : cannot open file " << fname <<". \n";
     return;
  }

  cout << " Reading binary record, size="<< n <<"\n";

  fread(dtmp, n, 1, fp);
  for (unsigned int i = 0; i < size(); i++) 
     data[i] = DataType_t(dtmp[i]);
  fclose(fp);
  delete [] dtmp;
}

//: resizes data array to a new length n.
template<class DataType_t> 
void wavearray<DataType_t>::resize(unsigned int n)
{
   DataType_t *p = NULL;
   if(n==0){
     if(data) free(data);
     data = NULL;
     Size = 0;
     Slice = std::slice(0,0,0); 
     return;
   }

   p = (DataType_t *)realloc(data,n*sizeof(DataType_t));

   if(p){ 
      data = p;
      Size = n;
      Slice = std::slice(0,n,1);
   }
   else
      cout<<"wavearray::resize(): memory allocation failed.\n";
}

/************************************************************************
 * Creates new data set by resampling the original data from "a"        *
 * with new sample frequency "f". Uses polynomial interpolation scheme  *
 * (Lagrange formula) with np-points, "np" must be even number, by      *
 * default np=6 (corresponds to 5-th order polynomial interpolation).   *
 * This function calls wavearray::resize() function to adjust            *
 * current data array if it's size does not exactly match new number    *
 * of points.                                                           *
 ************************************************************************/

template<class DataType_t> 
void wavearray<DataType_t>::
resample(wavearray<DataType_t> const &a, double f, int nF)
{
   int nP = nF;
   if(nP<=1) nP = 6;
   if(int(a.size())<nP) nP = a.size();
   nP = (nP>>1)<<1;

   int N;
   int nP2 = nP/2;

   const DataType_t *p = a.data;

   register int i;
   register int iL;
   register double x;
   double *temp=new double[nF];
   
   rate(f);
   double ratio = a.rate() / rate();
   N = int(a.size()/ratio + 0.5); 
   
   if ( (int)size() != N )  resize(N); 
   
// left border
   
   int nL = int(nP2/ratio);

   for (i = 0; i < nL; i++)
      data[i] = (DataType_t)Nevill(i*ratio, nP, p, temp);

// in the middle of array

   int nM = int((a.size()-nP2)/ratio);  
   if(nM < nL) nM = nL;

   if(nM&1 && nM>nL) { 
      x = nL*ratio;
      iL = int(x) - nP2 + 1 ;
      data[i] = (DataType_t)Nevill(x-iL, nP, p+iL, temp);
      nL++;
   }
   for (i = nL; i < nM; i+=2) { 
      x = i*ratio;
      iL = int(x) - nP2 + 1 ;
      data[i] = (DataType_t)Nevill(x-iL, nP, p+iL, temp);
      x += ratio;
      iL = int(x) - nP2 + 1 ;
      data[i+1] = (DataType_t)Nevill(x-iL, nP, p+iL, temp);
   }

// right border

   int nR = a.size() - nP;
   p += nR;
   for (i = nM; i < N; i++)
      data[i] = (DataType_t)Nevill(i*ratio-nR, nP, p, temp);

   delete [] temp;
}

template<class DataType_t> 
void wavearray<DataType_t>::resample(double f, int nF)
{
   wavearray<DataType_t> a;
   a = *this;
   resample(a,f,nF);
}


template<class DataType_t> 
void wavearray<DataType_t>::
Resample(wavearray<DataType_t> const &a, double f, int np)
{
  int i1, N, n1, n2, n, np2 = np/2;
  double s, x, *c, *v;
  c=new double[np];
  v=new double[np];

  rate(f);
  double ratio = a.rate() / rate();
  n = a.size();
  N = int(n/ratio  + 0.5); 

  if ( (int)size() != N )  resize(N); 

// calculate constant filter part c(k) = -(-1)^k/k!/(np-k-1)!
  for (int i = 0; i < np; i++) {
    int m = 1;
    for (int j = 0; j < np; j++)  if (j != i) m *= (i - j);
    c[i] = 1./double(m);
  }

  for (int i = 0; i < N; i++)
  { 
    x = i*ratio;        

    i1 = int(x);
    x = x - i1 + np2 - 1;

// to treat data boundaries we need to calculate critical numbers n1 and n2
    n1 = i1 - np2 + 1;
    n2 = i1 + np2 + 1 - n;

// here we calculate the sum of products like h(k)*a(k), k=0..np-1,
// where h(k) are filter coefficients, a(k) are signal samples
// h(k) = c(k)*prod (x-i), i != k, c(k) = -(-1)^k/k!/(np-k-1)! 

// get signal part multiplied by constant filter coefficients
    if ( n1 >= 0 )
      if ( n2 <= 0 ) 
// regular case - far from boundaries
        for (int j = 0; j < np; j++) v[j] = c[j]*a.data[i1 + j - np2 + 1];
      
      else {
// right border case
        x = x + n2;
        for (int j = 0; j < np; j++) v[j] = c[j]*a.data[n + j - np];
      }  
    else {
// left border case
      x = x + n1;
      for (int j = 0; j < np; j++)  v[j] = c[j]*a.data[j];
    }

// multiply resulted v[j] by x-dependent factors (x-k), if (k != j)
    for (int k = 0; k < np; k++) {
      for (int j = 0; j < np; j++) if (j != k) v[j] *= x; 
      x -= 1.; }

    s = 0.;
    for (int j = 0; j < np; j++)  s += v[j]; 
        
    data[i] = (DataType_t)s;
  }

  delete [] c;
  delete [] v;
}

/******************************************************************************
 * copies data from data array of the object wavearray a to *this
 * object data array. Procedure starts from the element "a_pos" in
 * the source array which is copied to the element "pos" in the
 * destination array. "length" is the number of data elements to copy.
 * By default "length"=0, which means copy all source data or if destination
 * array is shorter, copy data until the destination is full.
 *****************************************************************************/
template<class DataType_t> void wavearray<DataType_t>::
cpf(const wavearray<DataType_t> &a, int length, int a_pos, int pos)
{ 
   if (rate() != a.rate()){
      cout << "wavearray::cpf() warning: sample rate mismatch.\n";
      cout<<"rate out: "<<rate()<<"  rate in: "<<a.rate()<<endl;
   }

   if (length == 0 ) 
      length = ((size() - pos) < (a.size() - a_pos))? 
	 (size() - pos) : (a.size() - a_pos);

   if( length > (int)(size() - pos) ) length = size() - pos; 
   if( length > (int)(a.size() - a_pos) ) length = a.size() - a_pos; 

   for (int i = 0; i < length; i++)
     data[i + pos] = a.data[i + a_pos];

   rate(a.rate());
}

/******************************************************************************
 * Adds data from data array of the object wavearray a to *this
 * object data array. Procedure starts from the element "a_pos" in
 * the source array which is added starting from the element "pos" in the
 * destination array. "length" is the number of data elements to add.
 * By default "length"=0, which means add all source data or if destination
 * array is shorter, add data up to the end of destination.
 *****************************************************************************/
template<class DataType_t> void wavearray<DataType_t>::
add(const wavearray<DataType_t> &a, int length, int a_pos, int pos)
{
   if (rate() != a.rate())
      cout << "wavearray::add() warning: sample rate mismatch.\n";

   if (length == 0 )
      length = ((size() - pos) < (a.size() - a_pos))? 
	 (size() - pos) : (a.size() - a_pos);

   if( length > (int)( size()- pos) ) length = size() - pos; 
   if( length > (int)(a.size() - a_pos) ) length = a.size() - a_pos; 
   
   for (int i = 0; i < length; i++)
      data[i + pos] += a.data[i + a_pos];
}


/******************************************************************************
 * Subtracts data array of the object wavearray a from *this
 * object data array. Procedure starts from the element "a_pos" in
 * the source array which is subtracted starting from the element "pos" in
 * the destination array. "length" is number of data elements to subtract.
 * By default "length"=0, which means subtract all source data or if the 
 * destination array is shorter, subtract data up to the end of destination.
 ******************************************************************************/
template<class DataType_t> void wavearray<DataType_t>::
sub(const wavearray<DataType_t> &a, int length, int a_pos, int pos)
{
   if (rate() != a.rate())
      cout << "wavearray::sub() warning: sample rate mismatch.\n";

   if ( length == 0 )
      length = ((size() - pos) < (a.size() - a_pos))? 
	 (size() - pos) : (a.size() - a_pos);

   if( length > (int)(size() - pos) ) length = size() - pos; 
   if( length > (int)(a.size() - a_pos) ) length = a.size() - a_pos; 

   for (int i = 0; i < length; i++)
      data[i + pos] -= a.data[i + a_pos];
}

/*****************************************************************
 * append two wavearrays with the same rate ignoring start time 
 * of input array.
 *****************************************************************/
template<class DataType_t> size_t wavearray<DataType_t>::
append(const wavearray<DataType_t> &a)
{
   size_t n = this->size();
   size_t m = a.size();
   
   if (this->rate() != a.rate())
      cout << "wavearray::append() warning: sample rate mismatch.\n";

   if(m == 0 ) return this->size();
   this->resize(n+m);
   this->cpf(a,m,0,n);

   return n+m;
}

/*****************************************************************
 * append a data point 
 *****************************************************************/
template<class DataType_t> size_t wavearray<DataType_t>::
append(DataType_t a)
{
   size_t n = this->size();
   this->resize(n+1);
   this->data[n] = a;
   return n+1;
}


/*****************************************************************
 * Calculates Fourier Transform for real signal using
 * Fast Fourier Transform (FFT) algorithm. Packs resulting 
 * complex data into original array of real numbers.
 * Calls wavefft() function which is capable to deal with any
 * number of data points. FFT(1) means forward transformation,
 * which is default, FFT(-1) means inverse transformation.
 *****************************************************************/
template<class DataType_t> 
void wavearray<DataType_t>::FFT(int direction) 
{ 
  double *a, *b;
  int N = size();
  int n2 = N/2;

  a=new double[N];
  b=new double[N];

  switch (direction)
  { case 1:

// direct transform
    
    for (int i=0; i<N; i++) { a[i] = data[i]; b[i]=0.;}

    wavefft(a, b, N, N, N, -1);

// pack complex numbers to real array
    for (int i=0; i<n2; i++)
    { data[2*i]   = (DataType_t)a[i]/N;
      data[2*i+1] = (DataType_t)b[i]/N; 
    }

// data[1] is not occupied because imag(F[0]) is always zero and we
// store in data[1] the value of F[N/2] which is pure real for even "N"
      data[1] = (DataType_t)a[n2]/N;

// in the case of odd number of data points we store imag(F[N/2])
// in the last element of array data[N]
      if ((N&1) == 1) data[N-1] = (DataType_t)b[n2]/N;

      break;

    case -1:
// inverse transform

// unpack complex numbers from real array
      for (int i=1;i<n2;i++)
      { a[i]=data[2*i];
        b[i]=data[2*i+1];
        a[N-i]=data[2*i];
        b[N-i]=-data[2*i+1];
       }

      a[0]=data[0];
      b[0]=0.;

      if ((N&1) == 1)
        { a[n2]=data[1]; b[n2]=data[N-1]; }  // for odd n
      else
        { a[n2]=data[1]; b[n2]=0.; }

      wavefft(a, b, N, N, N, 1);             // call FFT for inverse tranform

      for (int i=0; i<N; i++)  
	 data[i] = (DataType_t)a[i];         // copy the result from array "a"
  }

  delete [] b;
  delete [] a;
}

/*************************************************************************
 * Stack generates wavearray *this by stacking data from wavearray td.
 * Input data are devided on subsets with with samples "length"
 * then they are added together. Average over the all subsets is saved in
 * *this. 
 *************************************************************************/
template<class DataType_t> 
double wavearray<DataType_t>::
Stack(const wavearray<DataType_t> &td, int length)
{
  double ss, s0, s2;
  rate(td.rate());
  int k = td.size()/length;
  int n = k*length;

  if (k == 0) {
    cout <<" Stack() error: data length too short to contain \n"
         << length << " samples\n";
    return 0.;
  }

  if (size() != (unsigned int)length) resize(length);

// sum (stack) all k periods of frequency f to produce 1 cycle
  s0 = 0.;
  for (int i = 0; i < length; i++) {
    ss = 0.;
    for (int j = i; j < n; j += length) ss += td.data[j];
    data[i] = (DataType_t)ss/k;
    s0 += ss;
  }
  s0 /= (k*length);

// remove constant displacement (zero frequency component) 
  s2 = 0.;
  for (int i = 0; i < length; i++) {
    data[i] -= (DataType_t)s0;
    s2 += data[i]*data[i];
  }
   s2 /= length;

   return s2;        // return stacked signal power (energy/N)
}

/*************************************************************************
 * Another version of Stack:
 * Input data (starting from sample "start") are devided on "k" subsets 
 * with sections "length"
 * Average over the all sections is saved in *this. 
 *************************************************************************/
template<class DataType_t> 
double wavearray<DataType_t>::
Stack(const wavearray<DataType_t> &td, int length, int start)
{
  double avr, rms;
  rate(td.rate());

  if(start+length > (int)td.size()) length = td.size()-start;

  int k = (size()<1) ? 0 : length/size();
  if (k == 0) {
    cout <<" Stack() error: data length too short to contain \n"
         << length << " samples\n";
    return 0.;
  }

  *this = 0;
  for (int i = 0; i<k; i++) add(td, size(), start+i*size());
  *this *= DataType_t(1./k);
  getStatistics(avr,rms);                                
  *this -= DataType_t(avr);

  return rms*rms;        // return stacked signal power
}

// Stack generates wavearray *this by stacking data from wavearray td.
// Input data are devided on subsets with with samples "length"
// then they are added together. Average over the all subsets is saved in
// *this.
template<class DataType_t> 
double wavearray<DataType_t>::
Stack(const wavearray<DataType_t> &td, double window)
{ 
   return this->Stack(td, int(td.rate() * window)); 
}

// wavearray mean
template<class DataType_t>
double wavearray<DataType_t>::mean(double sigma) const 
{
   register size_t i;
   register double x = 0;
   size_t N = (size()>>2)<<2;
   register DataType_t *p = data + size() - N;

   if(!size()) return 0.;

   if(sigma<=0.){

      for(i=0; i<size()-N; i++) { x += data[i]; }
      for(i=0; i<N; i+=4) 
	 x += p[i] + p[i+1] + p[i+2] + p[i+3];
      return double(x)/size();
   }
   else{
      int m=0;
      double a, y=0.;

      for(i=0; i<N; i+=4){ 
	 x += p[i] + p[i+1] + p[i+2] + p[i+3];
	 y += p[i]*p[i] + p[i+1]*p[i+1] + p[i+2]*p[i+2] + p[i+3]*p[i+3];
     }
      a = x/size();
      sigma *= sqrt(y/size()-x*x);
      x = 0.;
      for(i=0; i<N; i+=4){ 
	 if(fabs(p[i]-a) < sigma) { x += p[i]; m++; }
	 if(fabs(p[i]-a) < sigma) { x += p[i]; m++; }
	 if(fabs(p[i]-a) < sigma) { x += p[i]; m++; }
	 if(fabs(p[i]-a) < sigma) { x += p[i]; m++; }
      }
      return (m > 0) ? x/m : a;
   }

}


template<class DataType_t>
double wavearray<DataType_t>::mean(const std::slice &s)
{
   register size_t i;
   register double x = 0.;
   register DataType_t *p = data + s.start();
   size_t N = s.size();
   size_t m = (s.stride()<=0) ? 1 : s.stride();
   if(size()<limit(s)) N = (limit(s) - s.start() - 1)/m;
   for(i=0; i<N; i++) { x += *p; p += m; }
   return (N==0) ? 0. : x/N;
}

// running mean
template<class DataType_t>
void wavearray<DataType_t>::mean(double t, 
				 wavearray<DataType_t> *pm,
                                 bool clean, 
				 size_t skip)
{
   
   DataType_t* p=NULL;
   DataType_t* q=NULL;
   DataType_t* xx;
   double sum = 0.;
   
   size_t i,last;
   size_t step = Slice.stride();
   size_t N = Slice.size();            // number of samples in wavearray
   size_t n = size_t(t*rate()/step);   // # of samples in the window
   
   if(n<4) {
      cout<<"wavearray<DataType_t>::mean() short time window"<<endl;
      return;
   }   

   if(n&1) n--;                        // # of samples in the window - 1

   size_t nM = n/2;                    // index of median sample
   size_t nL = N-nM-1;

   if(pm){
      pm->resize(N/skip);
      pm->start(start());
      pm->rate(rate());
   }

   xx = (DataType_t  *)malloc((n+1)*sizeof(DataType_t));

   p = data+Slice.start();
   q = data+Slice.start();
   for(i=0; i<=n; i++) { 
      xx[i] = *p; 
      sum += xx[i]; 
      p += step;
   }
   last = 0;

   for(i=0; i<N; i++){

      if(pm) {
	    pm->data[i/skip]  = DataType_t(sum/(n+1.));
	 if(clean) q[i*step] -= DataType_t(sum/(n+1.));
      }
      else {
	 if(clean) q[i*step] -= DataType_t(sum/(n+1.));
	 else      q[i*step]  = DataType_t(sum/(n+1.));
      }

      if(i>=nM && i<nL) {              // copy next sample into last
	 sum -= xx[last]; 
	 sum += *p; 
	 xx[last++] = *p; 
	 p += step;
      }

      if(last>n) last = 0;
      
   } 

   free(xx);
   return;
}


template<class DataType_t>
double wavearray<DataType_t>::rms() 
{
   register size_t i;
   register double x = 0.;
   register double y = 0.;
   size_t N = (size()>>2)<<2;
   register DataType_t *p = data + size() - N;

   if(!size()) return 0.;

   for(i=0; i<size()-N; i++) { x += data[i]; y += data[i]*data[i]; }
   for(i=0; i<N; i+=4){ 
      x += p[i] + p[i+1] + p[i+2] + p[i+3];
      y += p[i]*p[i] + p[i+1]*p[i+1] + p[i+2]*p[i+2] + p[i+3]*p[i+3];
   }
   x /= size();
   return sqrt(y/size()-x*x);
}

// running 50% percentile rms
template<class DataType_t>
void wavearray<DataType_t>::rms(double t, 
				wavearray<DataType_t> *pm, 
				bool clean,
				size_t skip)
{
   
   DataType_t*  p=NULL;
   DataType_t*  q=NULL;
   DataType_t** pp;
   DataType_t*  xx;
   DataType_t rm = 1;
   
   size_t i,last;
   size_t step = Slice.stride();
   size_t N = Slice.size();            // number of samples in wavearray
   size_t n = size_t(t*rate()/step);
   
   if(n<4) {
      cout<<"wavearray<DataType_t>::median() short time window"<<endl;
      return;
   }   

   if(n&1) n--;                        // # of samples - 1

   size_t nM = n/2;                    // index of median sample
   size_t nL = N-nM-1;

   if(pm){
      pm->resize(N/skip);
      pm->start(start());
      pm->rate(rate());
   }

   pp = (DataType_t **)malloc((n+1)*sizeof(DataType_t*));
   xx = (DataType_t  *)malloc((n+1)*sizeof(DataType_t));

   p = data+Slice.start();
   q = data+Slice.start();
   for(i=0; i<=n; i++) { 
      xx[i] = *p>0 ? *p : -(*p); 
      pp[i] = xx+i; 
      p += step;
   }
   last = 0;

   for(i=0; i<N; i++){

      if(i==(i/skip)*skip) {
	 waveSplit(pp,0,n,nM);   // median split
	 rm=*pp[nM];
      }

      if(pm) {
	     pm->data[i/skip] = DataType_t(rm/0.6745);
	 if(clean) q[i*step] *= DataType_t(0.6745/rm);
      }
      else {
	 if(clean) q[i*step] *= DataType_t(0.6745/rm);
	 else      q[i*step]  = DataType_t(rm/0.6745);
      }

      if(i>=nM && i<nL) {              // copy next sample into last
	 xx[last++] = *p>0 ? *p : -(*p); 
	 p += step;
      }

      if(last>n) last = 0;
      
   } 

   
   free(pp);
   free(xx);
   return;
}


template<class DataType_t>
double wavearray<DataType_t>::rms(const std::slice &s) 
{
   register size_t i;
   register double a = 0.;
   register double x = 0.;
   register double y = 0.;
   register DataType_t *p = data + s.start();
   size_t n = s.size();
   size_t m = (s.stride()<=0) ? 1 : s.stride();

   if(size()<limit(s)) n = (limit(s) - s.start() - 1)/m;
   if(!n) return 0.;
   size_t N = (n>>2)<<2; 

   for(i=0; i<n-N; i++) 
      a = *p; x += a; y += a*a; p += m;
   for(i=0; i<N; i+=4) { 
      a = *p; x += a; y += a*a; p += m; 
      a = *p; x += a; y += a*a; p += m; 
      a = *p; x += a; y += a*a; p += m; 
      a = *p; x += a; y += a*a; p += m; 
   }
   x /= N;
   return sqrt(y/N-x*x);
}

template<class DataType_t>
DataType_t wavearray<DataType_t>::max() const
{
   if(!size()) return 0;
   register unsigned int i;
   register DataType_t x = data[0];
   for(i=1; i<size(); i++) { if(x<data[i]) x=data[i]; }
   return x;
}

template<class DataType_t>
DataType_t wavearray<DataType_t>::min() const
{
   if(!size()) return 0;
   register size_t i;
   register DataType_t x = data[0];
   for(i=1; i<size(); i++) { if(x>data[i]) x=data[i]; }
   return x;
}


template<class DataType_t>
int wavearray<DataType_t>::getSampleRank(size_t n, size_t l, size_t r) const
{
   DataType_t v;
   register int i = l-1;           // save left boundary 
   register int j = r;             // save right boundary 
   
   v = data[n];                    // pivot
   data[n]=data[r]; data[r]=v;     // store pivot
   
   while(i<j)
   {
      while(data[++i]<v && i<j);
      while(data[--j]>v && i<j);
   }
   data[r]=data[n]; data[n]=v;     // put pivot back
   
   return i-int(l)+1;              // rank ranges from 1 to r-l+1
}

template<class DataType_t>
int wavearray<DataType_t>::getSampleRankE(size_t n, size_t l, size_t r) const
{
   DataType_t v,vv;
   register int i = l-1;           // save left boundary 
   register int j = r;             // save right boundary 
   
   v = data[n];                    // pivot
   data[n]=data[r]; data[r]=v;     // store pivot
   vv = v>0 ? v : -v;              // sort absolute value   

   while(i<j)
   {
      while((data[++i]>0 ? data[i] : -data[i])<vv && i<j);
      while((data[--j]>0 ? data[j] : -data[j])>vv && i<j);
   }
   data[r]=data[n]; data[n]=v;     // put pivot back
   
   return i-int(l)+1;              // rank ranges from 1 to r-l+1
}


template<class DataType_t>
void wavearray<DataType_t>::waveSplit(DataType_t** pp, size_t l, size_t r, size_t m) const
{
   DataType_t v;
   DataType_t* p;
   size_t i = (r+l)/2;
   size_t j = r-1;

// sort l,i,r
   if(*pp[l] > *pp[i]) {p=pp[l]; pp[l]=pp[i]; pp[i]=p;}
   if(*pp[l] > *pp[r]) {p=pp[l]; pp[l]=pp[r]; pp[r]=p;}
   if(*pp[i] > *pp[r]) {p=pp[i]; pp[i]=pp[r]; pp[r]=p;}
   if(r-l < 3) return;                // all sorted
   
   v = *pp[i];                        // pivot
   p=pp[i]; pp[i]=pp[j]; pp[j]=p;     // store pivot
   i = l;
   
   for(;;)
   {
      while(*pp[++i] < v);
      while(*pp[--j] > v);
      if(j<i) break;
      p=pp[i]; pp[i]=pp[j]; pp[j]=p;  // swap i,j
   }
   p=pp[i]; pp[i]=pp[r-1]; pp[r-1]=p; // put pivot  
   
        if(i > m) waveSplit(pp,l,i,m);
   else if(i < m) waveSplit(pp,i,r,m);

   return;
}


template<class DataType_t>
double wavearray<DataType_t>::median(size_t l, size_t r) const
{
   if(!r) r = size()-1;
   if(r<=l) return 0.;

   size_t i;
   size_t N = r-l+1;
   size_t m = N/2+(N&1);  // median

   double x = 0.;

   DataType_t **pp = (DataType_t **)malloc(N*sizeof(DataType_t*));
   for(i=l; i<=r; i++) pp[i] = data + i;
   
   waveSplit(pp,0,N,m);
   x = *pp[m];

   free(pp);
   return x;
}


template<class DataType_t>
void wavearray<DataType_t>::median(double t, 
				   wavearray<DataType_t> *pm, 
				   bool clean,
				   size_t skip)
{
   
   DataType_t*  p=NULL;
   DataType_t*  q=NULL;
   DataType_t** pp;
   DataType_t*  xx;
   DataType_t   am=0;
   
   size_t i,last;
   size_t step = Slice.stride();
   size_t N = Slice.size();            // number of samples in wavearray
   size_t n = size_t(t*rate()/step);   // # of samples in running window
   
   if(n<4) {
      cout<<"wavearray<DataType_t>::median() short time window"<<endl;
      return;
   }   

   if(n&1) n--;                        // # of samples - 1

   size_t nM = n/2;                    // index of median sample
   size_t nL = N-nM-1;

   if(pm){
      pm->resize(N/skip);
      pm->start(start());
      pm->rate(rate()/skip);
   }

   pp = (DataType_t **)malloc((n+1)*sizeof(DataType_t*));
   xx = (DataType_t  *)malloc((n+1)*sizeof(DataType_t));

   p = data+Slice.start();
   q = data+Slice.start();
   for(i=0; i<=n; i++) { 
      xx[i] = *p; 
      pp[i] = xx+i; 
      p += step;
   }
   last = 0;

   for(i=0; i<N; i++){

      if(i==(i/skip)*skip) {
	 waveSplit(pp,0,n,nM);      // median split
	 am=*pp[nM];
      }   

      if(pm) {
	     pm->data[i/skip] = am;
	 if(clean) q[i*step] -= am;
      }
      else {
	 if(clean) q[i*step] -= am;
	 else      q[i*step]  = am;
      }

      if(i>=nM && i<nL) {              // copy next sample into last
	 xx[last++] = *p; 
	 p += step;
      }

      if(last>n) last = 0;
      
   } 

   free(pp);
   free(xx);
   return;
}


template<class DataType_t>
void wavearray<DataType_t>::exponential(double t) 
{
   
   DataType_t*   p=NULL;
   DataType_t*   q=NULL;
   
   size_t i;
   double r;
   size_t last,next;   
   size_t N = Slice.size();            // number of samples in wavearray
   size_t step = Slice.stride();
   size_t n = size_t(t*rate()/step);
   
   if(n<4) {
      cout<<"wavearray<DataType_t>::median() short time window"<<endl;
      return;
   }   

   if(n&1) n--;                        // # of samples in running window

   size_t nM = n/2;                    // index of median sample
   size_t nL = N-nM-1;

   DataType_t** pp = (DataType_t **)malloc((n+1)*sizeof(DataType_t*));
   wavearray<DataType_t> xx(n+1);

   p = data+Slice.start();
   q = data+Slice.start();
   for(i=0; i<=n; i++) { 
      xx.data[i] = *p; 
      pp[i] = xx.data+i;
      p += step;
   }
   last = 0;
   next = 0;

   for(i=0; i<N; i++){

      r = (xx.getSampleRank(next,0,n)-double(nM)-1.)/(nM+1.);
      q[i*step]  = DataType_t(r>0. ? -log(1.-r) : log(1.+r));

      if(i>=nM && i<nL) {              // copy next sample into last
	 xx.data[last++] = *p; p+=step; 
      }
      
      next++;  
      if(next>n) next = 0;
      if(last>n) last = 0;
    
   } 

   free(pp);
   return;
}


template<class DataType_t>
DataType_t wavearray<DataType_t>::rank(double f) const
{
   int i;
   int n = size();
   DataType_t out = 0;
   DataType_t **pp = NULL;
   if(f<0.) f = 0.;
   if(f>1.) f = 1.;

   if(n)
      pp = (DataType_t **)malloc(n*sizeof(DataType_t*));
   else
      return out;

   for(i=0; i<n; i++) pp[i] = data + i;
   qsort(pp, n, sizeof(DataType_t *), 
         (qsort_func)&wavearray<DataType_t>::compare);

   i =int((1.-f)*n);
   if(i==0) out = *(pp[0]);
   else if(i>=n-1) out = *(pp[n-1]);
   else out = (*(pp[i])+*(pp[i+1]))/2;

   for(i=0; i<n; i++) *(pp[i]) = n-i;

   free(pp);
   return out;
}


// symmetric prediction error signal produced with 
// adaptive split lattice algoritm
template<class DataType_t>
void wavearray<DataType_t>::spesla(double T, double w, double oFFset)
{
   int k,j;
   double p,q,xp,xm;

   int K = int(rate()*T+0.5);              // filter duration
   int M = int(rate()*w+0.5);              // integration window
   if(M&1) M++;                            // make M even
   int m = M/2;

   int offset = int(oFFset*rate());        // data offset
   int N  = size();                        // data size
   int nf = offset+m;
   int nl = N-nf;

   wavearray<DataType_t> x0;
   wavearray<DataType_t> x1;

   x0 = *this;                   // previous X iteration
   x1 = *this; 
   x1.add(x0,x0.size()-1,0,1);   // current  X iteration
   x1 *= DataType_t(0.5);                   

//   cout<<nf<<" "<<nl<<" "<<offset<<" "<<K<<" "<<M<<" "<<x0.rms()<<" "<<x1.rms();

   for(k=1; k<K; k++) {

      p = q = 0.;
      for(j=offset; j<offset+M; j++) {  
	 p += x1.data[j]*x1.data[j];
	 q += x0.data[j]*x0.data[j];
      }

      for(j=1; j<N; j++) {
	 if(j>=nf && j<nl) {      // update p & q
	    xp = x1.data[j+m];
	    xm = x1.data[j-m];
	    p += (xp-xm)*(xp+xm);
	    xp = x0.data[j+m];
	    xm = x0.data[j-m];
	    q += (xp-xm)*(xp+xm);
	 }
	 xp = k==1 ? 2*p/q : p/q;
//	 if(fabs(xp)>10.) cout<<" "<<xp;
	 data[j]  = x1.data[j]+x1.data[j-1]-DataType_t(x0.data[j-1]*xp);
      }
      
      x0 = x1;
      x1 = *this;

   }      
//   cout<<" "<<x1.rms()<<endl;

   return;
}


// apply filter
template<class DataType_t>
void wavearray<DataType_t>::lprFilter(wavearray<double>& w)
{
   int i,j;
   int N = size();
   int m = w.size();
   wavearray<DataType_t> x;
   x = *this;

   for(i=0; i<N; i++) {
      for(j=1; j<m && (i-j)>=0; j++) {
	 data[i] += DataType_t(w.data[j]*x.data[i-j]);
      }
   }
   return;
}



// calculate and apply lpr filter
template<class DataType_t>
void wavearray<DataType_t>::lprFilter(double T, int mode, double stride, double oFFset)
{
   int i,j,l,n,nf,nl;
   double duration = size()/rate()-2*oFFset;
   double a=0.;
   double b=1.;

   int N = size();
   int M = int(rate()*T+0.5);              // filter duration
   int k = stride>0 ? int(duration/stride) : 1;   // number of intervals

   if(!k) k++;
   int m = int((N-2.*oFFset*rate())/k);
   if(m&1) m--;                            // make m even

   int offset = int(oFFset*rate()+0.5);    // data offset
   if(offset&1) offset++;

   wavearray<DataType_t> w(m);
   wavearray<DataType_t> x = *this;
   wavearray<double> f;
   w.rate(rate());

   for(l=0; l<k; l++) {

      n = l*m+(N-k*m)/2;
      w.cpf(x,m,n);
      f = w.getLPRFilter(M,0);
      
      nf = l==0   ? 0 : n;
      nl = l==k-1 ? N : n+m;
      n  = offset+M;

      if(mode == 0 || mode == 1) {            // forward LP
	 for(i=nf; i<nl; i++) {            
	    if(i < n) continue; 
	    b = (!mode && i<N-n) ? 0.5 : 1.;  // symmetric LP 
	    for(j=1; j<M; j++) {
	       a = double(f.data[j]*x.data[i-j]);
	       data[i] += DataType_t(a*b);
	    }
	 }
      }

      if(mode == 0 || mode == -1) {           // backward LP
	 for(i=nf; i<nl; i++) {
	    if(i >= N-n) continue; 
	    b = (!mode && i>=n) ? 0.5 : 1.;   // symmetric LP 
	    for(j=1; j<M; j++) {
	       a = double(f.data[j]*x.data[i+j]);
	       data[i] += DataType_t(a*b);  
	    }
	 }
      }
      
   }

   return;
}

//**************************************************************
// calculate autocorrelation function and
// solve symmetric Yule-Walker problem
//**************************************************************
template<class DataType_t>
wavearray<double> wavearray<DataType_t>::getLPRFilter(size_t M, size_t offset)
{
  size_t i,m;
  double f  = 0.03;                     // exclude tail
  size_t N  = size()-2*offset;
  size_t nL = size_t(f*N+0.5);         // left percentile
  size_t nR = N-nL-1;                  // right percentile
  size_t nn = N/2;                     // median
  size_t n  = size()-offset;

  if(size()<=offset) {
     cout<<"wavearray<DataType_t>::getLPRFilter() invalid input parameters\n";
     wavearray<double> a(1);
     return a;
  }

  wavearray<double> r(M);
  wavearray<double> a(M);
  wavearray<double> b(M);

  wavearray<DataType_t> x = *this;

  DataType_t ** pp = (DataType_t **)malloc(N*sizeof(DataType_t*));

  for(i=offset; i<n; i++) pp[i-offset] = x.data + i;

  waveSplit(pp,0,N-1,nn);
  waveSplit(pp,0,nn,nL);
  waveSplit(pp,nn,N-1,nR);

  x -= *pp[nn];                      // subtract median 
  for(i=0;  i<nL; i++) *pp[i] = 0;   // exclude tails
  for(i=nR; i<N;  i++) *pp[i] = 0;   // exclude tails

// autocorrelation

  offset += M;
  n = size()-offset;

  for (m=0; m<M; m++) {
    r.data[m] = 0.;
    for (i=offset; i<n; i++) {
      r.data[m] += x.data[i]*(x.data[i-m] + x.data[i+m])/2.;
    }
//    r.data[m] /= n;
  }


// Levinson algorithm: P.Delsarte & Y.Genin, IEEE, ASSP-35, #5 May 1987

  double p,q;
  double s = r.data[0];

  a = 0.; a.data[0] = 1.; 

  for (m=1; m<M; m++) {

    q = 0.;
    for (i=0; i<m; i++) q -= a.data[i] * r.data[m-i];

    p  = q/s;            // reflection coefficient
    s -= q*p;

    for (i=1; i<=m; i++) b.data[i] = a.data[i] + p*a.data[m-i];
    for (i=1; i<=m; i++) a.data[i] = b.data[i];

  }

/*
  double tmp, num, den;
  M--;

  a.data[1] = - r.data[1] / r.data[0];

  for (m=1; m<M; m++) {

    num = r.data[m + 1];
    den = r.data[0];

    for (i=1; i<=m; i++) {
      num += a.data[i] * r.data[m+1-i];
      den += a.data[i] * r.data[i];
    }

    a.data[m+1] = - num / den;

    for (i=1; i <= ((m+1)>>1); i++) {
      tmp = a.data[i] + a.data[m+1] * a.data[m+1-i];
      a.data[m+1-i] = a.data[m+1-i] + a.data[m+1] * a.data[i];
      a.data[i] = tmp;
    }

  }
  a.data[0] = 1;
*/

  free(pp);
  return a;
}


//**************************************************************
// normalize data by noise variance
// offset  | stride |      
// |*****|*X********X********X********X********X********X*|*****|
//         ^              ^                        ^
//     measurement        |<-      interval      ->|
//**************************************************************
template<class DataType_t>
wavearray<double> 
wavearray<DataType_t>::white(double t, size_t mode, double oFFset, double stride) const
{
   int i,j;
   int N = size();

   double segT = size()/rate();            // segment duration
   if(t<=0.) t = segT-2.*oFFset;
   int  offset = int(oFFset*rate()+0.5);   // offset in samples
   if(offset&1) offset--;                  // make offset even

   if(stride > t || stride<=0.) stride = t;                  
   int K = int((segT-2.*oFFset)/stride);   // number of noise measurement minus 1
   if(!K) K++;                             // special case   

   int n =  N-2*offset;                  
   int k =  n/K;                           // shift in samples
   if(k&1) k--;                            // make k even

   int m = int(t*rate()+0.5);              // number of samples in one interval
   if(m>n || !m) m=n;                      // handle special cases                  
   if(!(m&1)) m++;                         // make m odd

   int mm = m/2;                           // median m
   int mL = int(0.15865*m+0.5);            // -sigma index (CL=0.31732)
   int mR = m-mL-1;                        // +sigma index
   int jL = (N-k*K)/2;                     // array index for first measurement
   int jR = N-offset-m;                    // array index for last measurement
   int jj = jL-mm;                         // start of the first interval

   wavearray<double> meDIan(K+1);
   wavearray<double> norm50(K+1);

   meDIan.rate(rate()/k);
   meDIan.start(start()+jL/rate());
   norm50.rate(rate()/k);
   norm50.start(start()+jL/rate()); 

   if(m<3 || mL<2 || mR>m-2) {
     cout<<"wavearray::white(): too short input array."<<endl;
     return mode!=0 ? norm50 : meDIan;
   }

   register DataType_t *p = data;
   wavearray<DataType_t> w(m);
   double x;
   double r;

   DataType_t ** pp = (DataType_t **)malloc(m*sizeof(DataType_t*));

   for(j=0; j<=K; j++) {

      if(jj < offset)   p  = data + offset;     // left boundary
      else if(jj >= jR) p  = data + jR;         // right boundary
      else  p  = data + jj;
      jj += k;               // update jj

      if(p+m>data+N) cout<<"wavearray::white(): error1\n";

      for(i=0; i<m; i++) pp[i] = p + i;
      waveSplit(pp,0,m-1,mm);
      waveSplit(pp,0,mm,mL);
      waveSplit(pp,mm,m-1,mR);
      meDIan[j] = *pp[mm];
      norm50[j] = (*pp[mR] - *pp[mL])/2.;

   }

   p = data;

   if(mode) {

     mm = jL;
     for(i=0; i<mm; i++){
       x  = double(*p)-meDIan.data[0];
       r  = norm50.data[0];
       *(p++) = mode==1 ? DataType_t(x/r) : DataType_t(x/r/r);
     }
     
     for(j=0; j<K; j++) {
       for(i=0; i<k; i++) {
	 x  = double(*p)-(meDIan.data[j+1]*i + meDIan.data[j]*(k-i))/k;
	 r  = (norm50.data[j+1]*i + norm50.data[j]*(k-i))/k;
	 *(p++) = mode==1 ? DataType_t(x/r) : DataType_t(x/r/r);
       }
     }
     
     mm = (data+N)-p;
     for(i=0; i<mm; i++){
       x  = double(*p)-meDIan.data[K];
       r  = norm50.data[K];
       *(p++) = mode==1 ? DataType_t(x/r) : DataType_t(x/r/r);
     }
     
   }

   free(pp);
   return mode!=0 ? norm50 : meDIan;
}



template<class DataType_t>
void wavearray<DataType_t>::waveSort(DataType_t** pin, size_t l, size_t r) const
{
   size_t k;

   register DataType_t v;
   register DataType_t* p;
   register DataType_t** pp = pin;

   if(pp==NULL) return;      

   register size_t i = (r+l)/2;         // median
   register size_t j = r-1;             // pivot storage index

// sort l,i,r
   if(*pp[l] > *pp[i]) {p=pp[l]; pp[l]=pp[i]; pp[i]=p;}
   if(*pp[l] > *pp[r]) {p=pp[l]; pp[l]=pp[r]; pp[r]=p;}
   if(*pp[i] > *pp[r]) {p=pp[i]; pp[i]=pp[r]; pp[r]=p;}
   if(r-l < 3) return;                  // all sorted
   
   v = *pp[i];                          // pivot
   p=pp[i]; pp[i]=pp[j]; pp[j]=p;       // store pivot
   i = l;
   
   for(;;)
   {
      while(*pp[++i] < v);
      while(*pp[--j] > v);
      if(j<i) break;
      p=pp[i]; pp[i]=pp[j]; pp[j]=p;    // swap i,j
   }
   
   p=pp[i]; pp[i++]=pp[r-1]; pp[r-1]=p; // return pivot back
   
   if(j-l > 2) waveSort(pp,l,j);
   else if(j>l){                        // sort l,k,j
      k = l+1;
      if(*pp[l] > *pp[k]) {p=pp[l]; pp[l]=pp[k]; pp[k]=p;}
      if(*pp[l] > *pp[j]) {p=pp[l]; pp[l]=pp[j]; pp[j]=p;}
      if(*pp[k] > *pp[j]) {p=pp[k]; pp[k]=pp[j]; pp[j]=p;}
   }

   if(r-i > 2) waveSort(pp,i,r);
   else if(r>i){                        // sort i,k,r
      k = i+1;
      if(*pp[i] > *pp[k]) {p=pp[i]; pp[i]=pp[k]; pp[k]=p;} 
      if(*pp[i] > *pp[r]) {p=pp[i]; pp[i]=pp[r]; pp[r]=p;}
      if(*pp[k] > *pp[r]) {p=pp[k]; pp[k]=pp[r]; pp[r]=p;}
   }

   return;
}

//: returns mean and root mean square of the signal.
template<class DataType_t>
double wavearray<DataType_t>::getStatistics(double &m, double &r) const
{
   register size_t i;
   register double a;
   register double b;
   DataType_t *p = const_cast<DataType_t *>(data);
   double y = 0.;
   size_t N = size() - 1 + size_t(size()&1);

   if(!size()) return 0.;

   m = p[0];
   r = p[0]*p[0];
   if(N < size()){
      m += p[N];
      r += p[N]*p[N];
      y += p[N]*p[N-1];
   }

   for(i=1; i<N; i+=2) { 
      a = p[i]; b = p[i+1];
      m += a + b; 
      r += a*a + b*b; 
      y += a*(p[i-1]+b);
   }

   N  = size();
   m = m/double(N);
   r = r/double(N) - m*m;

   y  = 4.*(y/N - m*m + m*(p[0]+p[i]-m)/N);
   y /= 4.*r - 2.*((p[0]-m)*(p[0]-m)+(p[i]-m)*(p[i]-m))/N;
   r = sqrt(r);

   a = (fabs(y) < 1) ? sqrt(0.5*(1.-fabs(y))) : 0.;

   return y>0 ? -a : a;
}


#ifdef _USE_DMT

template<class DataType_t>
wavearray<DataType_t>& wavearray<DataType_t>::
operator=(const TSeries &ts)
{
   Interval ti = ts.getTStep();
   double Tsample = ti.GetSecs();

   unsigned int n=ts.getNSample();
   if(n != size()) resize(n);

   if ( Tsample > 0. )
      rate(double(int(1./Tsample + 0.5)));
   else {
      cout <<" Invalid sampling interval = 0 sec.\n";
   }

   start(double(ts.getStartTime().totalS()));

   TSeries r(ts.getStartTime(), ti, ts.getNSample());
   r = ts;
   float *vec_ref;
   vec_ref= (float*)(r.refData());

   for ( unsigned int i=0; i<n; i++ ) 
      data[i] = (DataType_t) (vec_ref[i]);
   return *this;
}

#endif

// instantiations

#define CLASS_INSTANTIATION(class_) template class wavearray< class_ >;

CLASS_INSTANTIATION(int)
CLASS_INSTANTIATION(short)
CLASS_INSTANTIATION(long)
CLASS_INSTANTIATION(long long)
CLASS_INSTANTIATION(unsigned int)
CLASS_INSTANTIATION(float)
CLASS_INSTANTIATION(double)
//CLASS_INSTANTIATION(std::complex<float>)
//CLASS_INSTANTIATION(std::complex<double>)

#undef CLASS_INSTANTIATION


#if !defined (__SUNPRO_CC)
template wavearray<double>::
wavearray(const double *, unsigned int, double);

template wavearray<double>::
wavearray(const float *, unsigned int, double );

template wavearray<double>::
wavearray(const short *, unsigned int, double );

template wavearray<float>::
wavearray(const double *, unsigned int, double );

template wavearray<float>::
wavearray(const float *, unsigned int, double );

template wavearray<float>::
wavearray(const short *, unsigned int, double );

#else
// FAKE CALLS FOR SUN CC since above instatinations are 
// not recognized
static void fake_instatination_SUN_CC () 
{
   double x;
   float  y;
   short  s;
   wavearray<double> wvdd (&x, 1, 0);
   wavearray<double> wvdf (&y, 1, 0);
   wavearray<double> wvds (&s, 1, 0);
   wavearray<float>  wvfd (&x, 1, 0);
   wavearray<float>  wvff (&y, 1, 0);
   wavearray<float>  wvfs (&s, 1, 0);
}

#endif









