// $Id: wseries.cc,v 1.5 2005/07/01 02:25:58 klimenko Exp $

#define WSeries_CC
#include <time.h>
#include <iostream>
#include <stdexcept>
#include "wseries.hh"
#include "wavearray.hh"

using namespace std;

// constructors

template<class DataType_t>
WSeries<DataType_t>::WSeries() : wavearray<DataType_t>()
{
   this->pWavelet = new WaveDWT<DataType_t>();
   this->pWavelet->allocate(this->size(),this->data);
   this->bpp = 1.;
   this->f_low = 0.;
   this->f_high = 0.;
   this->w_mode = 0;
}

template<class DataType_t>
WSeries<DataType_t>::WSeries(const Wavelet &w) : 
wavearray<DataType_t>()
{ 
   this->pWavelet = NULL;
   this->setWavelet(w);
   this->bpp = 1.;
   this->f_low = 0.;
   this->f_high = 0.;
   this->w_mode = 0;
}  

template<class DataType_t>
WSeries<DataType_t>::
WSeries(const wavearray<DataType_t>& value, const Wavelet &w) : 
wavearray<DataType_t>(value)
{   
   this->pWavelet = NULL;
   this->setWavelet(w);
   this->bpp = 1.;
   this->f_low = 0.;
   this->f_high = value.rate()/2.;
   this->w_mode = 0;
}

template<class DataType_t>
WSeries<DataType_t>::
WSeries(const WSeries<DataType_t>& value) : 
wavearray<DataType_t>(value)
{
   this->pWavelet = NULL;
   this->setWavelet(*(value.pWavelet));
   this->bpp = value.getbpp();
   this->f_low = value.getlow();
   this->f_high = value.gethigh();
   this->w_mode = value.w_mode;
}

// destructor

template<class DataType_t>
WSeries<DataType_t>::~WSeries()
{
   this->pWavelet->release();
   if(this->pWavelet) delete this->pWavelet;
}

// metadata methods

template<class DataType_t>
int WSeries<DataType_t>::
getMaxLevel()
{ 
   int maxlevel = 0;

   if(pWavelet->allocate())
      maxlevel = pWavelet->getMaxLevel();

   return maxlevel;
}

// Accessors

// access to data in wavelet domain
// get wavelet coefficients from a layer with specified frequency index
// if index<0 - get coefficients from the layer = |index|
template<class DataType_t>
int WSeries<DataType_t>::
getLayer(wavearray<DataType_t> &value, int index)
{
   if(index > maxLayer()) index = maxLayer();   
   slice s = pWavelet->getSlice(index);
 
   if(this->limit(s) <= this->size()){
      value.resize(s.size());
      value.rate(this->rate()/s.stride());
      value.start(this->start());
      value.Slice = slice(0,s.size(),1);
      value << (*this)[s];               // get slice of wavelet valarray
      return index;
   }
   else{
      cout<<"WSeries::getLayer(): data length mismatch: "<<this->limit(s)<<" "<<this->size()<<"\n";
      return -1;
   }
}

// access to data in wavelet domain
// put wavelet coefficients into layer with specified frequency index
// if index<0 - put coefficients into layer = |index|
template<class DataType_t>
void WSeries<DataType_t>::
putLayer(wavearray<DataType_t> &value, int index)
{
   slice s = this->pWavelet->getSlice(index);

   if( (s.size() < value.size()) || (this->limit(s) > this->size()) ){
      cout<<"WSeries::putLayer(): invalid array size.\n";
   }      
   else{
      (*this)[s] << value;    // put slice into wavelet valarray
   }
}

// mutators

template<class DataType_t>
void WSeries<DataType_t>::setWavelet(const Wavelet &w)
{ 
   if(pWavelet){              // delete old wavelet object
      pWavelet->release();
      delete pWavelet; 
   } 

   pWavelet = (WaveDWT<DataType_t> *)w.Clone();
   pWavelet->allocate(this->size(), this->data);
}

template<class DataType_t>
void WSeries<DataType_t>::
Forward(int k)
{
   if(pWavelet->allocate()){
      pWavelet->t2w(k);
   }
   else{
      throw std::invalid_argument
      ("WSeries::Forward(): data is not allocated");
   }
}

template<class DataType_t>
void WSeries<DataType_t>::
Forward(wavearray<DataType_t> &x, int k)
{
   wavearray<DataType_t>* p = this;
   if(pWavelet->allocate()) pWavelet->release();
   *p = x;
   f_high = x.rate()/2.;
   pWavelet->allocate(this->size(), this->data);
   pWavelet->reset();
   Forward(k);
}

template<class DataType_t>
void WSeries<DataType_t>::
Forward(wavearray<DataType_t> &x, Wavelet &w, int k)
{
   wavearray<DataType_t>* p = this;
   if(pWavelet->allocate()) pWavelet->release();
   *p = x;
   f_high = x.rate()/2.;
   setWavelet(w);
   Forward(k);
}

template<class DataType_t>
void WSeries<DataType_t>::
Inverse(int k)
{ 
   if(pWavelet->allocate()){
      pWavelet->w2t(k); 
   }
   else{
      throw std::invalid_argument
      ("WSeries::Inverse(): data is not allocated");
   }
}


//: operators =

template<class DataType_t>
WSeries<DataType_t>& WSeries<DataType_t>::
operator=(const wavearray<DataType_t>& a)
{
   wavearray<DataType_t>* p = this;
   if( pWavelet->allocate() ) pWavelet->release();
   if(p->size() != a.size()) pWavelet->reset();
   *p = a;
   f_high = a.rate()/2.;
   pWavelet->allocate(this->size(), this->data);
   return *this;
}

template<class DataType_t>
WSeries<DataType_t>& WSeries<DataType_t>::
operator=(const WSeries<DataType_t>& a)
{
   const wavearray<DataType_t>* p = &a;
   wavearray<DataType_t>* q = this;
   *q = *p;
   setWavelet(*(a.pWavelet));
   bpp = a.getbpp();
   f_low = a.getlow();
   f_high = a.gethigh();
   w_mode = a.w_mode;
   return *this;
}

template<class DataType_t>
WSeries<DataType_t>& WSeries<DataType_t>::
operator[](const std::slice& s)
{
   this->Slice = s;
   if(this->limit() > this->size()){
      cout << "WSeries::operator[]: Illegal argument: "<<this->limit()<<" "<<this->size()<<"\n";
      this->Slice = std::slice(0,this->size(),1);
   }
   return *this;
}

template<class DataType_t>
WSeries<DataType_t>& WSeries<DataType_t>::
operator=(const DataType_t a)
{ this->wavearray<DataType_t>::operator=(a); return *this; }

template<class DataType_t>
WSeries<DataType_t>& WSeries<DataType_t>::
operator*=(const DataType_t a)
{ this->wavearray<DataType_t>::operator*=(a); return *this; }

template<class DataType_t>
WSeries<DataType_t>& WSeries<DataType_t>::
operator-=(const DataType_t a)
{ this->wavearray<DataType_t>::operator-=(a); return *this; }

template<class DataType_t>
WSeries<DataType_t>& WSeries<DataType_t>::
operator+=(const DataType_t a)
{ this->wavearray<DataType_t>::operator+=(a); return *this; }

//template<class DataType_t>
//WSeries<DataType_t>& WSeries<DataType_t>::
//operator*=(wavearray<DataType_t> &a)
//{ this->wavearray<DataType_t>::operator*=(a); return *this; }

template<class DataType_t>
WSeries<DataType_t>& WSeries<DataType_t>::
operator-=(wavearray<DataType_t> &a)
{ this->wavearray<DataType_t>::operator-=(a); return *this; }

template<class DataType_t>
WSeries<DataType_t>& WSeries<DataType_t>::
operator+=(wavearray<DataType_t> &a)
{ this->wavearray<DataType_t>::operator+=(a); return *this; }

template<class DataType_t>
WSeries<DataType_t>& WSeries<DataType_t>::
operator*=(WSeries<DataType_t>& a)
{
   size_t i;
   wavearray<DataType_t> x;
   wavearray<DataType_t>* p  = (wavearray<DataType_t>*)this;
   wavearray<DataType_t>* pa = (wavearray<DataType_t>*)&a;
   size_t max_layer = (maxLayer() > a.maxLayer()) ? a.maxLayer() : maxLayer();

   if(pWavelet->m_TreeType != a.pWavelet->m_TreeType){
     cout<<"WSeries::operator* : wavelet tree type mismatch."<<endl;
     return *this;
   }

   if(this->size()==a.size()) { 
     this->wavearray<DataType_t>::operator*=(*pa); 
     return *this; 
   }

   for(i=0; i<= max_layer; i++)
     (*p)[pWavelet->getSlice(i)] *= (*pa)[a.pWavelet->getSlice(i)];

   return *this;
}

template<class DataType_t>
WSeries<DataType_t>& WSeries<DataType_t>::
operator+=(WSeries<DataType_t>& a)
{
   size_t i;
   wavearray<DataType_t>* p  = (wavearray<DataType_t>*)this;
   wavearray<DataType_t>* pa = (wavearray<DataType_t>*)&a;
   size_t max_layer = (maxLayer() > a.maxLayer()) ? a.maxLayer() : maxLayer();

   if(pWavelet->m_TreeType != a.pWavelet->m_TreeType){
     cout<<"WSeries::operator+ : wavelet tree type mismatch."<<endl;
     return *this;
   }

   if(this->size()==a.size()) { 
     this->wavearray<DataType_t>::operator+=(*pa); 
     return *this; 
   }

   for(i=0; i<= max_layer; i++)
       (*p)[pWavelet->getSlice(i)] += (*pa)[a.pWavelet->getSlice(i)];

   return *this;
}

template<class DataType_t>
WSeries<DataType_t>& WSeries<DataType_t>::
operator-=(WSeries<DataType_t>& a)
{
   size_t i;
   wavearray<DataType_t>* p  = (wavearray<DataType_t>*)this;
   wavearray<DataType_t>* pa = (wavearray<DataType_t>*)&a;
   size_t max_layer = (maxLayer() > a.maxLayer()) ? a.maxLayer() : maxLayer();

   if(pWavelet->m_TreeType != a.pWavelet->m_TreeType){
     cout<<"WSeries::operator- : wavelet tree type mismatch."<<endl;
     return *this;
   }

   if(this->size()==a.size()) { 
     this->wavearray<DataType_t>::operator-=(*pa); 
     return *this; 
   }

   for(i=0; i<= max_layer; i++)
       (*p)[pWavelet->getSlice(i)] -= (*pa)[a.pWavelet->getSlice(i)];

   return *this;
}

template<class DataType_t>
WSeries<DataType_t>& WSeries<DataType_t>::
operator*=(wavearray<DataType_t>& a)
{
   size_t i;
   wavearray<DataType_t>* p  = (wavearray<DataType_t>*)this;
   size_t max_layer = maxLayer()+1;
   
   if(max_layer == a.size()) {
     for(i=0; i< max_layer; i++) {
       (*p)[pWavelet->getSlice(i)] *= a.data[i];
     }     
   }
   else if(this->size()==a.size()) { 
     this->wavearray<DataType_t>::operator*=(a); 
     return *this; 
   }
   else cout<<"WSeries::operator* - no operation is performed"<<endl;

   return *this;
}


//: Dumps data array to file *fname in ASCII format.
template<class DataType_t>
void WSeries<DataType_t>::Dump(const char *fname, int app)
{
  wavearray<DataType_t> a;
  int i,j;
  int n = this->size();
  int m = this->maxLayer()+1;
  char mode[3] = "w";
  if (app == 1) strcpy (mode, "a");

  FILE *fp;

  if ( (fp = fopen(fname, mode)) == NULL ) {
     cout << " Dump() error: cannot open file " << fname <<". \n";
     return;
  };

  if(app == 0) {
    fprintf( fp,"# start time: -start %lf \n", this->Start );
    fprintf( fp,"# sampling rate: -rate %lf \n", this->Rate );
    fprintf( fp,"# number of samples: -size %d \n", (int)this->Size );
    fprintf( fp,"# number of layers: -n %d \n", m );
  }

  for (i = 0; i < m; i++) {
    this->getLayer(a,i);
    n = (int)a.size();
    for(j = 0; j < n; j++) fprintf( fp,"%e ", (float)a.data[j]);
    fprintf( fp,"\n");
  }
  fclose(fp); 
}


template<class DataType_t>
void WSeries<DataType_t>::resize(unsigned int n)
{
   wavearray<DataType_t>* p = this;
   if( pWavelet->allocate() ) pWavelet->release();
   p->wavearray<DataType_t>::resize(n);
   pWavelet->allocate(this->size(), this->data);
   pWavelet->reset();
   bpp = 1.;
   f_low = 0.;
   f_high = this->rate()/2.;
}

template<class DataType_t>
void WSeries<DataType_t>::resample(double f, int nF)
{
   wavearray<DataType_t>* p = this;
   if( pWavelet->allocate() ) pWavelet->release();
   p->wavearray<DataType_t>::resample(f,nF);
   pWavelet->allocate(this->size(), this->data);
   pWavelet->reset();
   bpp = 1.;
   f_low = 0.;
   f_high = 0.;
}

template<class DataType_t>
double WSeries<DataType_t>::coincidence(WSeries<DataType_t>& a, int t, int f, double threshold)
{
#if !defined (__SUNPRO_CC)
   register int i;
   register int j;
   register int u;
   register int v;
   register float* q = NULL;
   register float* p = NULL; 

   int is, ie, js, je;
   
   wavearray<DataType_t> x;
   wavearray<DataType_t> y;

   float energy;

   if(!pWavelet->BinaryTree()) return 1.;

   int ni = maxLayer()+1;
   int nj = this->size()/ni;
   int n = ni-1;
   int m = nj-1;

   bool CROSS = t<0 || f<0;

   t = abs(t);
   f = abs(f);

   float A[ni][nj];
   float B[ni][nj];

   for(i=0; i<=n; i++){
      p  = A[i]; 
      q  = B[i];
      a.getLayer(x,i);	
        getLayer(y,i);	

      for(j=0; j<=m; j++){
	 p[j] = (float)x.data[j];
	 q[j] = (float)y.data[j];
      }
   }

   for(i=0; i<=n; i++){
      p  = A[i]; 
      q  = B[i];
      a.getLayer(x,i);
        getLayer(y,i);
      
      for(j=0; j<=m; j++){

	 if(p[j]==0. && q[j]==0.) continue;

	 is = i-f<0 ? 0 : i-f;
	 js = j-t<0 ? 0 : j-t;
	 ie = i+f>n ? n : i+f;
	 je = j+t>m ? m : j+t;

	 energy = 0.;
	 if(x.data[j]!=0.) {
	    for(u=is; u<=ie; u++)
	       for(v=js; v<=je; v++){
		  if(CROSS && !(i==u || j==v)) continue;
		  if(B[u][v]!=0.) energy += log(fabs(B[u][v])); 
	       }
	    if(energy < threshold) x.data[j]=0;
	 }

	 energy = 0.;
	 if(y.data[j]!=0.) {
	    for(u=is; u<=ie; u++)
	       for(v=js; v<=je; v++){
		  if(CROSS && !(i==u || j==v)) continue;
		  if(A[u][v]!=0.) energy += log(fabs(A[u][v])); 
	       }
	    if(energy < threshold) y.data[j]=0;
	 }
	 
	 if(y.data[j]==0. && x.data[j]!=0.) y.data[j]=DataType_t(a.size());

      }

      putLayer(y,i);

   }
#endif
   return 0.;
}

template<class DataType_t>
double WSeries<DataType_t>::Coincidence(WSeries<DataType_t>& a, double w, double So)
{
   size_t i,k,m;
   size_t event = 0;
   size_t step;
   size_t N = a.size();               
   int n;
   
   double E,S;
   bool pass;

   slice x,y;
   DataType_t* p = NULL;
   DataType_t* q = NULL;
   DataType_t* P = NULL;
   DataType_t* Q = NULL;

   if(pWavelet->m_TreeType != a.pWavelet->m_TreeType){
     cout<<"WSeries::operator- : wavelet tree type mismatch."<<endl;
     return 0.;
   }

   size_t max_layer = (maxLayer() > a.maxLayer()) ? a.maxLayer() : maxLayer();

   for(k=0; k<= max_layer; k++){

      x =   getSlice(k);
      y = a.getSlice(k);

      if(x.size()   != y.size())   continue;
      if(x.stride() != y.stride()) continue;
      if(x.start()  != y.start())  continue;

      step = x.stride();
      n = int(w*a.rate()/2./step);     // 2*(n+1) - coincidence window in pixels 
      if(n < 0) n = 0;                 // not negative n
      if(!n && w>=0.) n++;             // min 0 vs min 3 pixels
      S = log(float(n));               // threshold correction
      S = So + 2.*S/3.;                // threshold
//      S = So + 0.5*log(2*n+1.);      // threshold
//      S = So + S*S/(S+1);            // threshold
      P = a.data+x.start();
      Q = a.data+x.start()+step*(x.size()-1);
      n *= step;

      for(i=x.start(); i<N; i+=step)
      {
	 if(this->data[i] == 0.) continue;

	 p = a.data+i-n;               // start pointer
	 if(p<P) p = P;
	 q = a.data+i+n;               // stop pointer
	 if(q>Q) q = Q;

// calculate total likelihood and number of black pixels
// and set threshold on significance

	 pass = false;
	 E = 0.; m = 0;

	 while(p <= q) {
//	    if(*p>S) { pass = true; break; }
	    if(*p>0) { E += *p; m++; }
	    p += step;
	 }

	 if(m>0 && !pass) {
	    if(gammaCL(E,m) > S-log(double(m))) pass = true;
	 }

	 if(pass) event++;
	 else     this->data[i] = 0;
      }
   }

// handle wavelets with different number of layers

   if(size_t(maxLayer())>max_layer){
      wavearray<DataType_t>* pw  = (wavearray<DataType_t>*)this;
      for(k=max_layer+1; k<= size_t(maxLayer()); k++) {
	(*pw)[getSlice(k)] = 0;
      }
   }

   return double(event)/this->size();
}


template<class DataType_t>
void WSeries<DataType_t>::median(double t, bool r)
{
   int i;
   int M = maxLayer()+1;

   for(i=0; i<M; i++){
      this->setSlice(getSlice(i));
      wavearray<DataType_t>::median(t,NULL,r);
   }

   std::slice S(0,this->size(),1);
   this->setSlice(S);
   return;
}


template<class DataType_t>
void WSeries<DataType_t>::lprFilter(double T, int mode, double stride, double offset)
{
   if(offset<T) offset = T;
   size_t i;
   size_t M = maxLayer()+1;

   wavearray<DataType_t> a;

   for(i=0; i<M; i++){
      getLayer(a,i);
      if(mode<2) a.lprFilter(T,mode,stride,offset);
      else       a.spesla(T,stride,offset);
      putLayer(a,i);
   } 
   return;
}



template<class DataType_t>
WSeries<double> WSeries<DataType_t>::white(double t, size_t mode, double offset, double stride)
{
   int i;
   double segT = this->size()/this->rate();
   if(t <= 0.) t = segT-2.*offset;
   int M = maxLayer()+1;

   if(mode) w_mode = mode;                  
   if(stride > t || stride<=0.) stride = t;                  
   size_t K = size_t((segT-2.*offset)/stride);   // number of noise measurement minus 1
   if(!K) K++;                                   // special case   

   Wavelet* pw = pWavelet->Clone();
   wavearray<DataType_t> a;
   wavearray<double>     b(M*(K+1));
   WSeries<double>       ws(b,*pw);

   for(i=0; i<M; i++){
      getLayer(a,i); 
      b = a.white(t,mode,offset,stride);
      if(b.size() != K+1) cout<<"wseries::white(): array size mismatch\n";
      ws.putLayer(b,i);
      putLayer(a,i);
   } 

   ws.start(b.start());
   ws.rate(b.rate());
   ws.setlow(0.);
   ws.sethigh(this->rate()/2.);

   delete pw;
   return ws;
}


//**************************************************************************
// whiten wavelet series by using TF noise rms array
//**************************************************************************
template<class DataType_t>
bool WSeries<DataType_t>::white(WSeries<double> nRMS)
{
  size_t i,j,k,J;
  int n,m;

  if(!nRMS.size()) return false;

  slice S,N;

  int M     = nRMS.size()-1;             // nRMS boundary
  size_t Io = nRMS.maxLayer()+1;
  size_t I  = this->maxLayer()+1;
  double T  = nRMS.start()-this->start();
  double Ro = nRMS.rate();               // !rate in a single noise layer!
  double Rn = nRMS.gethigh()*2;
  double R  = this->rate();
  double x;

  cout<<Ro<<" "<<Rn<<" "<<T<<endl;

  if(int(Rn+0.1) != int(R+0.1)) {
    cout<<"wseries::white warning: different bandwidth of WSeries and nRMS\n";
    return false;
  }
  
  k = I>Io ? I/Io : Io/I;               
  k = size_t(log(k)/log(2)+0.1);        // number of transformation steps
  
  if(I>Io) this->Inverse(k);
  else     this->Forward(k);
  
  for(i=0; i<I; i++){                   // loop over layers
    S = this->getSlice(i);              // slice of WSeries 
    N = nRMS.getSlice(i);               // slice of noise array
    
    for(j=0; j<S.size(); j++) {
      J = S.start()+j*S.stride();             // index in WSeries
      n = int((J/R-T)*Ro);                    // time index in the noise array
      m = int(N.start())+n*int(N.stride());   // index in nRSM
      x = m<0 || m>M ? 0. : 1./nRMS.data[m];  // noise RMS
      this->data[J] *= x;                     // whiten WSeries 
    }
  }
  
  cout<<this->rms()*this->size()<<endl;

  if(I<Io) this->Inverse(k);
  else     this->Forward(k);
  
  return true;
}

template<class DataType_t>
wavearray<double> WSeries<DataType_t>::filter(size_t n)
{
   size_t i;
   size_t M = maxLayer()+1;
   size_t k = 1<<n;
   double x;

   wavearray<DataType_t> a;
   wavearray<double>     b;
   wavearray<double>     v(M);      // effective variance

   if(!pWavelet->BinaryTree()) { v = 1.; return v; }
   else                        { v = 0.;           }

   Forward(n);                      // n decomposition steps 

   for(i=0; i<M; i++){
      getLayer(a,i); 
      b = a.white(1);
      x = b.data[0];
      v.data[i/k] += x>0. ? 1./x/x : 0.;
      putLayer(a,i);
   } 

   Inverse(n);                      // n reconstruction steps 

   for(i=0; i<v.size(); i++) 
     v.data[i] = sqrt(k/v.data[i]);

   v.start(this->start());

   return v;
}


template<class DataType_t>
WSeries<float> WSeries<DataType_t>::variability(double t, double fl, double fh)
{
   if(fl<0.) fl = this->getlow();
   if(fh<0.) fh = this->gethigh();

   size_t i,j;
   size_t  M = maxLayer()+1;        // number of layers
   size_t  N = this->size()/M;      // zero layer length
   size_t ml = size_t(2.*M*fl/this->rate()+0.5);    // low layer 
   size_t mh = size_t(2.*M*fh/this->rate()+0.5);    // high layer

   if(mh>M) mh = M;

   size_t nL = ml+int((mh-int(ml))/4.+0.5); // left sample (50%)
   size_t nR = mh-int((mh-int(ml))/4.+0.5); // right sample (50%)
   size_t n  = size_t(fabs(t)*this->rate()/M);

   DataType_t* p;
   DataType_t* pp[M];               // sorting
   size_t inDex[M];                 // index(layer)
   size_t laYer[M];                 // layer(index)
   WSeries<float> v;                // effective variance
   WSeries<float> V;                // variance correction
   slice S;
   v.resize(N);

   if(!pWavelet->BinaryTree() || mh<ml+8 || nL<1) 
        { v = 1.; return v; }
   else { v = 0.;           }

// calculate frequency order of layers

   for(j=0; j<M; j++){ 
      S = getSlice(j);
      inDex[j] = S.start();   // index in the array (layer)
      laYer[inDex[j]] = j;    // layer (index in the array)
   }

// calculate variability for each wavelet collumn

   for(i=0; i<N; i++){
     p = this->data + i*M;
     for(j=0; j<M; j++){ pp[j] = &(p[inDex[j]]); }
     this->waveSplit(pp,ml,mh-1,nL-1);          // left split
     this->waveSplit(pp,nL,mh-1,nR);            // right split
     v.data[i]  = float(*pp[nR] - *pp[nL-1])/2./0.6745;
   } 

   v.start(this->start());
   v.rate(this->rate()/M);
   v.setlow(fl);
   v.sethigh(fl);

   if(n<2) return v;

// calculate variability correction

   V=v; V-=1.;
   V.lprFilter(fabs(t),0,0.,fabs(t)+1.);
   v -= V;

// correct variability
   
   p = this->data;

   for(i=0; i<N; i++){             
     for(j=0; j<M; j++){
	if(laYer[j]>=ml && laYer[j]<mh) *p /= (DataType_t)v.data[i];
	p++; 
     }
   }

   return v;
}


template<class DataType_t>
double WSeries<DataType_t>::fraction(double t, double f, int mode)
{
   slice S;
   DataType_t*  p=NULL;
   register DataType_t*  P=NULL;
   register DataType_t** pp;
   DataType_t A, aL, aR;
   size_t i,j,k;
   size_t nS, kS, nL, nR, lS;
   size_t nZero = 0;
   size_t n0 = 1;
   size_t nsub = t>0. ? size_t(this->size()/this->rate()/t+0.1) : 1;
   long r;

   if(!nsub) nsub++;

   f = fabs(f);
   if((f>1. || bpp!=1.) && mode) { 
      cout<<"WSeries fraction(): invalid bpp: "<<bpp<<" fraction="<<f<<endl;
      return bpp;
   }
   if(f>0.) bpp = f;

   size_t M = maxLayer()+1;

   n0 = 1;
   pp = (DataType_t **)malloc(sizeof(DataType_t*));
   wavearray<DataType_t> a(n0);

// percentile fraction

   if(mode && f>0.){  
      for(i=0; i<M; i++){

	  S = getSlice(i);	      
	 nS = S.size()/nsub;                           // # of samles in sub-interval
	 kS = S.stride();                              // stride for this layer
	 lS = nS*nsub<S.size() ? S.size()-nS*nsub : 0; // leftover

// loop over subintervals

	 for(k=0; k<nsub; k++) {

	    p = this->data + nS*k*kS + S.start();  // beginning of subinterval

	    if(k+1 == nsub) nS += lS;     // add leftover to last interval

	    nL = nS&1 ? nS/2 : nS/2-1; 
	    nL = size_t(f*nL);            // set left boundary
	    nR = nS - nL - 1;             // set right boundary
	    if(nL<1 || nR>nS-1) { 
	       cout<<"WSeries::fraction() error: too short wavelet layer"<<endl; 
	       return 0.;
	    }

	    if(nS!=n0) {                      // adjust array length
	       pp = (DataType_t **)realloc(pp,nS*sizeof(DataType_t*));
	       a.resize(nS);       
	       n0 = nS;
	    }

	    for(j=0; j<nS; j++) pp[j] = p + j*kS;

	    this->waveSplit(pp,0,nS-1,nL);           // left split
	    this->waveSplit(pp,nL,nS-1,nR);          // right split
	    aL = *pp[nL]; aR = *pp[nR];

	    for(j=0; j<nS; j++){
	       P =  pp[j]; A = *P;

	            if(j<nL) *P = (DataType_t)fabs(A - aL);
	       else if(j>nR) *P = (DataType_t)fabs(A - aR);
	            else   { *P = 0; nZero++; }
	    
	       if(mode > 1) {                   // do initialization for scrambling
		  a.data[j] = *P;               // save data
		  *P = 0;                       // zero sub-interval
	       }
	    }

	    if(mode == 1) continue;             // all done for mode=1

// scramble	 

	    for(j=0; j<nS; j++){
	       if(a.data[j] == 0.) continue;
	       do{ r = int(nS*drand48()-0.1);}
	       while(p[r*kS] != 0);
	       p[r*kS] = a.data[j]; 
	    }
	 }	
      }
   }
   
   else if(f>0.){                            // random fraction
      M = this->size();
      for(i=0; i<M; i++)
	 if(drand48() > f) { this->data[i] = 0; nZero++; }
   }

   else{                                     // calculate zero coefficients
      M = this->size();
      for(i=0; i<M; i++) {
	 if(this->data[i]==0.) nZero++;
      }      
   }

   free(pp);
   return double(this->size()-nZero)/double(this->size());
}


template<class DataType_t>
double WSeries<DataType_t>::significance(double T, double f)
{
   slice S;
   DataType_t*  p=NULL;
   DataType_t** pp;
   size_t i,j,k,l,m;
   size_t nS,nP,nB;
   size_t M = maxLayer()+1;
   size_t il = size_t(2.*M*getlow()/this->rate());
   size_t ih = size_t(2.*M*gethigh()/this->rate()+0.5);
   int nZero = 0;
   double ratio = double(this->size());

   if(ih>M) ih = M;
   if(il>=ih) { 
      cout<<"WSeries::significance(): invalid low and high:  ";
      cout<<"low = "<<il<<"  high = "<<ih<<endl;
      il = 0;
      ih = M;
   }

// zero unused layers   

   for(i=0; i<M; i++){

      if(i>=il && i<=ih) continue;

      S = getSlice(i);
      k = S.size();
      m = S.stride();
      p = this->data+S.start();
      ratio -= double(k);

      for(j=0; j<k; j++) p[j*m] = 0; 
   }
   ratio /= this->size();            // fraction of pixels between high and low

// calculate number of sub-intervals

   S = getSlice(0);                   // layer 0

   size_t n = size_t(fabs(T)*this->rate()/S.stride()/ratio+0.1); // number of towers
   if(n<1) n = S.size();

   k = S.size()/n;                     // number of sub-intervals
   m = this->size()/S.size();                // number of samples in each tower

   f = fabs(f);
   if(f>1.) f = 1.;
   if(f>0. && f<bpp) bpp = f; 
   nS = n*m;                           // # of samples in one sub-interval
   nB = size_t(bpp*nS*ratio);          // expected number of black pixels

   if(!nS || !nB || this->rate()<=0. || m*S.size()!=this->size()) {
      cout<<"WSeries::significance() error: invalid parameters"<<endl; 
      return 0.;
   } 
   
   l = (S.size()-k*n)*m;               // leftover
   if(l) k++;                          // add one more subinterval

//   cout<<"k="<<k<<" m="<<m<<" bpp="<<bpp<<" nS="<<nS<<endl;

   pp = (DataType_t **)malloc(nS*sizeof(DataType_t*));

   p = this->data;

   for(i=0; i<k; i++){

// fill pp and a

      nP = 0;

      for(j=0; j<nS; j++) {
	 if(*p == 0.) {p++; continue;}
	 *p = (DataType_t)fabs(double(*p));
	 pp[nP++] = p++;
	 nZero++;                                // count non-Zero pixels
      }

      if(nP>2) this->waveSort(pp,0,nP-1);              // sort black pixels

      for(j=0; j<nP; j++) {
	 if(!i && l && pp[j]>=this->data+l) continue;  // handle leftover
	 *pp[j] = nP<nB ? (DataType_t)log(double(nP)/(nP-j)) :  
	                  (DataType_t)log(double(nB)/(nP-j));
	 if(*pp[j] < 0) { 
	    *pp[j] = 0;
	    nZero--;
	 } 
      }

      p = this->data+i*nS+l;
      if(!l) p += nS;
   }
   
   free(pp);
   return double(nZero)/ratio/this->size();
}


template<class DataType_t>
double WSeries<DataType_t>::rsignificance(size_t n, double f)
{
   DataType_t*   p=NULL;
   DataType_t*  px=NULL;
   DataType_t** pp;
   DataType_t** qq;
   DataType_t*  xx;
   DataType_t*  yy;

   double aL, aR;

   size_t i,j,m;
   size_t last, next;
   size_t nS,nB,nL,nR;
   size_t nBlack = 0;
   size_t index;
 
   slice S=getSlice(0);                // layer 0
   size_t N = S.size();                // number of towers in WSeries

   m = this->size()/S.size();                // number of samples in each tower

   f = fabs(f);
   if(f>1.) f = 1.;
   if(f>0. && f<bpp) bpp = f; 
   nS = (2*n+1)*m;                     // # of samples in one sub-interval
   nB = size_t(bpp*nS);                // expected number of black pixels
   if(nB&1) nB++;
   nL = nB/2;                          // left bp boundary
   nR = nS - nL;                       // right bp boundary

   if(!nS || !nB || this->rate()<=0. || m*S.size()!=this->size()) {
      cout<<"WSeries::significance() error: invalid WSeries"<<endl; 
      return 0.;
   } 
   
   pp = (DataType_t **)malloc(nS*sizeof(DataType_t*));
   xx = (DataType_t  *)malloc(nS*sizeof(DataType_t));
   qq = (DataType_t **)malloc(nS*sizeof(DataType_t*));
   yy = (DataType_t  *)malloc(nS*sizeof(DataType_t));

   p = this->data;
   for(j=0; j<nS; j++){ 
      xx[j] = *p; 
      pp[j] = xx+j;
      qq[j] = yy+j;      
      *p++ = 0;
   }
   last = 0;
   next = 0;

   for(i=0; i<N; i++){

      this->waveSplit(pp,0,nS-1,nL-1);       // left split
      this->waveSplit(pp,nL,nS-1,nR);        // right split
      aL = *pp[nL]; aR = *pp[nR];

      for(j=0;  j<nL; j++) yy[j]       = (DataType_t)fabs(*pp[j] - aL);
      for(j=nR; j<nS; j++) yy[j+nL-nR] = (DataType_t)fabs(*pp[j] - aR);

      this->waveSort(qq,0,nB-1);             // sort black pixels

      for(j=0; j<nB; j++) {
	 index = qq[j]-yy;             // index in yy
	 if(index>nL) index+=nR-nL;    // index in pp  
	 index = pp[index]-xx;         // index in xx
	 if(next != index/m) continue; // compare with current tower index in xx
	 this->data[index-next*m+i*m] = (DataType_t)log(double(nB)/(nB-j));
	 nBlack++;
      }

      if(i>=n && i<N-n) {              // copy next tower into last
	 px = xx+last*m;
	 for(j=0; j<m; j++) { *(px++) = *p; *p++ = 0;}
	 last++;                       // update last tower index in array xx
      }

      next++;                          // update current tower index in array xx
      if(next>2*n) next = 0;
      if(last>2*n) last = 0;
      
   } 

   free(pp);
   free(qq);
   free(xx);
   free(yy);

   return double(nBlack)/double(this->size());
}


template<class DataType_t>
double WSeries<DataType_t>::rSignificance(double T, double f, double t)
{
   wavearray<DataType_t> wa;

   DataType_t*   p=NULL;
   DataType_t*  xx;    // buffer for data
   DataType_t*  yy;    // buffer for black pixels
   DataType_t** px;    // pointers to xx
   DataType_t** py;    // pointers to yy
   DataType_t** pp;    // pointers to this->data, so *pp[j] == xx[j]

   double aL, aR;

   size_t i,j,m,l,J;
   size_t last, next;
   size_t nS,nB,nL,nR;
   size_t nBlack = 0;
   size_t index;

   size_t M = maxLayer()+1;            // number of samples in a tower  
   size_t il = size_t(2.*M*getlow()/this->rate());
   size_t ih = size_t(2.*M*gethigh()/this->rate()+0.5);

   if(ih>M) ih = M;
   if(il>=ih) { 
      cout<<"WSeries::significance(): invalid low and high:  ";
      cout<<"low = "<<il<<"  high = "<<ih<<endl;
      il = 0;
      ih = M;
   }

   m = ih-il;                          // # of analysis samples in a tower
 
   for(j=0; j<il; j++){ this->getLayer(wa,j); wa=1234567891.; this->putLayer(wa,j); }
   for(j=ih; j<M; j++){ this->getLayer(wa,j); wa=1234567891.; this->putLayer(wa,j); }

   t *= double(M)/m; 
   T *= double(M)/m; 

   slice S=getSlice(0);                // layer 0
   size_t N = S.size();                // number of towers in WSeries
   size_t k = size_t(t*this->rate()/M);      // sliding step in towers
   size_t n = size_t(T*this->rate()/M/2.);   // 1/2 sliding window in towers

   if(t<=0. || k<1) k = 1;
   if(T<=0. || n<1) n = 1;

   size_t Nnk = N-n-k;
   size_t nW  = (2*n+k)*M;             // total # of samples in the window

   f = fabs(f);
   if(f>1.) f = 1.;
   if(f>0. && f<bpp) bpp = f; 
   nS = (2*n+k)*m;                     // # of analysis samples in the window
   nB = size_t(bpp*nS);                // expected number of black pixels
   if(nB&1) nB++;
   nL = nB/2;                          // left bp boundary
   nR = nS - nL;                       // right bp boundary

   if(!nS || !nB || this->rate()<=0. || M*S.size()!=this->size()) {
      cout<<"WSeries::significance() error: invalid WSeries"<<endl; 
      return 0.;
   } 
   
   pp = (DataType_t **)malloc(nS*sizeof(DataType_t*));
   px = (DataType_t **)malloc(nS*sizeof(DataType_t*));
   py = (DataType_t **)malloc(nS*sizeof(DataType_t*));
   xx = (DataType_t  *)malloc(nS*sizeof(DataType_t));
   yy = (DataType_t  *)malloc(nS*sizeof(DataType_t));

   p = this->data;
   J = 0;
   for(i=0; i<nW; i++){ 
      if(*p != 1234567891.){
	 xx[J] = *p; 
	 pp[J] =  p; 
	 px[J] = xx+J;
	 py[J] = yy+J;
	 J++;
      }      
      *p++ = 0;
   }
   last = 0;
   next = 0;

   if(J != nS) {
     cout<<"wseries::rSignificance() error 1 - illegal sample count"<<endl;
     exit(0);
   }

   for(i=0; i<N; i+=k){

      this->waveSplit(px,0,nS-1,nL-1);       // left split
      this->waveSplit(px,nL,nS-1,nR);        // right split
      aL = *px[nL]; aR = *px[nR];

      for(j=0;  j<nL; j++) yy[j]       = (DataType_t)fabs(*px[j] - aL);
      for(j=nR; j<nS; j++) yy[j+nL-nR] = (DataType_t)fabs(*px[j] - aR);

      if(nB != nS-nR+nL) {
	cout<<"wseries::rSignificance:  nB="<<nB<<",  N="<<nS-nR+nL<<endl;
	nB = nS-nR+nL;
      }

      this->waveSort(py,0,nB-1);             // sort black pixels

      for(j=0; j<nB; j++) {            // save rank in *this
	 index = py[j]-yy;             // index in yy
	 if(index>nL) index+=nR-nL;    // index in xx  
	 index = px[index]-xx;         // index in pp
         index = pp[index]-this->data; // index in WS data array

         if(index>=i*M && index<(i+k)*M) {  // update pixels in window t
           if(this->data[index]!=0) {
             cout<<"WSeries::rSignificance error: "<<this->data[index]<<endl;
           }
           this->data[index] = DataType_t(log(double(nB)/(nB-j))); // save rank
           nBlack++;
         }
      }

      for(l=i; l<i+k; l++){             // copy towers
	 if(l>=n && l<Nnk) {            // copy next tower into last
	    J = last*m;                 // pointer to last tower storage at xx
	    for(j=0; j<M; j++) { 
	       if(*p != 1234567891.) { xx[J] = *p; pp[J++]=p; }
	       *p++ = 0;
	    }
	    last++;                    // update last tower index in array xx
	    if(J != last*m) { 
	       cout<<"wseries::rSignificance() error 2 - illegal sample count"<<endl;
	       exit(0);
	    }
	 }

	 next++;                       // update current tower index in array xx
	 if(next==2*n+k) next = 0;
	 if(last==2*n+k) last = 0;
      }
   } 

   free(pp);
   free(px);
   free(py);
   free(xx);
   free(yy);

   return double(nBlack)/double(this->size());
}



template<class DataType_t>
double WSeries<DataType_t>::gSignificance(double T, double f, double t)
{
   wavearray<DataType_t> wa;

   DataType_t*   p=NULL;
   DataType_t*  xx;    // buffer for data
   DataType_t*  yy;    // buffer for black pixels
   DataType_t** px;    // pointers to xx
   DataType_t** py;    // pointers to yy
   DataType_t** pp;    // pointers to this->data, so *pp[j] == xx[j]

   double aL, aR;

   size_t i,j,m,l,J;
   size_t last, next;
   size_t nS,nB,nL,nR;
   size_t nBlack = 0;
   size_t index;

   size_t M = maxLayer()+1;            // number of samples in a tower  
   size_t il = size_t(2.*M*getlow()/this->rate());
   size_t ih = size_t(2.*M*gethigh()/this->rate()+0.5);

   if(ih>M) ih = M;
   if(il>=ih) { 
      cout<<"WSeries::significance(): invalid low and high:  ";
      cout<<"low = "<<il<<"  high = "<<ih<<endl;
      il = 0;
      ih = M;
   }

   m = ih-il;                          // # of analysis samples in a tower
 
   for(j=0; j<il; j++){ this->getLayer(wa,j); wa=1234567891.; this->putLayer(wa,j); }
   for(j=ih; j<M; j++){ this->getLayer(wa,j); wa=1234567891.; this->putLayer(wa,j); }

   t *= double(M)/m; 
   T *= double(M)/m; 

   slice S=getSlice(0);                // layer 0
   size_t N = S.size();                // number of towers in WSeries
   size_t k = size_t(t*this->rate()/M);      // sliding step in towers
   size_t n = size_t(T*this->rate()/M/2.);   // 1/2 sliding window in towers

   if(t<=0. || k<1) k = 1;
   if(T<=0. || n<1) n = 1;

   size_t Nnk = N-n-k;
   size_t nW  = (2*n+k)*M;             // total # of samples in the window

   f = fabs(f);
   if(f>1.) f = 1.;
   if(f>0. && f<bpp) bpp = f; 
   nS = (2*n+k)*m;                     // # of analysis samples in the window
   nB = size_t(bpp*nS);                // expected number of black pixels
   if(nB&1) nB++;
   nL = nB/2;                          // left bp boundary
   nR = nS - nL;                       // right bp boundary

   if(!nS || !nB || this->rate()<=0. || M*S.size()!=this->size()) {
      cout<<"WSeries::gSignificance() error: invalid WSeries"<<endl; 
      return 0.;
   } 
   
   pp = (DataType_t **)malloc(nS*sizeof(DataType_t*));
   px = (DataType_t **)malloc(nS*sizeof(DataType_t*));
   py = (DataType_t **)malloc(nS*sizeof(DataType_t*));
   xx = (DataType_t  *)malloc(nS*sizeof(DataType_t));
   yy = (DataType_t  *)malloc(nS*sizeof(DataType_t));

   p = this->data;
   J = 0;
   for(i=0; i<nW; i++){ 
      if(*p != 1234567891.){
	 xx[J] = *p; 
	 pp[J] =  p; 
	 px[J] = xx+J;
	 py[J] = yy+J;
	 J++;
      }      
      *p++ = 0;
   }
   last = 0;
   next = 0;

   if(J != nS) {
     cout<<"wseries::gSignificance() error 1 - illegal sample count"<<endl;
     exit(0);
   }

   for(i=0; i<N; i+=k){

      this->waveSplit(px,0,nS-1,nL-1);       // left split
      this->waveSplit(px,nL,nS-1,nR);        // right split
      aL = *px[nL]; aR = *px[nR];

      for(j=0;  j<nL; j++) yy[j]       = (DataType_t)fabs(*px[j]);
      for(j=nR; j<nS; j++) yy[j+nL-nR] = (DataType_t)fabs(*px[j]);

      if(nB != nS-nR+nL) {
	cout<<"wseries::gSignificance:  nB="<<nB<<",  N="<<nS-nR+nL<<endl;
	nB = nS-nR+nL;
      }

      this->waveSort(py,0,nB-1);             // sort black pixels

      for(j=0; j<nB; j++) {            // save significance in *this
	 index = py[j]-yy;             // index in yy
	 if(index>nL) index+=nR-nL;    // index in pp  
	 index = px[index]-xx;         // index in xx
	 *(pp[index]) = pow(*py[j]+1.11/2,2)/2./1.07 + log(bpp); // save significance
      }
      nBlack += nB;

      for(l=i; l<i+k; l++){             // copy towers
	 if(l>=n && l<Nnk) {            // copy next tower into last
	    J = last*m;                 // pointer to last tower storage at xx
	    for(j=0; j<M; j++) { 
	       if(*p != 1234567891.) { xx[J] = *p; pp[J++]=p; }
	       *p++ = 0;
	    }
	    last++;                    // update last tower index in array xx
	    if(J != last*m) { 
	       cout<<"wseries::gSignificance() error 2 - illegal sample count"<<endl;
	       exit(0);
	    }
	 }

	 next++;                       // update current tower index in array xx
	 if(next==2*n+k) next = 0;
	 if(last==2*n+k) last = 0;
      }
   } 

   free(pp);
   free(px);
   free(py);
   free(xx);
   free(yy);

   return double(nBlack)/double(this->size());
}



template<class DataType_t>
double WSeries<DataType_t>::pixclean(double S)
{
   size_t k;
   size_t event = 0;
   int i, j, n;
   int nm, np, mp, mm;
   
   bool one;

   wavearray<DataType_t>  a;
   wavearray<DataType_t>  am;
   wavearray<DataType_t>  ac;
   wavearray<DataType_t>  ap;
   wavearray<DataType_t>* p;
   wavearray<DataType_t>* pm;
   wavearray<DataType_t>* pc;
   wavearray<DataType_t>* pp;

   size_t max_layer = maxLayer()+1;

   pc = &ac; pp = &ap; pm = p = NULL;
   mp = mm = 1;
   getLayer(a,0);
   ac = a;

   for(k=1; k<=max_layer; k++){

      if(k<max_layer) getLayer(*pp,k);  // next layer
      else pp = NULL;
      
      if(pp!=NULL) mp = pp->size()/pc->size();  // scale for upper layer
      if(pm!=NULL) mm = pc->size()/pm->size();  // scale for lower layer

      n  = pc->size()-1; 

      for(i=0; i<=n; i++) {
	 one = true;

	 if(pc->data[i] == 0.)        continue;
	 if(pc->data[i] > 9.7) cout<<"pixclean: "<<pc->data[i]<<endl;
	 event++;
	 if(i>0 && pc->data[i-1]!=0.) continue;
	 if(i<n && pc->data[i+1]!=0.) continue;

// work on upper (+) layer

	 if(pp!=NULL) {
	    nm = i*mp-1;              // left index for + layer
	    np = i*mp+2;              // right index for + layer
	    if(nm < 0) nm = 0; 
	    if(np > n) np = n;
	    
	    for(j=nm; j<np; j++) {
	       if(pp->data[j] != 0) {
		  one = false;
		  break;
	       }
	    } 
	 }
	 if(!one) continue;

// work on lower (-) layer

	 if(pm!=NULL) {
	    nm = i/mm-1;                // left index for + layer
	    np = i/mm+2;              // right index for + layer
	    if(nm < 0) nm = 0; 
	    if(np > n) np = n;
	    
	    for(j=nm; j<np; j++) {
	       if(pm->data[j] != 0) {
		  one = false;
		  break;
	       }
	    } 
	 }
	 if(!one) continue;

	 if(pc->data[i]<S) {a.data[i]=0; event--;}
      }

      putLayer(a,k-1);

// shaffle layers

      if(pp==NULL) break; 

       a = *pp;
       p = pm==NULL ? &am : pm;
      pm = pc;
      pc = pp;
      pp = p;
   }

   return double(event)/this->size();
}


template<class DataType_t>
double WSeries<DataType_t>::percentile(double f, int mode, WSeries<DataType_t>* pin)
{
   slice S;
   DataType_t*  p=NULL;
   register DataType_t*  P=NULL;
   register DataType_t** pp;
   double A, aL, aR;
   double x;
   size_t i,j;
   size_t nS, kS, mS, nL, nR;
   size_t nZero = 0;
   long r;

   f = fabs(f);
   if((f>=1. || bpp!=1.) && mode) { 
      cout<<"WSeries percentile(): invalid bpp: "<<bpp<<" fraction="<<f<<endl;
      return bpp;
   }
   bpp = f;

   if(pin) *this = *pin;                    // copy input wavelet if specified

   size_t M = maxLayer()+1;
   WaveDWT<DataType_t>* pw = pWavelet;

   S=pw->getSlice(0);	      
   size_t n0 = S.size();
   if(n0) pp = (DataType_t **)malloc(n0*sizeof(DataType_t*));
   else return 0.;

   wavearray<DataType_t> a(n0);
   wavearray<DataType_t> b;

   if(mode && f>0.){                                 // percentile fraction
      for(i=0; i<M; i++){

	  S = pw->getSlice(i);	      
	 nS = S.size();
	 kS = S.stride();
	 mS = S.start();
	  p = this->data+S.start();
	 nL = size_t(f*nS/2.+0.5);
	 nR = nS - nL;

	 if(nL<2 || nR>nS-2) { 
	   cout<<"WSeries::percentile() error: too short wavelet layer"<<endl; 
	   return 0.;
	 }
	 
	 if(nS!=n0) {
	    pp = (DataType_t **)realloc(pp,nS*sizeof(DataType_t*));
	    a.resize(nS);       
	 }

	 for(j=0; j<nS; j++) pp[j] = p + j*kS;

	 this->waveSplit(pp,0,nS-1,nL-1);            // left split
	 this->waveSplit(pp,nL,nS-1,nR);             // right split
	 aL = double(*pp[nL-1]); aR = double(*pp[nR]);

	 for(j=0; j<nS; j++){
	    P =  pp[j]; A = double(*P);

//   	         if(j<nL) *P = -sqrt(A*A - aL*aL);
//	    else if(j>nR) *P =  sqrt(A*A - aR*aR);

   	         if(j<nL) *P = DataType_t(fabs(A - aL));
	    else if(j>nR) *P = DataType_t(fabs(A - aR));
		 else   { *P = 0; nZero++; }
	    
	    if(mode == -1) continue;            // all done for mode = -1
	    if(pin) pin->data[mS+P-p] = *P;     // update pin slice
	    if(j>nL && j<nR) continue;          // skip zero amplitudes
	    a.data[(P-p)/kS] = *P;              // save data
	    if(j<nL) *P *= -1;                  // absolute value
	    if(j>=nR) pp[nL+j-nR] = P;          // copy right sample
	 }

	 if(mode == -1) continue;               // keep wavelet amplitudes

	 nL *= 2;
	 this->waveSort(pp,0,nL-1);                   // sort black pixels

	 if(abs(mode)!=1) b = a;

	 for(j=0; j<nL; j++){
	    r = (pp[j]-p)/kS;
//	    x = a.data[r]<0. ? -double(nL)/(nL-j) : double(nL)/(nL-j);
	    x = log(double(nL)/(nL-j));
	    *pp[j] = mode==1 ? DataType_t(x) : 0;
	    if(mode>1) a.data[r]=DataType_t(x);  // save data for random sample
	 }

	 if(abs(mode)==1) continue;

// scramble	 

	 for(j=0; j<nL; j++){
	    P = pp[j];
	    do{ r = int(nS*drand48()-0.1);}
	    while(p[r*kS] != 0);
	    p[r*kS] = a.data[(P-p)/kS]; 
	    if(pin) pin->data[mS+r*kS] = b.data[(P-p)/kS]; 
	 }
      }	
   }
   
   else if(f>0.){                            // random fraction
      M = this->size();
      for(i=0; i<M; i++)
	 if(drand48() > f) { this->data[i] = 0; nZero++; }
   }

   else{                                     // calculate zero coefficients
      M = this->size();
      for(i=0; i<M; i++)
	 if(this->data[i]==0) nZero++;
   }

   free(pp);
   return double(this->size()-nZero)/double(this->size());
}


template<class DataType_t>
WSeries<double> WSeries<DataType_t>::calibrate(size_t n, double df,
					       d_complex* R, d_complex* C,
					       wavearray<double> &a,
					       wavearray<double> &g,
					       size_t channel)
{
   size_t i,k,m,N;
   size_t count, step;
   size_t M = maxLayer()+1;
  
   double left, righ;
   double tleft, trigh;
   double reH, imH;
   double c, t, dt;
   double cstep = 1./a.rate();
   double sTart = this->start();
   double sTop  = this->start()+this->size()/this->rate();
   DataType_t* p;
   slice S;

   wavecomplex* pR=R;
   wavecomplex* pC=C;

   Wavelet* pw = pWavelet->Clone();
   wavearray<double> alp;
   wavearray<double> gam;
   wavearray<double> reR(M);
   wavearray<double> reC(M);
   wavearray<double> imR(M);
   wavearray<double> imC(M);

   alp = a; alp.start(0.);
   gam = a; gam.start(0.);

// select alpha to be in WB segment
   count=0;
   for(i=0; i<a.size(); i++) {
      if(a.start()+i/a.rate() < sTart) continue;  // skip samples in front
      if(a.start()+i/a.rate() > sTop)  break;     // skip samples in the back
      if(alp.start()==0.) alp.start(a.start()+i/a.rate());
      alp.data[count++] = a.data[i];
   }
   alp.resize(count);

// select gamma to be in WB segment
   count=0;
   for(i=0; i<g.size(); i++) {
      if(g.start()+i/g.rate() < sTart) continue;  // skip samples in front
      if(g.start()+i/g.rate() > sTop)  break;     // skip samples in the back
      if(gam.start()==0.) gam.start(g.start()+i/g.rate());
      gam.data[count++] = g.data[i];
   }
   gam.resize(count);

   if(gam.size()>alp.size()) gam.resize(alp.size());
   if(gam.size()<alp.size()) alp.resize(gam.size());

   wavearray<double>     x(M*alp.size());
   WSeries<double>       cal(x,*pw);

   if(!alp.size() || a.rate()!=g.rate()) { 
      cout<<"WSeries<DataType_t>::calibrate() no calibration data\n";
      return cal;
   }

   cal = 0.;
   left = righ = 0.;

   reR = 0.; reC = 0.;
   imR = 0.; imC = 0.;
   for(k=0; k<M; k++){

      S = getSlice(k);
      left  = righ;                       // left border
      righ += this->rate()/2./S.stride();       // right border
      if(righ > n*df) break;

// average R and C

      count = 0;
      while(left+count*df < righ){
	 reR.data[k] += pR->real();
	 imR.data[k] += pR->imag();
	 reC.data[k] += pC->real();
	 imC.data[k] += pC->imag();
	 count++; pR++; pC++;
      }
      reR.data[k] /= count; reC.data[k] /= count;
      imR.data[k] /= count; imC.data[k] /= count;

// calculate calibration constants

      cal.getLayer(x,k);
      for(i=0; i<alp.size(); i++){

	 if(alp.data[i]<=0. || gam.data[i]<=0.) {
	    cout<<"WSeries<DataType_t>::calibrate() zero alpha error\n";
	    alp.data[i] = 1.;
	    gam.data[i] = 1.;
	 }

	 reH = reR.data[k]*reC.data[k]-imR.data[k]*imC.data[k];
	 imH = reR.data[k]*imC.data[k]+imR.data[k]*reC.data[k];
	 reH = 1.+(reH-1.)*gam.data[i];
	 imH*= gam.data[i];
	 x.data[i]  = sqrt(reH*reH+imH*imH);
	 x.data[i] /= sqrt(reC.data[k]*reC.data[k]+imC.data[k]*imC.data[k]);
	 x.data[i] /= channel==0 ? alp.data[i] : gam.data[i];
      }
      cal.putLayer(x,k);

// apply energy calibration

      S    = getSlice(k);
      step = S.stride();
      N    = S.size();
      p    = this->data + S.start();
      dt   = step/this->rate();        // sampling time interval
      t    = this->start();            // time stamp
      sTart = alp.start();       // first calibration point
      sTop  = alp.start()+(alp.size()-1)*cstep;         // last calibration sample

      tleft = sTart;             // left and right borders defining beginning and end
      trigh = sTart+cstep;       // of the current calibration cycle.
      m = 0;

      for(i=0; i<N; i++){
	 t += dt;
	 if(t <= sTart) *p *= (DataType_t)x.data[0];
	 else if(t >= sTop) *p *= (DataType_t)x.data[alp.size()-1];
	 else {
	    if(t>trigh) { tleft=trigh; trigh+=cstep; m++; }
	    c = (t-tleft)/cstep;
	    *p *= DataType_t(x.data[m]*(1-c) + x.data[m+1]*c);
	 }	    
	 p += step;
      }

   }

   return cal;
}


// instantiations

#define CLASS_INSTANTIATION(class_) \
template class WSeries< class_ >;

CLASS_INSTANTIATION(float)
CLASS_INSTANTIATION(double)

#undef CLASS_INSTANTIATION











