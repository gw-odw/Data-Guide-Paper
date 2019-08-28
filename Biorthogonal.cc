// Wavelet Analysis Tool
//--------------------------------------------------------------------
// Implementation of 
// Bi-othogonal wavelet transforms using lifting scheme 
// References:
//   A.Cohen, I.Daubechies, J.Feauveau Bases of compactly supported wavelets
//   Comm. Pure. Appl. Math. 45, 485-560, 1992
//   W. Sweldens - Building your own wavelets at home
//--------------------------------------------------------------------
//$Id: Biorthogonal.hh,v 0.2 2001/08/06 19:37:00 klimenko Exp $

#define BIORTHOGONAL_CC

#include "Biorthogonal.hh"

//namespace datacondAPI {
//namespace wat {

// constructors

template<class DataType_t> Biorthogonal<DataType_t>::
Biorthogonal(const Wavelet &w) : 
WaveDWT<DataType_t>(w) 
{ 
   setFilter();
}

template<class DataType_t> Biorthogonal<DataType_t>::
Biorthogonal(const Biorthogonal<DataType_t> &w) : 
WaveDWT<DataType_t>(w) 
{ 
   setFilter();
}

template<class DataType_t> Biorthogonal<DataType_t>::
Biorthogonal(int m, int tree, enum BORDER border) :
WaveDWT<DataType_t>(m,m,tree,border) 
{
   setFilter();
}

// destructor
template<class DataType_t>
Biorthogonal<DataType_t>::~Biorthogonal()
{ 
   if(PForward) delete [] PForward;
   if(PInverse) delete [] PInverse;
   if(UForward) delete [] UForward;
   if(UInverse) delete [] UInverse;
}

// clone
template<class DataType_t>
Biorthogonal<DataType_t>* Biorthogonal<DataType_t>::Clone() const
{
  return new Biorthogonal<DataType_t>(*this);
}

// set filter and wavelet type
template<class DataType_t>
void Biorthogonal<DataType_t>::setFilter()
{ 
   int n = this->m_H;

   n = (n>>1)<<1;
   if(n < 2) n=4;
//   if(n > 30) n=30;   // limit is due to the unrolled code length

   PForward=new double[n];
   PInverse=new double[n];
   UForward=new double[n];
   UInverse=new double[n];

   for(int i=0; i<n; i++) 
   {
      PForward[i] = Lagrange(n,i,0.);
      UForward[i] = 0.5*PForward[i];
      PInverse[i] = -PForward[i];
      UInverse[i] = -UForward[i];
//      printf("%8.3e ",PForward[i]);
   }
   this->m_H = n;
   this->m_L = n;
   this->m_WaveType = BIORTHOGONAL;
}

// forward function does one step of forward transformation.
// <level> input parameter is the level to be transformed
// <layer> input parameter is the layer to be transformed.
template<class DataType_t>
void Biorthogonal<DataType_t>::forward(int level,int layer)
{
   int i;
   int step = 1<<level;
   this->predict(level, layer, PForward);
   this->update(level, layer, UForward);
   const int nS = this->nWWS>>level;                       // number of samples in the layer
   DataType_t *pData = this->pWWS+this->getOffset(level,layer);  // pointer to the first sample in the layer 
   if(this->m_Heterodine){
      for(i=1; i<nS; i+=4) {pData[i*step] *= -1;}  // heterodine detailes coefficients
   }
}

// inverse function does one step of inverse transformation.
// <level> input parameter is the level to be reconstructed
// <layer> input parameter is the layer to be reconstructed.
template<class DataType_t>
void Biorthogonal<DataType_t>::inverse(int level,int layer)
{
   int i;
   int step = 1<<level;
   const int nS = this->nWWS>>level;                       // number of samples in the layer
   DataType_t *pData = this->pWWS+this->getOffset(level,layer);  // pointer to the first sample in the layer 

   if(this->m_Heterodine){
      for(i=1; i<nS; i+=4) {pData[i*step] *= -1;}    // heterodine detailes coefficients
   }

   this->update(level, layer, UInverse);
   this->predict(level, layer, PInverse);
}


// instantiations

#define CLASS_INSTANTIATION(class_) template class Biorthogonal< class_ >;

CLASS_INSTANTIATION(float)
CLASS_INSTANTIATION(double)
//CLASS_INSTANTIATION(std::complex<float>)
//CLASS_INSTANTIATION(std::complex<double>)

#undef CLASS_INSTANTIATION

//template Biorthogonal<float>::
//Biorthogonal(const Biorthogonal<float> &);
//template Biorthogonal<double>::
//Biorthogonal(const Biorthogonal<double> &);

//}  // end namespace wat
//}  // end namespace datacondAPI









