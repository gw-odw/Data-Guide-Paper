// $Id: WaveDWT.cc,v 1.6 2002/06/27 06:30:43 klimenko Exp $

#define WAVEDWT_CC

#include <sstream>
#include <stdexcept>
#include "WaveDWT.hh"

//namespace datacondAPI {
//namespace wat {

using namespace std;

// constructors

template<class DataType_t>
WaveDWT<DataType_t>::
WaveDWT(int mH,int mL,int tree, enum BORDER border) :
   Wavelet(mH, mL, tree, border), pWWS(NULL), nWWS(0)
{ }

template<class DataType_t>
WaveDWT<DataType_t>::WaveDWT(const Wavelet &w) :
   Wavelet(w), pWWS(NULL), nWWS(0)
{ }

template<class DataType_t>
WaveDWT<DataType_t>::WaveDWT(const WaveDWT<DataType_t> &w) :
   Wavelet(w), pWWS(NULL), nWWS(0)
{ }

// destructor
template<class DataType_t>
WaveDWT<DataType_t>::~WaveDWT()
{
   release();
}

template<class DataType_t>
WaveDWT<DataType_t>* WaveDWT<DataType_t>::Clone() const
{
  return new WaveDWT<DataType_t>(*this);
}

// - Virtual procedure for forward wavelet transform.
// - Real code for decompose will appear in a derivative
//   classe, since it depends of the wavelet type.
// - Only ldeep steps of decomposition will be done.
// - By default ldeep=1, which means do one step.
template<class DataType_t>
void WaveDWT<DataType_t>::t2w(int ldeep)
{
   int maxLevel = getMaxLevel();

   int levs = m_Level;
   int levf = m_Level+ldeep;
   if((ldeep == -1) || (levf > maxLevel)) levf = maxLevel;

   int layf;
   for(int level=levs; level<levf; level++)
   {
      layf = (m_TreeType) ? 1<<level : 1;

      for(int layer=0; layer<layf; layer++)
         forward(level,layer);

      m_Level=level+1;

   }

   m_Level=levf;
}

// - Virtual procedure for inverse wavelet transform.
// - Real code for reconstruct will appear in a derivative
//   classe, since it depends of the type of wavelet.
// - Only ldeep steps of reconstruction will be done.
// - By default ldeep=1, which means do one step.
template<class DataType_t>
void WaveDWT<DataType_t>::w2t(int ldeep)
{
   int levs = m_Level;
   int levf = m_Level-ldeep;
   if((ldeep == -1) || (levf < 0)) levf = 0;

   int layf;
   for(int level=levs-1; level>=levf; level--)
   {
      layf = (m_TreeType) ? 1<<level : 1;

      for(int layer=0; layer<layf; layer++)
         inverse(level,layer);

      m_Level=level;

   }
   m_Level=levf;
}


// access function
// convert layer index or frequency index to slice
template<class DataType_t>
slice WaveDWT<DataType_t>::getSlice(const int index)
{
   int level = m_Level;
   int layer = abs(index);

   int maxlayer = (BinaryTree()) ? (1<<level)-1 : level;

   if(layer > maxlayer){
      layer = maxlayer;

	  stringstream oss;
	  oss << "WaveDWT::getSlice(): argument "<<index<<" is set to " 
		<< layer << endl;

      std::invalid_argument exception(oss.str());
      throw exception;
   }

   if(BinaryTree()){
      layer = convertF2L(level,layer);
   }
   else {

      if(layer) {               // detail layers
         level -= layer-1;
         layer = 1;
      }
      else{                     // approximation layer
         layer = 0;
      }
   }

   return getSlice(level,layer);
}

// convert (level,layer) to slice
template<class DataType_t>
slice WaveDWT<DataType_t>::getSlice(const int level, const int layer)
{
   if(!allocate()){
      std::invalid_argument("WaveDWT::getSlice(): data is not allocated");
      return slice(0,1,1);
   }

   size_t m = nWWS>>level;                   // number of elements
   size_t s = 1<<level;                      // slice step
   size_t i = getOffset(level,layer);        // first element

   if(i+(m-1)*s+1 > nWWS){
      std::invalid_argument("WaveDWT::getSlice(): invalide arguments");
      return slice(0,1,1);
   }

   return slice(i,m,s);
}


// Calculate maximal level of wavelet decomposition,
// which depends on the length of wavelet filters.
template<class DataType_t>
int WaveDWT<DataType_t>::getMaxLevel()
{
   if(!allocate()) return 0;
   return Wavelet::getMaxLevel((int)nWWS);
}

// allocate input data
template<class DataType_t>
bool WaveDWT<DataType_t>::
allocate(size_t n, DataType_t *p)
{
   bool allocate = false;
   if(pWWS == NULL && n>0 && p != NULL){
      allocate = true;
      pWWS = p;
      nWWS = n;
   }
   return allocate;
}

// check allocation status
template<class DataType_t>
bool WaveDWT<DataType_t>::allocate()
{
   return (pWWS == NULL || nWWS == 0) ? false : true;
}

// release input data
template<class DataType_t>
void WaveDWT<DataType_t>::release()
{
   pWWS = NULL;
   nWWS = 0;
}

// forward does one step of forward Fast Wavelet Transform.
// It's implemented for FWT with even number of coefficients.
// Also the lenght of high and low pass filters is the same,
// It is used for Daubechies, Symlet and Meyer wavelets.
//
// <level> input parameter is the level to be transformed
// <layer> input parameter is the layer to be transformed.
// <pF>    wavelet filter, the pF length is m_H=m_L
//
// note: borders are handled in B_CYCLE mode
// note: even wavelet mode is standard

template<class DataType_t>
void WaveDWT<DataType_t>::forwardFWT(int level, int layer,
                                  const double* pLPF,
                                  const double* pHPF)
{
   int VM = (m_H/2&1);           // 1 if Odd Vanishing Moments
   int nS = nWWS>>level;         // number of samples in the layer
   int kL = -(m_H/2-VM);         // k left limit for left border
   int iL = nS-m_H+1;            // i left limit for regular case

// switch to odd wavelet mode
//   if(m_Parity && layer&1) kL -= VM ? 1 : -1;
   if(m_Parity) kL -= VM ? 1 : -1;

   int iR = nS+kL;               // i right limit for right border

//EVM--------------odd wavelet mode-----------------------
//
// LP       [a0]       [h0] | h1  h2  h3   0   0   0   0  h0 |        [s0]
// HP       [d0]       [h3] |-h2  h1 -h0   0   0   0   0  h3 |        [s1]
//          [a1]            |  0  h0  h1  h2  h3   0   0   0 |        [s2]
//   for    [d1] =          |  0  h3 -h2  h1 -h0   0   0   0 |     X  [s3]
//   DB2    [a2] =          |  0   0   0  h0  h1  h2  h3   0 |     X  [s4]
//          [d2]            |  0   0   0  h3 -h2  h1 -h0   0 |        [s5]
//          [a3]            |  0   0   0   0   0  h0  h1  h2 | [ h3]  [s6]
//          [d3]            |-h0   0   0   0   0  h3 -h2  h1 | [-h0]  [s7]
//
//EVM--------------even border handling---------------------
//
// LP       [a0]   [h0  h1] | h2  h3   0   0   0   0  h0  h1 |   [s0]
// HP       [d0]   [h3 -h2] | h1 -h0   0   0   0   0  h3 -h2 |   [s1]
//          [a1]            | h0  h1  h2  h3   0   0   0   0 |   [s2]
//   for    [d1] =          | h3 -h2  h1 -h0   0   0   0   0 | X [s3]
//   DB2    [a2] =          |  0   0  h0  h1  h2  h3   0   0 | X [s4]
//          [d2]            |  0   0  h3 -h2  h1 -h0   0   0 |   [s5]
//          [a3]            |  0   0   0   0  h0  h1  h2  h3 |   [s6]
//          [d3]            |  0   0   0   0  h3 -h2  h1 -h0 |   [s7]
//
//
// temp array:  a d a d a d ..... a d a d a d
// a - approximations, d - details
//---------------------------------------------------------


//OVM---------------odd border handling-----------------------
//
// index  i:  -3  -2  -1 |   0  1   2  3   4  5   6  7   8
//   limits:                                 iL     iR
//
// index  k:  -3  -2  -1 |   0  1   2  3   4  5   6  7   8  9  10 11 |  12 13
//   limits:  kL                                                           kR
//
//   LP       h0  h1  h2 |  h3 h4  h5  0   0  0   0  0   0 h0  h1 h2 |
//   HP       h5 -h4  h3 | -h2 h1 -h0  0   0  0   0  0   0 h5 -h4 h3 |
//                    h0 |  h1 h2  h3 h4  h5  0   0  0   0  0   0 h0 |
//   for              h5 | -h4 h3 -h2 h1 -h0  0   0  0   0  0   0 h5 |
//   DB3                 |   0 h0  h1 h2  h3 h4  h5  0   0  0   0  0 |
//                       |   0 h5 -h4 h3 -h2 h1 -h0  0   0  0   0  0 |
//                       |   0  0   0 h0  h1 h2  h3 h4  h5  0   0  0 |
//                       |   0  0   0 h5 -h4 h3 -h2 h1 -h0  0   0  0 |
//                       |   0  0   0  0   0 h0  h1 h2  h3 h4  h5  0 |
//                       |   0  0   0  0   0 h5 -h4 h3 -h2 h1 -h0  0 |
//                       |  h5  0   0  0   0  0   0 h0  h1 h2  h3 h4 |  h5
//                       | -h0  0   0  0   0  0   0 h5 -h4 h3 -h2 h1 | -h0
//
//OVM---------------even border handling-----------------------
//
// index  i:      -2  -1 |  0   1  2   3  4   5  6   7  8   9 10
//   limits:                                    iL            iR
//
// index  k:      -2  -1 |  0   1  2   3  4   5  6   7  8   9 10  11 | 12  13 14
//   limits:      kL                                                          kR
//
//   LP           h0  h1 | h2  h3 h4  h5  0   0  0   0  0   0 h0  h1 |
//   HP           h5 -h4 | h3 -h2 h1 -h0  0   0  0   0  0   0 h5 -h4 |
//                       | h0  h1 h2  h3 h4  h5  0   0  0   0  0   0 |
//   for                 | h5 -h4 h3 -h2 h1 -h0  0   0  0   0  0   0 |
//   DB3                 |  0   0 h0  h1 h2  h3 h4  h5  0   0  0   0 |
//                       |  0   0 h5 -h4 h3 -h2 h1 -h0  0   0  0   0 |
//                       |  0   0  0   0 h0  h1 h2  h3 h4  h5  0   0 |
//                       |  0   0  0   0 h5 -h4 h3 -h2 h1 -h0  0   0 |
//                       |  0   0  0   0  0   0 h0  h1 h2  h3 h4  h5 |
//                       |  0   0  0   0  0   0 h5 -h4 h3 -h2 h1 -h0 |
//                       | h4  h5  0   0  0   0  0   0 h0  h1 h2  h3 | h4  h5
//                       | h1 -h0  0   0  0   0  0   0 h5 -h4 h3 -h2 | h1 -h0
//
// temp array:  a d a d a d ..... a d a d a d
// a - approximations, d - details
//---------------------------------------------------------

   if(pLPF==NULL || pHPF==NULL) return;

   register int i,j,k;
   register double sumA, sumD, data;
   register const double *p = pLPF;
   register const double *q = pHPF;
   register DataType_t *pD;
   register int stride = 1<<level; // stride parameter

// pointer to the first sample in the layer
   DataType_t *pData = pWWS+getOffset(level,layer);

   double *temp = new double[nS]; // temporary array

// left border

   i = kL;

   while(i<0) {
      sumA=0.; sumD=0.;

      for(j=0; j<m_H; j++) {
         k = i+j;
         if(k < 0) k += nS;
         data = pData[k<<level];
         sumA += *p++ * data;
         sumD += *q++ * data;
      }

      *temp++ = sumA;
      *temp++ = sumD;
      i += 2;
      p -= m_H;
      q -= m_H;
   }

// processing data in the middle of array

   while(i<iL) {
      pD = pData + (i<<level) - stride;
      sumA=0.; sumD=0.;

      for(j=0; j<m_H; j+=2) {
         data = *(pD += stride);
         sumA += *p++ * data;
         sumD += *q++ * data;
         data = *(pD += stride);
         sumA += *p++ * data;
         sumD += *q++ * data;
      }

      *temp++ = sumA;
      *temp++ = sumD;
      i += 2;
      p -= m_H;
      q -= m_H;
   }

// right border

   while(i<iR) {
      sumA=0.; sumD=0.;

      for(j=0; j<m_H; j++) {
         k = i+j;
         if(k >= nS) k -= nS;
         data = pData[k<<level];
         sumA += *p++ * data;
         sumD += *q++ * data;
      }

      *temp++ = sumA;
      *temp++ = sumD;
      i += 2;
      p -= m_H;
      q -= m_H;
   }

// writing data back from temporary storage
   for(i=nS-1; i>=0; i--) {pData[i<<level] = *(--temp);}

   if(m_Heterodine){
      for(i=1; i<nS; i+=4) {pData[i<<level] *= -1;}     // heterodine detailes coefficients
   }

   delete [] temp;
}


// inverse does one step of inverse Fast Wavelet Transform.
// It's implemented for FWT with even number of coefficients.
// Also the lenght of high and low pass filters is the same
// It is used for Daubechies and Symlet wavelets.
//
// <level> input parameter is the level to be transformed
// <layer> input parameter is the layer to be transformed.
// <pF>    wavelet filter, the pF length is m_H=m_L
//
// note: borders are handled in B_CYCLE mode
// note: even wavelet mode is standard

template<class DataType_t>
void WaveDWT<DataType_t>::inverseFWT(int level, int layer,
                                  const double* pLPF,
                                  const double* pHPF)
{
   if(pLPF==NULL || pHPF==NULL) return;

   int VM = (m_H/2&1);           // 1 if Odd Vanishing Moments
   int nS = nWWS>>level;         // number of samples in the layer
   int kL = -(m_H/2-2+VM);       // k left limit
   int iL = nS-m_H+1;            // i left limit

   //   bool ODD = m_Parity && layer&1;
   bool ODD = m_Parity;
   if(ODD) { kL -= 2-2*VM; }

   int iR = nS+kL;               // i right limit

// iLP filter for db2:  h3 -h0  h1 -h2
// iHP filter for db2:  h2  h1  h0  h3
// inverse matrix is transpose of forward matrix
//EVM------------------odd border handling-----------------------
//
// index  i:       -2  -1    0   1   2   3   4   5   6
//   limits:                                iL      iR
//
// index  k:       -2  -1    0   1   2   3   4   5   6   7   8   9
//   limits:       kL
// iHP      [s0]   h3 -h0 | h1 -h2   0   0   0   0  h3 -h0 |          [a0]
// iLP      [s1]          | h2  h1  h0  h3   0   0   0   0 |          [d0]
//          [s2]          | h3 -h0  h1 -h2   0   0   0   0 |          [a1]
//   for    [s3] =        |  0   0  h2  h1  h0  h3   0   0 |        X [d1]
//   DB2    [s4] =        |  0   0  h3 -h0  h1 -h2   0   0 |        X [a2]
//          [s5]          |  0   0   0   0  h2  h1  h0  h3 |          [d2]
//          [s6]          |  0   0   0   0  h3 -h0  h1 -h2 |          [a3]
//          [s7]          | h0  h3   0   0   0   0  h2  h1 | h0  h3   [d3]
//
//EVM---------------even border handling-----------------------
//
// iLP      [s0]          | h2  h1  h0  h3   0   0   0   0 |   [a0]
// iHP      [s1]          | h3 -h0  h1 -h2   0   0   0   0 |   [d0]
//          [s2]          |  0   0  h2  h1  h0  h3   0   0 |   [a1]
//   for    [s3] =        |  0   0  h3 -h0  h1 -h2   0   0 | X [d1]
//   DB2    [s4] =        |  0   0   0   0  h2  h1  h0  h3 | X [a2]
//          [s5]          |  0   0   0   0  h3 -h0  h1 -h2 |   [d2]
//          [s6]          | h0  h3   0   0   0   0  h2  h1 |   [a3]
//          [s7]          | h1 -h2   0   0   0   0  h3 -h0 |   [d3]
//
// temp array:  a d a d a d ..... a d a d a d
// a - approximations, d - details
//---------------------------------------------------------

// iLP filter for db3:  h4  h1  h2  h3  h0  h5
// iHP filter for db3:  h5 -h0  h3 -h2  h1 -h4
//OVM---------------odd border handling-----------------------
//
// index  i:      -2  -1 |  0   1  2   3  4   5  6   7  8   9 10
//   limits:                                    iL            iR
//
// index  k:      -2  -1 |  0   1  2   3  4   5  6   7  8   9 10  11 |  12 13
//   limits:      kL                                                           kR
//
//   HP           h5 -h0 | h3 -h2 h1 -h4  0   0  0   0  0   0 h5 -h4 |
//   LP                  | h4  h1 h2  h3 h0  h5  0   0  0   0  0   0 |
//                       | h5 -h0 h3 -h2 h1 -h4  0   0  0   0  0   0 |
//   for                 |  0   0 h4  h1 h2  h3 h0  h5  0   0  0   0 |
//   DB3                 |  0   0 h5 -h0 h3 -h2 h1 -h4  0   0  0   0 |
//                       |  0   0  0   0 h4  h1 h2  h3 h0  h5  0   0 |
//                       |  0   0  0   0 h5 -h0 h3 -h2 h1 -h4  0   0 |
//                       |  0   0  0   0  0   0 h4  h1 h2  h3 h0  h5 |
//                       |  0   0  0   0  0   0 h5 -h0 h3 -h2 h1 -h4 |
//                       | h0  h5  0   0  0   0  0   0 h4  h1 h2  h3 |  h0  h5
//                       | h1 -h4  0   0  0   0  0   0 h5 -h0 h3 -h2 |  h1 -h4
//                       | h2  h3 h0  h5  0   0  0   0  0   0 h4  h1 |  h2  h3 h0 h5
//
//OVM---------------even border handling-----------------------
//
// index  i:      -2  -1 |  0   1  2   3  4   5  6   7  8
//   limits:                                    iL     iR
//
// index  k:      -2  -1 |  0   1  2   3  4   5  6   7  8   9 10  11 | 12  13 14
//   limits:      kL                                                          kR
//
//   LP           h4  h1 | h2  h3 h0  h5  0   0  0   0  0   0 h4  h1 |
//   HP           h5 -h0 | h3 -h2 h1 -h4  0   0  0   0  0   0 h5 -h0 |
//                       | h4  h1 h2  h3 h0  h5  0   0  0   0  0   0 |
//   for                 | h5 -h0 h3 -h2 h1 -h4  0   0  0   0  0   0 |
//   DB3                 |  0   0 h4  h1 h2  h3 h0  h5  0   0  0   0 |
//                       |  0   0 h5 -h0 h3 -h2 h1 -h4  0   0  0   0 |
//                       |  0   0  0   0 h4  h1 h2  h3 h0  h5  0   0 |
//                       |  0   0  0   0 h5 -h0 h3 -h2 h1 -h4  0   0 |
//                       |  0   0  0   0  0   0 h4  h1 h2  h3 h0  h5 |
//                       |  0   0  0   0  0   0 h5 -h0 h3 -h2 h1 -h4 |
//                       | h0  h5  0   0  0   0  0   0 h4  h1 h2  h3 | h0  h5
//                       | h1 -h4  0   0  0   0  0   0 h5 -h0 h3 -h2 | h1 -h4
//
// temp array:  a d a d a d ..... a d a d a d
// a - approximations, d - details
//---------------------------------------------------------


   register long i,j,k;
   register double sumA, sumD, data;
   register const double *p = pLPF;
   register const double *q = pHPF;
   register DataType_t *pD;
   register int stride = 1<<level; // stride parameter

// pointer to the first sample in the layer
   DataType_t *pData = pWWS+getOffset(level,layer);

   double *temp = new double[nS]; // temporary array

   if(m_Heterodine){
      for(i=1; i<nS; i+=4) {pData[i<<level] *= -1;}     // heterodine detailes coefficients
   }

// left border

   i = kL;

   if(ODD){               // handle ODD wavelet mode on left border
      p = pHPF;
      *temp = 0;

      for(j=0; j<m_H; j++) {
         k = i+j;
         if(k < 0) k += nS;
         *temp += *p++ * pData[k<<level];
      }
      temp++;
      i += 2;
      p = pLPF;
   }

   while(i<0) {
      sumA = 0.; sumD = 0.;

      for(j=0; j<m_H; j++) {
         k = i+j;
         if(k < 0) k += nS;
         data = pData[k<<level];
         sumA += *p++ * data;
         sumD += *q++ * data;
      }

      *temp++ = sumA;
      *temp++ = sumD;
      i += 2;
      p -= m_H;
      q -= m_H;
   }

// processing data in the middle of array

   while(i<iL) {
      pD = pData + (i<<level) - stride;
      sumA = 0.; sumD = 0.;

      for(j=0; j<m_H; j+=2) {
         data = *(pD += stride);
         sumA += *p++ * data;
         sumD += *q++ * data;
         data = *(pD += stride);
         sumD += *q++ * data;
         sumA += *p++ * data;
      }

      *temp++ = sumA;
      *temp++ = sumD;
      i += 2;
      p -= m_H;
      q -= m_H;
   }

// right border

   while(i<iR) {
      sumA = 0.; sumD = 0.;

      for(j=0; j<m_H; j++) {
         k = i+j;
         if(k >= nS) k -= nS;
         data = pData[k<<level];
         sumA += *p++ * data;
         sumD += *q++ * data;
      }

      *temp++ = sumA;
      *temp++ = sumD;
      i += 2;
      p -= m_H;
      q -= m_H;
   }

   if(ODD){         // handle odd wavelet mode on right border
      q = pLPF;
      *temp = 0.;

      for(j=0; j<m_H; j++) {
         k = i+j;
         if(k >= nS) k -= nS;
         *temp += *q++ * pData[k<<level];
      }
      temp++;
   }

// writing data back from temporary storage
   for(i=nS-1; i>=0; i--)
      pData[i<<level] = *(--temp);

   delete [] temp;
}

// predict function does one lifting prediction step
// <level> input parameter is the level to be transformed
// <layer> input parameter is the layer to be transformed.
// <p_H>   pointer to prediction filter. It has length m_H
template<class DataType_t>
void WaveDWT<DataType_t>::predict(int level,
                                  int layer,
                                  const double* p_H)
{
   level++;                          // increment level (next level now)

//------------------predict border handling-----------------------
// use even samples to predict odd samples
// an example for m_H=8 and 20 samples
// i index limits     nL............nM....nR-1
// i index         :  -3 -2 -1 0 1 2 3 4 5 6
// j index (approx):  -3 -2 -1 0 1 2 3 4 5 6 7 8 9 | 10 11 12 13
// odd samples:                 0 1 2 3 4 5 6 7 8 9
//                              L L L M M M R R R R
// L,R - samples affected by borders
//   M - not affected samples
//---------------------------------------------------------

   int nS = nWWS>>level;         // nS - number of samples in the layer
   int nL = 1-m_H/2;             // nL - left limit of predict i index
   int nR = nS+nL;               // nR - right limit of predict i index
   int nM = nS-m_H+1;            // nM - number of M samples (first aR sample)
   int mM = nM<<level;           // (number of M samples)<<level

   double data;
   double hsum = 0.;             // filter sum

   register int i,j,k;
   register double sum = 0.;

   register const double *h;             // pointer to filter coefficient
   register const DataType_t *dataL;     // pointer to left data sample
   register const DataType_t *dataR;     // pointer to right data sample
   register const int stride = 1<<level; // stride parameter

   double *pBorder=new double[2*(m_H-nL)]; // border array
   double *pB;

   DataType_t *dataA, *dataD;

   dataA=pWWS+getOffset(level,layer<<1);     // pointer to approximation layer
   dataD=pWWS+getOffset(level,(layer<<1)+1); // pointer to detail layer

   for(k=0; k<m_H; k++) hsum += p_H[k];

// left border

   pB = pBorder;
   dataL = dataA;                            // first (left) sample
   for(k=0; k<(m_H-nL); k++){
      j = k + nL;
      pB[k] = *(dataL + abs(j<<level));

      if(j>=0) continue;

      data = *(dataL + ((j+nS)<<level));
      switch (m_Border) {
         case B_PAD_ZERO : pB[k] = 0.;     break;  // pad zero
         case B_PAD_EDGE : pB[k] = *dataL; break;  // pad by edge value
         case B_CYCLE    : pB[k] = data;   break;  // cycle data
              default    :                 break;  // mirror or interpolate
      }
   }

   for(i=nL; i<0; i++) {                    // i index

      if(m_Border != B_POLYNOM){
         sum = 0.;
         for(k=0; k<m_H/2; k++)
            sum += p_H[k] * (pB[k] + pB[m_H-1-k]);
         pB++;
      }
      else{
//       pB = pBorder - nL;                                 // point to dataA
//       sum = hsum*Nevill(i+0.5-nL, m_H+i, pB, pB+m_H);    // POLYNOM1
         pB = pBorder - nL;                                 // point to dataA
         sum = hsum*Nevill(i+0.5-nL, m_H+2*i, pB, pB+m_H);  // POLYNOM2
      }

      *dataD -= sum;
       dataD += stride;
   }

// regular case (no borders)

  k = (m_H-1)<<level;

  for(i=0; i<mM; i+=stride) {
    dataL = dataA+i;
    dataR = dataL+k;
    h = p_H;

    sum=0.;
    do sum  += *(h++) * (*dataL + *dataR);
    while((dataL+=stride) < (dataR-=stride));
/*
    sum  = *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto P0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto P0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto P0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto P0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto P0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto P0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto P0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto P0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto P0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto P0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto P0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto P0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto P0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto P0;
    sum  += *(h++) * (*dataL + *dataR);

P0: *dataD -= sum; */
    *dataD -= sum;
     dataD += stride;
  }

// right border

   pB = pBorder;
   dataR = dataA + ((nS-1)<<level);     // last (right) sample

   for(k=0; k<(m_H-nL+1); k++){
      j = m_H - 1 - k;
      pB[k] = *(dataR - abs(j<<level));

      if(j>=0) continue;

      data = *(dataR - ((j+nS)<<level));
      switch (m_Border) {
         case B_PAD_ZERO : pB[k] = 0.;     break;  // pad zero
         case B_PAD_EDGE : pB[k] = *dataR; break;  // pad by edge value
         case B_CYCLE    : pB[k] = data;   break;  // cycle data
              default    :                 break;  // mirror or interpolate
      }
   }

   k = 0;
   for(i=nM; i<nR; i++) {

      if(m_Border != B_POLYNOM){
         sum = 0.; pB++;

         for(k=0; k<m_H/2; k++)
            sum += p_H[k] * (pB[k] + pB[m_H-1-k]);
      }
      else{
//       sum = hsum*Nevill(0.5-nL, nS-i, ++pB, pBorder+m_H+1);
         k += 2;
         sum = hsum*Nevill((m_H-k-1)/2., m_H-k, pB+k, pBorder+m_H+1);
         if(k == m_H) sum = *(pB+k-1) * hsum;
      }

      *dataD -= sum;
       dataD += stride;
   }

  delete [] pBorder;
}


// update function does one lifting update step
// <level> input parameter is the level to be transformed.
// <layer> input parameter is the layer to be transformed.

template<class DataType_t>
void WaveDWT<DataType_t>::update(int level,
                                 int layer,
                                 const double* p_L)
{
   level++;                      // current level

//------------------update border handling-----------------------
// use odd samples to update even samples
// an example for m_H=8 and 20 samples
//                               L L L L M M M R R R
// even samples  :               0 1 2 3 4 5 6 7 8 9
// j index (detail):  -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 | 10 11 12
// i index         :  -4 -3 -2 -1 0 1 2 3 4 5
// i index limits     nL...............nM..nR-1
// L,R - samples affected by borders
//   M - not affected samples
//---------------------------------------------------------

   int nS = nWWS>>level;         // nS - number of samples in the layer
   int nL = -m_L/2;              // nL - left limit of update i index
   int nR = nS+nL;               // nR - right limit of update i index
   int nM = nS-m_L+1;            // nM - number of M samples
   int mM = nM<<level;           // (number of M samples)<<level

   double data;
   double hsum = 0.;             // filter sum

   register int i,j,k;
   register double sum = 0.;

   register const double *h;             // pointer to filter coefficient
   register const DataType_t *dataL;     // pointer to left data sample
   register const DataType_t *dataR;     // pointer to right data sample
   register const int stride = 1<<level; // stride parameter

   DataType_t *dataA, *dataD;

   double *pBorder=new double[2*(m_L-nL)]; // border array
   double *pB;

   dataA=pWWS+getOffset(level,layer<<1);     // pointer to approximation layer
   dataD=pWWS+getOffset(level,(layer<<1)+1); // pointer to detail layer

   for(k=0; k<m_L; k++) hsum += p_L[k];

// left border

   pB = pBorder;
   dataL = dataD;                            // first (left) sample
   for(k=0; k<(m_L-nL); k++){
      j = k + nL;
      pB[k] = *(dataL + abs(j<<level));

      if(j>=0) continue;

      data = *(dataL + ((j+nS)<<level));
      switch (m_Border) {
         case B_PAD_ZERO : pB[k] = 0.;     break;  // pad zero
         case B_PAD_EDGE : pB[k] = *dataL; break;  // pad by edge value
         case B_CYCLE    : pB[k] = data;   break;  // cycle data
              default    :                 break;  // mirror or interpolate
      }
   }

   for(i=nL; i<0; i++) {                    // i index

      if(m_Border != B_POLYNOM){
         sum = 0.;
         for(k=0; k<m_L/2; k++)
            sum += p_L[k] * (pB[k] + pB[m_L-1-k]);
         pB++;
      }
      else{
//       pB = pBorder - nL;                                 // point to dataD
//       sum = hsum*Nevill(i-0.5-nL, m_L+i, pB, pB+m_L);    // POLYNOM1
         pB = pBorder - nL;                                 // point to dataD
         sum = hsum*Nevill(i-0.5-nL, m_L+2*i, pB, pB+m_L);  // POLYNOM2
      }

      *dataA += sum;
       dataA += stride;
   }

// regular case (no borders)

  k = (m_L-1)<<level;

  for(i=0; i<mM; i+=stride) {
    dataL = dataD+i;
    dataR = dataL+k;
    h = p_L;

    sum=0.;
    do sum  += *(h++) * (*dataL + *dataR);
    while((dataL+=stride) < (dataR-=stride));
/*
    sum  = *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto U0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto U0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto U0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto U0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto U0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto U0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto U0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto U0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto U0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto U0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto U0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto U0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto U0;
    sum  += *(h++) * (*dataL + *dataR);
    if ((dataL+=stride) > (dataR-=stride)) goto U0;
    sum  += *(h++) * (*dataL + *dataR);

U0: *dataA += sum; */
    *dataA += sum;
     dataA += stride;
  }

// right border

   pB = pBorder;
   dataR = dataD + ((nS-1)<<level);     // last detail sample

   for(k=0; k<(m_L-nL); k++){
      j = m_L - 1 - k;
      pB[k] = *(dataR - abs(j<<level));

      if(j>=0) continue;

      data = *(dataR - ((j+nS)<<level));
      switch (m_Border) {
         case B_PAD_ZERO : pB[k] = 0.;     break;  // pad zero
         case B_PAD_EDGE : pB[k] = *dataR; break;  // pad by edge value
         case B_CYCLE    : pB[k] = data;   break;  // cycle data
              default    :                 break;  // mirror or interpolate
      }
   }

   k = 0;
   for(i=nM; i<nR; i++) {

      if(m_Border != B_POLYNOM){
         sum = 0.; pB++;
         for(k=0; k<m_L/2; k++)
            sum += p_L[k] * (pB[k] + pB[m_L-1-k]);
      }
      else{
//       sum = hsum*Nevill(-0.5-nL, nS-i, ++pB, pBorder+m_L+1);
         k += 2;
//      sum = hsum*Nevill((m_L-k-1)/2., m_L-k, pB+k, pBorder+m_L+1); // wat version
         sum = hsum*Nevill((m_H-k-1)/2., m_H-k, pB+k, pBorder+m_H+1); // datacond version
      }

      *dataA += sum;
       dataA += stride;
   }

  delete [] pBorder;
}


// instantiations

#define CLASS_INSTANTIATION(class_) template class WaveDWT< class_ >;

CLASS_INSTANTIATION(float)
CLASS_INSTANTIATION(double)
//CLASS_INSTANTIATION(std::complex<float>)
//CLASS_INSTANTIATION(std::complex<double>)

#undef CLASS_INSTANTIATION

//} // namespace wat
//} // namespace datacondAPI
