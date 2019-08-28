// Wavelet Analysis Tool
//--------------------------------------------------------------------
// Implementation of 
// the Haar wavelet transform using lifting scheme 
// References:
//   A.Cohen, I.Daubechies, J.Feauveau Bases of compactly supported wavelets
//   Comm. Pure. Appl. Math. 45, 485-560, 1992
//--------------------------------------------------------------------

//$Id: Haar.hh,v 0.2 2001/08/06 19:37:00 klimenko Exp $
#ifndef HAAR_HH
#define HAAR_HH

#include "WaveDWT.hh"

//namespace datacondAPI {
//namespace wat {

template<class DataType_t>
class Haar : public WaveDWT<DataType_t>
{
   public:

      //: construct from wavelet parameters
      Haar(int tree=0);
      
      //: construct from the base class
      Haar(const Wavelet &);

      //: copy constructors
      Haar(const Haar<DataType_t> &);

      //: destructor

      virtual ~Haar();

      //: Duplicate on heap
      virtual Haar* Clone() const;

      //: decomposition method
      void forward(int level, int layer);
      //: reconstruction method      
      void inverse(int level, int layer);

}; // class Haar

//}; // namespace wat
//}; // namespace datacondAPI

#endif // HAAR_HH












