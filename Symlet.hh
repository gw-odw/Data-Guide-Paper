// Wavelet Analysis Tool
//--------------------------------------------------------------------
// Implementation of 
// Symlet wavelets using Fast Wavelet Transform 
// References:
//   I.Daubechies, Ten lectures on wavelets
//   ISBN 0-89871-274-2, 1992
//--------------------------------------------------------------------

//$Id: Symlet.hh,v 0.2 2001/08/06 19:37:00 klimenko Exp $
#ifndef SYMLET_HH
#define SYMLET_HH

#include "WaveDWT.hh"

//namespace datacondAPI {
//namespace wat {

template<class DataType_t>
class Symlet : public WaveDWT<DataType_t>
{
   private:

      //: forward LP filter coefficients.
      double *pLForward;
      //: inverse LP filter coefficients.
      double *pLInverse;
      //: forward LP filter coefficients.
      double *pHForward;
      //: inverse LP filter coefficients.
      double *pHInverse;

      void setFilter();

   public:
      
      //: construct from base class
      Symlet(const Wavelet &);

      //: copy constructors
      Symlet(const Symlet<DataType_t> &);

      //: construct from wavelet parameters
      Symlet(int order=4, int tree=0, enum BORDER border=B_CYCLE);

      //: destructor
      virtual ~Symlet();

      //: Duplicate on heap
      virtual Symlet* Clone() const;

      //: decomposition method
      virtual void forward(int level, int layer);
      //: reconstruction method      
      virtual void inverse(int level, int layer);

}; // class Symlet

//}; // namespace wat
//}; // namespace datacondAPI

#endif // SYMLET_HH












