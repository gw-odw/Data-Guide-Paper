// Wavelet Analysis Tool
//--------------------------------------------------------------------
// Implementation of
// Meyer wavelets using Fast Wavelet Transform
// References:
//--------------------------------------------------------------------

//$Id: Meyer.hh,v 0.2 2001/08/06 19:37:00 klimenko Exp $
#ifndef MEYER_HH
#define MEYER_HH

#include "WaveDWT.hh"

//namespace datacondAPI {
//namespace wat {

template<class DataType_t>
class Meyer : public WaveDWT<DataType_t>
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
      Meyer(const Wavelet &);

      //: copy constructors
      Meyer(const Meyer<DataType_t> &);

      //: construct from wavelet parameters
      Meyer(int m, int tree=0, enum BORDER border=B_CYCLE);

      //: destructor
      virtual ~Meyer();

      //: Duplicate on heap
      virtual Meyer* Clone() const;

      //: calculate wavelet filter
      //!param: taper function order n
      //!param: beta(n,n) - value of Euler's beta function
      //!param: integration step
      double filter(int, double, double=1.e-6);

      // get maximum possible level of decomposition
      int getMaxLevel();

      //: decomposition method
      virtual void forward(int level, int layer);
      //: reconstruction method
      virtual void inverse(int level, int layer);

}; // class Meyer

//}; // namespace wat
//}; // namespace datacondAPI

#endif // MEYER_HH

