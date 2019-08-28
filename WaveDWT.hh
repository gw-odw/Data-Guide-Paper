// Wavelet Analysis Tool
//$Id: WaveDWT.hh,v 1.3 2001/12/15 03:27:29 jzweizig Exp $
#ifndef WAVEDWT_HH
#define WAVEDWT_HH

//#include "slice.h"
//#include <valarray>
#ifndef __CINT__
#include <valarray>
#else
namespace std {
  class slice;
}
#endif
#include "Wavelet.hh"

//namespace datacondAPI {
//namespace wat {

template<class DataType_t> 
class WaveDWT : public Wavelet
{
   public: 

      //: constructor
      WaveDWT(int mH=1, int mL=1, int tree=0, enum BORDER border=B_CYCLE);

      //: construct from the base class
      WaveDWT(const Wavelet &);

      //: copy constructor
      WaveDWT(const WaveDWT<DataType_t> &);

      //: Destructor
      virtual ~WaveDWT();
    
      //: Duplicate on heap
      virtual WaveDWT* Clone() const;

      //: get maximum possible level of wavelet decompostion
      virtual int getMaxLevel();
      virtual int getMaxLevel (int i) {
         return Wavelet::getMaxLevel(i); }

      //: make slice for layer with specified index
      virtual std::slice getSlice(const int);

      //: make slice for (level,layer)
      virtual std::slice getSlice(const int, const int);

      //: Allocate data (set pWWS)
      bool allocate(size_t, DataType_t *);

      //: return allocate status (true if allocated)
      bool allocate();

      //: Release data
      void release();

      //: forward wavelet transform
      virtual void t2w(int=1);
      //: inverse wavelet transform
      virtual void w2t(int=1);

      //: makes one FWT decomposition step
      virtual void forwardFWT(int, int,
			   const double*,
			   const double*);
      //: makes one FWT reconstruction step
      virtual void inverseFWT(int, int,
			   const double*,
			   const double*);

      //: makes one prediction step for Lifting Wavelet Transform
      virtual void predict(int,int,const double*);
      //: makes one update step for Lifting Wavelet Transform
      virtual void update(int,int,const double*);

      //: virtual functions for derived wavelet classes

      //: makes one FWT decomposition step
      virtual void forward(int,int){}
      //: makes one FWT reconstruction step
      virtual void inverse(int,int){}

      DataType_t *pWWS;     // pointer to wavelet work space      
      unsigned int nWWS;    // size of the wavelet work space      

}; // class WaveDWT


//}; // namespace wat
//}; // namespace datacondAPI

#endif // WAVEDWT_HH















