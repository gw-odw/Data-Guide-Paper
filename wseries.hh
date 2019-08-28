// Wavelet Analysis Tool
// universal data container for wavelet transforms
// used with DMT and ROOT
//
//$Id: wseries.hh,v 1.1 2005/05/16 06:05:10 igor Exp $

#ifndef WSERIES_HH
#define WSERIES_HH

#ifndef WAVEARRAY_HH
#include "wavearray.hh"
#endif

#include <vector>
#include <list>
#include "WaveDWT.hh"
#include <complex>
#include "wavecomplex.hh"

typedef wavecomplex d_complex;
typedef std::vector<int> vector_int;


template<class DataType_t>
class WSeries : public wavearray<DataType_t>
{
    
   public:
      
      // constructors
      
      //: Default constructor
      WSeries();
    
      //: Construct WSeries for specific wavelet type 
      //+ default constructor
      explicit WSeries(const Wavelet &w);

      //: Construct from wavearray
      //!param: value - data to initialize the WSeries object
      explicit WSeries(const wavearray<DataType_t>& value,
		       const Wavelet &w);
    
      //: Copy constructor
      //!param: value - object to copy from 
      WSeries(const WSeries<DataType_t>& value);
      
      //: destructor
      virtual ~WSeries();
    
      // operators

      WSeries<DataType_t>& operator= (const wavearray<DataType_t> &);
      WSeries<DataType_t>& operator= (const WSeries<DataType_t> &);
      WSeries<DataType_t>& operator= (const DataType_t);

      // operator[](const slice &) sets the Slice object of the wavearray class
      virtual WSeries<DataType_t>& operator[](const std::slice &);

      //: operators for WSeries objects, which can have different length
      //: it is required they have the same type of transform (standard or binary)
      //: and the same size of approximation levels.
      //: warning: there is no check that approximation levels have the same 
      //: sampling rate.
      virtual WSeries<DataType_t>& operator+=(WSeries<DataType_t> &);
      virtual WSeries<DataType_t>& operator-=(WSeries<DataType_t> &);
      virtual WSeries<DataType_t>& operator*=(WSeries<DataType_t> &);

      // just to trick ANSI standard
      virtual WSeries<DataType_t>& operator+=(wavearray<DataType_t> &);
      virtual WSeries<DataType_t>& operator-=(wavearray<DataType_t> &);
      virtual WSeries<DataType_t>& operator*=(wavearray<DataType_t> &);
      virtual WSeries<DataType_t>& operator+=(const DataType_t);
      virtual WSeries<DataType_t>& operator-=(const DataType_t);
      virtual WSeries<DataType_t>& operator*=(const DataType_t);

      //: Dump data array to an ASCII file 
      virtual void Dump(const char*, int=0);

      // accessors

      //: Get maximum possible level of decompostion
      int getMaxLevel();

      //: Get level of decompostion
      inline int getLevel(){ return pWavelet->m_Level; }

      //: set level of decompostion
      inline void setLevel(size_t n){ pWavelet->m_Level = n; }

      //: Set black pixel probability
      inline void setbpp(double f) { bpp=f; return; }
      //: Get black pixel probability
      inline double getbpp() const { return bpp; }

      //: Set low frequency boundary
      inline void setlow(double f) { 
	 f_low=f>0. ? f : 0.; return; 
      }
      //: Get low frequency boundary
      inline double getlow() const { return f_low; }

      //: Set high frequency boundary
      inline void sethigh(double f) { 
	 f_high = f; return; 
      }
      //: get high frequency boundary
      inline double gethigh() const { return f_high; }

      //: Get max layer of decompostion
      inline int maxLayer(){ 
	 return pWavelet->BinaryTree() ? (1<<getLevel())-1 : getLevel();
      }

      //: Get slice structure for specified layer
      inline std::slice getSlice(size_t n=0){ return pWavelet->getSlice(n); }

      //: Get wavelet layer frequency resolution
      inline double resolution(size_t n){ 
	 return this->rate()*pWavelet->getSlice(n).size()/this->size()/2.; 
      }

      //: Extract wavelet coefficients from specified layer
      //!param: n - layer number
      int getLayer(wavearray<DataType_t> &, int n);

      //: replace wavelet data for specified layer with data from Sequence
      //!param: n - layer number
      void putLayer(wavearray<DataType_t> &, int n);

      // mutators

      virtual void   resize(unsigned int);
      virtual void   resample(double, int=6);

      //: initialize wavelet parameters from Wavelet object
      void setWavelet(const Wavelet &w);

      //: Perform n steps of forward wavelet transform
      //!param: n - number of steps (-1 means full decomposition)
      void Forward(int n = -1);
      void Forward(wavearray<DataType_t> &, int n = -1);
      void Forward(wavearray<DataType_t> &, Wavelet &, int n = -1);

      //: Perform n steps of inverse wavelet transform
      //!param: n - number of steps (-1 means full reconstruction)
      void Inverse(int n = -1);


  //++++++++++++++ wavelet data conditioning +++++++++++++++++++++++++++

      //: calculate running medians with window of t seconds
      //: and subtract the median from this (false key) or
      //: calculate median for abs(this) and normalize this (true)
      virtual void median(double, bool=false);


      //: apply linear predictor to each layer. 
      // param: filter length in seconds
      // param: filter mode: -1/0/1 - backward/symmetric/forward
      // param: stride for filter training (0 - train on whole TS)
      // param: boundary offset to account for wavelet artifacts
      virtual void lprFilter(double,int=0,double=0.,double=0.);

      //: tracking of noise non-stationarity and whitening. 
      // param 1 - time window dT. if = 0 - dT=T, where T is wavearray duration
      // param 2 - 0 - no whitening, 1 - single whitening, >1 - double whitening
      // param 3 - boundary offset 
      // param 4 - noise sampling interval (window stride)  
      //           the number of measurements is k=int((T-2*offset)/stride)
      //           if stride=0, then stride is set to dT
      // return: noise array if param2>0, median if param2=0
      //!what it does: each wavelet layer is devided into k intervals.
      //!The data for each interval is sorted and the following parameters 
      //!are calculated: median and the amplitude 
      //!corresponding to 31% percentile (wp). Wavelet amplitudes (w) are 
      //!normalized as  w' = (w-median(t))/wp(t), where median(t) and wp(t)
      //!is a linear interpolation between (median,wp) measurements for
      //!each interval. 
      virtual WSeries<double> white(double=0.,size_t=1,double=0.,double=0.);

      //: whiten TF map by using input noise array
      virtual bool white(WSeries<double>);

      //: local whitening, works only for binary wavelets. 
      //: returns array of noise rms for wavelet layers
      //!param: n - number of decomposition steps 
      //!algorithm: 
      //! 1) do forward wavelet transform with n decomposition steps
      //! 2) whiten wavelet layers and calculate noise rms as
      //!    1/Sum(1/var)
      //! 3) do inverse wavelet transform with n reconstruction steps
      virtual wavearray<double> filter(size_t);
    
      //: works only for binary wavelets. 
      //: calculates, corrects and returns noise variability 
      //!param: first  - time window to calculate normalization constants
      //!       second - low frequency boundary for correction
      //!       third  - high frequency boundary for correction
      //!algorithm: 
      //! 1) sort wavelet amplitudes with the same time stamp
      //! 2) calculate left(p) and right(p) amplitudes
      //!    put (right(p)-left(p))/2 into output array
      //! 3) if first parameter >0 - devide WSeries by average
      //!    variability
      virtual WSeries<float> variability(double=0., double=-1., double=-1.);

      //: Selection of a fixed fraction of pixels
      //: reduced wavelet amplitudes are stored in this   
      //: Returns fraction of non-zero coefficients.
      //!param: t - sub interval duration. If can not divide on integer
      //            number of sub-intervals then add leftover to the last
      //            one.
      //!param: f - black pixel fraction
      //!param: m - mode
      //!options: f = 0, m = 0 - returns black pixel occupancy  
      //!         m = 1 - set threshold f
      //!         m = 2 - random policy   
      //!         m = 0 - random pixel selection
      virtual double fraction(double=0.,double=0.,int=0);

      //: calculate rank logarithic significance of wavelet pixels
      //: reduced wavelet amplitudes are stored in this   
      //: Returns pixel occupancy for significance>0.
      //!param: n - sub-interval duration in seconds
      //!param: f - black pixel fraction
      //!options: f = 0 - returns black pixel occupancy  
      virtual double significance(double, double=1.);

      //: calculate running rank logarithic significance of wavelet pixels
      //: reduced wavelet amplitudes are stored in this   
      //: Returns pixel occupancy for significance>0.
      //!param: n - sub-interval duration in domain units
      //!param: f - black pixel fraction
      //!options: f = 0 - returns black pixel occupancy  
      virtual double rsignificance(size_t=0, double=1.);

      //: calculate running rank logarithic significance of wavelet pixels
      //: reduced wavelet amplitudes are stored in this   
      //: Returns pixel occupancy for significance>0.
      //!param: T - sliding window duration in seconds
      //!param: f - black pixel fraction
      //!param: t - sliding step in seconds
      //!options: f = 0 - returns black pixel occupancy  
      //!options: t = 0 - sliding step = wavelet time resolution.  
      virtual double rSignificance(double, double=1., double=0.);
    
      //: calculate running logarithic significance of wavelet pixels
      //: reduced wavelet amplitudes are stored in this   
      //: Returns pixel occupancy for significance>0.
      //!param: T - sliding window duration in seconds
      //!param: f - black pixel fraction
      //!param: t - sliding step in seconds
      //!options: f = 0 - returns black pixel occupancy  
      //!options: t = 0 - sliding step = wavelet time resolution.  
      virtual double gSignificance(double, double=1., double=0.);
    
      //: Selection of a fixed fraction of pixels for each wavelet layer 
      //: Returns fraction of non-zero coefficients.
      //!param: f - black pixel fraction
      //!param: m - mode
      //!options: f = 0 - returns black pixel occupancy  
      //!         m = 1 - set threshold f, returns percentile amplitudes 
      //!         m =-1 - set threshold f, returns wavelet amplitudes
      //!         m > 1 - random policy,returns percentile amplitudes   
      //!         m <-1 - random policy,returns wavelet amplitudes   
      //!         m = 0 - random pixel selection
      //! if m<0 return wavelet amplitudes instead of the percentile amplitude
      virtual double percentile(double=0.,int=0, WSeries<DataType_t>* = NULL);

      //: clean up single pixels   
      //!param: S - threshold on pixel significance
      //!return pixel occupancy.
      virtual double pixclean(double=0.);

      //: select pixels from *this which satisfy a coincidence rule 
      //: within specified window
      //!param: WSeries object used for coincidence 
      //!param: coincidence window in seconds
      //!return pixel occupancy 
      virtual double coincidence(WSeries<DataType_t> &, int=0, int=0, double=0.);

      //: select pixels from *this which satisfy a coincidence rule 
      //: within specified window w, above threshold T, 
      //!param: WSeries object used for coincidence 
      //!param: coincidence window in seconds
      //!param: threshold on significance
      //!return pixel occupancy 
      virtual double Coincidence(WSeries<DataType_t> &, double=0., double=0.);

      //: calculate calibration coefficients and apply energy calibration
      //: to wavelet data in this 
      //: AS_Q calibration: R(f)=(1 + gamma*(R*C-1))/alpha*C
      //: DARM calibration: R(f)=(1 + gamma*(R*C-1))/gamma*C
      //: input 
      //!param: number of samples in calibration arrays R & C 
      //!param: frequency resolution
      //!param: pointer to response function R in Fourier domain
      //!param: pointer to sensing function C in Fourier domain
      //!param: time dependent calibration coefficient alpha 
      //!param: time dependent calibration coefficient gamma 
      //!param: 0/1 - AS_Q/DARM_ERR calibration, by default is 0 
      //!return array with calibration constants for each wavelet layer
      virtual WSeries<double> calibrate(size_t, double,
					d_complex*, d_complex*,
					wavearray<double> &,
					wavearray<double> &,
					size_t ch=0);


      //: mask WSeries data with the pixel mask defined in wavecluster object
      //!param: int n 
      //! if n<0,  zero pixels defined in mask (regression)
      //! if n>=0, zero all pixels except ones defined in the mask
      //!param: bool  - if true, set WSeries data to be positive
      //! if pMask.size()=0, mask(0,true) is equivalent to abs(data)
      //!return core pixel occupancy 
      //virtual double mask(int=1, bool=false);

      //: calculate pixel occupancy for clusters passed selection cuts
      //!param: true - core pccupancy, false - total occupancy;
      //!return wavearray<double> with occupancy.
      //virtual wavearray<double> occupancy(bool=true);
      
// data members

      //: parameters of wavelet transform
      WaveDWT<DataType_t>* pWavelet;

      //: whitening mode
      size_t w_mode;

   private:
#ifndef __CINT__
      //: black pixel probability
      double bpp;
      //: low frequency boundary
      double f_low;
      //: high frequency boundary
      double f_high;
#endif

}; // class WSeries<DataType_t>


#endif // WSERIES_HH


















