/**********************************************************
 * Package: 	Wavelet Analysis Tool
 * this class defines complex numbers
 * File name: 	wavecomplex.h
 **********************************************************/

#ifndef WAVECOMPLEX_HH
#define WAVECOMPLEX_HH

#include <cstdio>
#include <cmath>

class wavecomplex
{

  public:

  wavecomplex(double,double);                      // Constructor

  wavecomplex();                                   // Default constructor

  wavecomplex(const wavecomplex&);                       // copy Constructor

  virtual ~wavecomplex();                          // Destructor

// operators

          wavecomplex& operator= (const wavecomplex &);
  virtual wavecomplex& operator+=(const wavecomplex &);
  virtual wavecomplex& operator-=(const wavecomplex &);
  virtual wavecomplex& operator*=(const wavecomplex &);
  virtual wavecomplex& operator/=(const wavecomplex &);
  virtual wavecomplex  operator+ (const wavecomplex &);
  virtual wavecomplex  operator- (const wavecomplex &);
  virtual wavecomplex  operator* (const wavecomplex &);
  virtual wavecomplex  operator/ (const wavecomplex &);

          wavecomplex& operator= (const double);
  virtual wavecomplex& operator+=(const double);
  virtual wavecomplex& operator-=(const double);
  virtual wavecomplex& operator*=(const double);
  virtual wavecomplex& operator/=(const double);
  virtual wavecomplex  operator+ (const double);
  virtual wavecomplex  operator- (const double);
  virtual wavecomplex  operator* (const double);
  virtual wavecomplex  operator/ (const double);

// member functions

  inline double real() const { return re; }
  inline double imag() const { return im; }
  inline double  arg() const { return atan2(im,re); }
  inline double  abs() const { return re*re+im*im; }
  inline double  mod() const { return sqrt(re*re+im*im); }
  inline void    set(double x, double y) { re=x; im=y; return; }
  inline wavecomplex conj() { wavecomplex z(re,-im); return z; }

//   private:

  double re;		        // real
  double im;		        // imagenary

};

#endif // WAVECOMPLEX_HH












