/*-------------------------------------------------------
 * Package: 	Wavelet Analysis Tool
 * complex class 
 * File name: 	wavecomplex.cc
 *-------------------------------------------------------
*/


#include <time.h>
#include <iostream>
#include "wavecomplex.hh"

//: constructors
wavecomplex::wavecomplex() 
{ re = im = 0.; }
wavecomplex::wavecomplex(double a, double b)
{ re = a; im = b; }

//: copy constructor 
wavecomplex::wavecomplex(const wavecomplex& a)
{ *this = a; }

//: destructor
wavecomplex::~wavecomplex() {} 

//: operators
wavecomplex& wavecomplex::operator=(const wavecomplex& a)
{ re=a.real(); im=a.imag(); return *this; }

wavecomplex& wavecomplex::operator+=(const wavecomplex &a)
{ re+=a.real(); im+=a.imag(); return *this; }
wavecomplex& wavecomplex::operator-=(const wavecomplex &a)
{ re-=a.real(); im-=a.imag(); return *this; }
wavecomplex& wavecomplex::operator*=(const wavecomplex &a)
{
  double x = re*a.real()-im*a.imag();
  im = re*a.imag()+im*a.real(); 
  re = x;
  return *this;
}
wavecomplex& wavecomplex::operator/=(const wavecomplex &a)
{
  double x = a.abs();
  double y = (re*a.real()+im*a.imag())/x;
  im = (im*a.real()-re*a.imag())/x; 
  re = y;
  return *this;
}

wavecomplex wavecomplex::operator+(const wavecomplex &a)
{ wavecomplex z = *this; z+=a; return z; }
wavecomplex wavecomplex::operator-(const wavecomplex &a)
{ wavecomplex z = *this; z-=a; return z; }
wavecomplex wavecomplex::operator*(const wavecomplex &a)
{ wavecomplex z = *this; z*=a; return z; }
wavecomplex wavecomplex::operator/(const wavecomplex &a)
{ wavecomplex z = *this; z/=a; return z; }


wavecomplex& wavecomplex::operator=(const double c)
{ re=c; im=0.; return *this; }

wavecomplex& wavecomplex::operator+=(const double c)
{ re+=c; return *this; }
wavecomplex& wavecomplex::operator-=(const double c)
{ re-=c; return *this;}
wavecomplex& wavecomplex::operator*=(const double c)
{ re*=c; im*=c; return *this;}
wavecomplex& wavecomplex::operator/=(const double c)
{ re/=c; im/=c; return *this;}

wavecomplex wavecomplex::operator+(const double c)
{ wavecomplex z = *this; z+=c; return z;}
wavecomplex wavecomplex::operator-(const double c)
{ wavecomplex z = *this; z-=c; return z;}
wavecomplex wavecomplex::operator*(const double c)
{ wavecomplex z = *this; z*=c; return z;}
wavecomplex wavecomplex::operator/(const double c)
{ wavecomplex z = *this; z/=c; return z;}



