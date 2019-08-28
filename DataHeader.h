#include <gsl/gsl_math.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

void cc(double *Hf, double *Lf, double t1, double t2, double T, int M,  double *cr);

void bwbpf(double *in, double *out, int fwrv, int M, int n, double s, double f1, double f2);
void Align(double *a, double *b, int n, double T, double *tx, double *px, double *nx);
void shift(double *b, double T, double ts, double ps, double norm, int n);
void max_array_element(double *max, int *index, double *array, int n);
void cross_correlate(double *Hf, double *Lf, double t1, double t2, double T, int M,  FILE *out);

double maxlike(double *a, double *b, int n, double T);
void massscale(double Mscale, double *time, double *tmpl, double *wave, double dt, int N);

void tukey(double *data, double alpha, int N);
void tukey_scale(double *s1, double *s2, double alpha, int N);

void tukeyF(double *data, double T, double frise, double fstart, double fend, int N);

void CubicSplineGSL(int N, double *x, double *y, int Nint, double *xint, double *yint);

void whiten(double *data, double *Sn, int N);

void SineGaussianF(double *hs, double *sigpar, double Tobs, int NMAX);
void SineGaussianC(double *hs, double *sigpar, double Tobs, int NMAX);
void phase_blind_time_shift(double *corr, double *corrf, double *data1, double *data2, int n);
double f_nwip(double *a, double *b, int n);
void Transform(double *a, double *freqs, double **tf, double **tfR, double **tfI, double Q, double Tobs, int n, int m);
void TransformC(double *a, double *freqs, double **tf, double **tfR, double **tfI, double Q, double Tobs, int n, int m);
double Getscale(double *freqs, double Q, double Tobs, double fmax, int n, int m);
double fb_nwip(double *a, double *b, int n, int imin, int imax);


