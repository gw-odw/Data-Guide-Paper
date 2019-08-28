/*************************************************************************
 
 Copyright (c) 2019 Neil Cornish
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 *************************************************************************/



#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_cdf.h>

#define PI 3.1415926535897932
#define TPI 6.2831853071795862319959269370884     // 2 Pi
#define SQPI 2.5066282746310002  // sqrt(TPI)

// gcc -o Qscans Qscans.c -lm -lgsl

double f_nwip(double *a, double *b, int n);
double fb_nwip(double *a, double *b, int n, int imin, int imax);
void phase_blind_time_shift(double *corr, double *corrf, double *data1, double *data2, int n);
void SineGaussianC(double *hs, double *sigpar, double Tobs, int NMAX);
void TransformC(double *a, double *freqs, double *tf, double *tfR, double *tfI, double Q, double Tobs, int n, int m);
void whiten(double *data, double *Sn, int N);
void tukey(double *data, double alpha, int N);


int main()
{
  int i, Nsamp;
  double SNR;
  double junk;
  double h1t, l1t, f, x, y, te;
  char filename[1024];
  int n, N, M;
  double *SH, *SL;
  double *SHR, *SLR;
  double *H, *L;
  double *HR, *LR;
  double Tobs, df, dt;
  FILE *h1;
  FILE *l1;
  FILE *time;
  FILE *out;
  FILE *in;
    
    // The BW files are from the TGR analysis from the time of the event
    // original data in /home/tyson/O1/GW150914/residuals_C01/job_1126259462_0.0/waveforms/
    // resisual analysis in CIT:/home/tyson/O1/GW150914/residuals_C01/IMRPhenomPv2
    
    Tobs = 4.0;
    df = 1.0/dt;
    dt = 1.0/1024.0;
    
    N = (int)(Tobs/dt);
    
    // frequency domain files don't start at zero frequency
    Nsamp = 2*1920;
    
    H=malloc(N*sizeof(double));
    L=malloc(N*sizeof(double));
    HR=malloc(N*sizeof(double));
    LR=malloc(N*sizeof(double));
    SH=malloc((N/2)*sizeof(double));
    SL=malloc((N/2)*sizeof(double));
    SHR=malloc((N/2)*sizeof(double));
    SLR=malloc((N/2)*sizeof(double));
    
   // set missing values of PSD to large number
    for(n=1; n< 128; n++)
    {
        SH[n] = 1.0;
        SL[n] = 1.0;
        SHR[n] = 1.0;
        SLR[n] = 1.0;
    }
    
    // Read in PSDs for the resdiual data
    
    h1 = fopen("clean_frequency_residual_199.dat.0","r");
    l1 = fopen("clean_frequency_residual_199.dat.1","r");
    for(n=128; n< N/2; n++)
    {
        fscanf(h1,"%lf%lf%lf%lf", &f, &x, &y, &SHR[n]);
        fscanf(l1,"%lf%lf%lf%lf", &f, &x, &y, &SLR[n]);
    }
    fclose(h1);
    fclose(l1);
    
    // Read in PSDs for the original data (the signal is removed by BW, but BL gives different PSD draw)
    
    h1 = fopen("clean_frequency_data_199.dat.0","r");
    l1 = fopen("clean_frequency_data_199.dat.1","r");
    for(n=128; n< N/2; n++)
    {
        fscanf(h1,"%lf%lf%lf%lf", &f, &x, &y, &SH[n]);
        fscanf(l1,"%lf%lf%lf%lf", &f, &x, &y, &SL[n]);
    }
    fclose(h1);
    fclose(l1);
    
    // Read in whitened time domain resdiual data
    
    h1 = fopen("clean_data_199.dat.0","r");
    l1 = fopen("clean_data_199.dat.1","r");
    for(n=0; n< N; n++)
    {
        fscanf(h1,"%lf%lf",  &x, &H[n]);
        fscanf(l1,"%lf%lf",  &x, &L[n]);
    }
    fclose(h1);
    fclose(l1);
    
    
    // Read in whitened time domain data
    
    h1 = fopen("clean_res_199.dat.0","r");
    l1 = fopen("clean_res_199.dat.1","r");
    for(n=0; n< N; n++)
    {
        fscanf(h1,"%lf%lf",  &x, &HR[n]);
        fscanf(l1,"%lf%lf",  &x, &LR[n]);
    }
    fclose(h1);
    fclose(l1);
    
    // FFT
    gsl_fft_real_radix2_transform(H, 1, N);
    gsl_fft_real_radix2_transform(L, 1, N);

    gsl_fft_real_radix2_transform(HR, 1, N);
    gsl_fft_real_radix2_transform(LR, 1, N);
    
    // FFT normaliztioon factor
    x = sqrt((double)(N)*16.0);
    
    for(n=0; n< N; n++)
    {
        HR[n] /= x;
        LR[n] /= x;
        H[n] /= x;
        L[n] /= x;
    }
    
    
    // re-color so raw data uses same PSD as the resdiual
    // Need to do this because BW stochastically samples the PSD - every draw is a little different
    for(n=1; n< N/2; n++)
    {
        x = sqrt(SH[n]/SHR[n]);
        H[n] *= x;
        H[N-n] *= x;
        x = sqrt(SL[n]/SLR[n]);
        L[n] *= x;
        L[N-n] *= x;
    }

    
    int subscale, octaves, Nf, j;
    double fmax, fmin;
    double *freqs, *sqf;
    double dx, dlnf, Q, t;
    
    fmin = 32.0;
    fmax = 512.0;
    
    
    // for logarithmic frequency spacing
    subscale = 64;  // number of semi-tones per octave
    octaves = (int)(rint(log(fmax/fmin)/log(2.0))); // number of octaves
    Nf = subscale*octaves+1;
    freqs = malloc(Nf*sizeof(double));   // frequencies used in the analysis
    sqf = malloc(Nf*sizeof(double));
    dx = log(2.0)/(double)(subscale);
    dlnf = dx;
    x = log(fmin);
    for(i=0; i< Nf; i++)
    {
        freqs[i] = exp(x);
        sqf[i] = sqrt(freqs[i]);
        x += dx;
    }
    
    Q = 8.0;
    
    double *tfH1R, *tfH1I;
    double *tfL1R, *tfL1I;
    double *tfH1, *tfL1;
    
    
    tfH1R = malloc((Nf*N)*sizeof(double));
    tfH1I = malloc((Nf*N)*sizeof(double));
    tfL1R = malloc((Nf*N)*sizeof(double));
    tfL1I = malloc((Nf*N)*sizeof(double));
    tfH1 = malloc((Nf*N)*sizeof(double));
    tfL1 = malloc((Nf*N)*sizeof(double));
    
    TransformC(HR, freqs, tfH1, tfH1R, tfH1I, Q, Tobs, N, Nf);
    TransformC(LR, freqs, tfL1, tfL1R, tfL1I, Q, Tobs, N, Nf);
    
    out = fopen("residualH1.dat","w");
    
    for(j = 0; j < Nf; j++)
    {
        f = freqs[j];
        
        for(i = 0; i < N; i++)
        {
            t = (double)(i)*dt;
            fprintf(out,"%e %e %e\n", t-Tobs/2.0, f, tfH1[j*N+i]);
        }
        
        fprintf(out,"\n");
    }
    fclose(out);
    
    
    out = fopen("residualL1.dat","w");
    for(j = 0; j < Nf; j++)
    {
        f = freqs[j];
        
        for(i = 0; i < N; i++)
        {
            t = (double)(i)*dt;
            fprintf(out,"%e %e %e\n", t-Tobs/2.0, f, tfL1[j*N+i]);
        }
        
        fprintf(out,"\n");
    }
    fclose(out);
    
    
    
    TransformC(H, freqs, tfH1, tfH1R, tfH1I, Q, Tobs, N, Nf);
    TransformC(L, freqs, tfL1, tfL1R, tfL1I, Q, Tobs, N, Nf);
    
    out = fopen("rawH1.dat","w");
    
    for(j = 0; j < Nf; j++)
    {
        f = freqs[j];
        
        for(i = 0; i < N; i++)
        {
            t = (double)(i)*dt;
            fprintf(out,"%e %e %e\n", t-Tobs/2.0, f, tfH1[j*N+i]);
        }
        
        fprintf(out,"\n");
    }
    fclose(out);
    
    
    out = fopen("rawL1.dat","w");
    for(j = 0; j < Nf; j++)
    {
        f = freqs[j];
        
        for(i = 0; i < N; i++)
        {
            t = (double)(i)*dt;
            fprintf(out,"%e %e %e\n", t-Tobs/2.0, f, tfL1[j*N+i]);
        }
        
        fprintf(out,"\n");
    }
    fclose(out);
    
    free(freqs);
    free(sqf);
    free(H);
    free(L);
    free(HR);
    free(LR);
    free(SH);
    free(SL);
    free(SHR);
    free(SLR);
    
    free(tfH1R);
    free(tfH1I);
    free(tfL1R);
    free(tfL1I);
    free(tfH1);
    free(tfL1);
    
}


void whiten(double *data, double *Sn, int N)
{
    double f, x, y, fix;
    int i;
    
    data[0] = 0.0;
    data[N/2] = 0.0;
    
    for(i=1; i< N/2; i++)
    {
        x = 1.0/sqrt(Sn[i]);
        data[i] *= x;
        data[N-i] *= x;
    }
    
}





void tukey(double *data, double alpha, int N)
{
    int i, imin, imax;
    double filter;
    
    imin = (int)(alpha*(double)(N-1)/2.0);
    imax = (int)((double)(N-1)*(1.0-alpha/2.0));
    
    for(i=0; i< N; i++)
    {
        filter = 1.0;
        if(i < imin) filter = 0.5*(1.0+cos(PI*( (double)(i)/(double)(imin)-1.0 )));
        if(i>imax) filter = 0.5*(1.0+cos(PI*( (double)(i)/(double)(imin)-2.0/alpha+1.0 )));
        data[i] *= filter;
    }
    
}


void TransformC(double *a, double *freqs, double *tf, double *tfR, double *tfI, double Q, double Tobs, int n, int m)
{
    double max, min;
    int i, j, k, p;
    int index;
    double tmp, mx;
    double f, t, dt, df, x, sqf;
    double *AC, *AF;
    double *b;
    double *params;
    double psh;
    double bmag, fix;
    
    // [0] t0 [1] f0 [2] Q [3] Amp [4] phi
    
    params= malloc(6*sizeof(double));
    
    params[0] = 0.0;
    params[2] = Q;
    params[3] = 1.0;
    params[4] = 0.0;
    
    dt = Tobs/(double)n;
    df = 1.0/Tobs;
    
    AC=malloc(n*sizeof(double));
    AF=malloc(n*sizeof(double));
    b = malloc(n*sizeof(double));
    
    fix = sqrt((double)(n/2));
    
    
    for(j = 0; j < m; j++)
    {
        
        f = freqs[j];
        
        params[1] = f;
        
        SineGaussianC(b, params, Tobs, n);
        
        bmag = sqrt(f_nwip(b, b, n)/(double)n);
        
        bmag /= fix;
        
        //printf("%f %f\n", f, bmag);
        
        phase_blind_time_shift(AC, AF, a, b, n);
        
        
        for(i = 0; i < n; i++)
        {
            
            tfR[j*n+i] = AC[i]/bmag;
            tfI[j*n+i] = AF[i]/bmag;
            tf[j*n+i] = tfR[j*n+i]*tfR[j*n+i]+tfI[j*n+i]*tfI[j*n+i];
        }
        
    }
    
    
    free(AC);
    free(AF);
    free(b);
    free(params);
    
}


void phase_blind_time_shift(double *corr, double *corrf, double *data1, double *data2, int n)
{
    int nb2, i, l, k, j;
    int imax, imin;
    
    nb2 = n / 2;
    
    corr[0] = 0.0;
    corrf[0] = 0.0;
    corr[nb2] = 0.0;
    corrf[nb2] = 0.0;
    
    for (i=1; i < nb2; i++)
    {
        l=i;
        k=n-i;
        
        corr[l]	= (data1[l]*data2[l] + data1[k]*data2[k]);
        corr[k]	= (data1[k]*data2[l] - data1[l]*data2[k]);
        corrf[l] = corr[k];
        corrf[k] = -corr[l];
    }
    
    gsl_fft_halfcomplex_radix2_inverse(corr, 1, n);
    gsl_fft_halfcomplex_radix2_inverse(corrf, 1, n);
    
    
}

double f_nwip(double *a, double *b, int n)
{
    int i, j, k;
    double arg, product;
    double test;
    double ReA, ReB, ImA, ImB;
    
    arg = 0.0;
    for(i=1; i<n/2; i++)
    {
        j = i;
        k = n-1;
        ReA = a[j]; ImA = a[k];
        ReB = b[j]; ImB = b[k];
        product = ReA*ReB + ImA*ImB;
        arg += product;
    }
    
    return(arg);
    
}

void SineGaussianC(double *hs, double *sigpar, double Tobs, int NMAX)
{
    double f0, t0, Q, sf, sx, Amp;
    double fmax, fmin, fac;
    double phi, f, t, x, y, z, zold;
    double tau, dt;
    int imin, imax;
    
    int i, id, N;
    
    // Torrence and Compo
    
    t0 = sigpar[0];
    f0 = sigpar[1];
    Q = sigpar[2];
    
    tau = Q/(TPI*f0);
    
    dt = Tobs/(double)(NMAX);
    
    fmax = f0 + 3.0/tau;  // no point evaluating waveform past this time (many efolds down)
    fmin = f0 - 3.0/tau;  // no point evaluating waveform before this time (many efolds down)
    
    fac = sqrt(sqrt(2.0)*PI*tau/dt);
    
    i = (int)(f0*Tobs);
    imin = (int)(fmin*Tobs);
    imax = (int)(fmax*Tobs);
    if(imax - imin < 10)
    {
        imin = i-5;
        imax = i+5;
    }
    
    if(imin < 0) imin = 0;
    if(imax > NMAX/2) imax = NMAX/2;
    
    hs[0] = 0.0;
    hs[NMAX/2] = 0.0;
    
    for(i = 1; i < NMAX/2; i++)
    {
        hs[i] = 0.0;
        hs[NMAX-i] = 0.0;
        
        if(i > imin && i < imax)
        {
            f = (double)(i)/Tobs;
            sf = fac*exp(-PI*PI*tau*tau*(f-f0)*(f-f0));
            hs[i] = sf;
            hs[NMAX-i] = 0.0;
        }
        
    }
    
    
}


double fb_nwip(double *a, double *b, int n, int imin, int imax)
{
    int i, j, k;
    double arg, product;
    double test;
    double ReA, ReB, ImA, ImB;
    
    arg = 0.0;
    for(i=1; i<n/2; i++)
    {
        if(i > imin && i < imax)
        {
            j = i;
            k = n-i;
            ReA = a[j]; ImA = a[k];
            ReB = b[j]; ImB = b[k];
            product = ReA*ReB + ImA*ImB;
            arg += product;
        }
    }
    
    return(arg);
    
}


