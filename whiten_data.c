/**************************************************************************
 
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
 
 ************************************************************************/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

// gcc -o whiten_data whiten_data.c -lm -lgsl


void spectrum(double *data, double *S, double *Sn, double *Smooth, double df, gsl_rng * r, int N);
void whiten(double *data, double *Sn, int N);
void tukey(double *data, double alpha, int N);

#include "Constants.h"

int main(int argc, char *argv[])
{
  int i, j, k, M, N, Nf, Nstep, Nclean, ii, m, rs, tsi, tti;
    int jj, kk, Nlines;
  int oflag, flag;
  int imin, imax;
  double SNR, max;
  double junk, Tobs, fix, f, t, t0, dt, dtm, df, x, y, z, dx;
  double fmax, fmin, dfx, Q, fny, delt, scale, dlnf;
  double Hmax, Lmax;
  double pshift;
  double *freqs, *data, *ref;
  double *inp, *oup, *slice;
  double *H1dat, *L1dat;
  double *Draw;
  double *D, *times;
  double *Dds;
  double *Sn;
  double *specD, *sspecD;
  double *sdata;
  double *intime, *sqf;
  double sigmean, sigmedian;
  int subscale, octaves;
    int mmax;
    double SNRsq, SNRold, pH, pL, pmax;
    double SNRH, SNRL, pw, alpha;
   double t_rise, s1, s2, ascale, fac;
    double av, var;
  double ttrig, tstart, tstart_clean, Tclean, starttime, endtime, Dfmax;
    
    int Oflag;
    
  int modelprint;
    
  double *linef, *linew, *lineh, *lineQ;
	
  char filename[1024];
  char command[1024];
    char Dname[1024];
   

  int n;
    
    const gsl_rng_type * P;
    gsl_rng * r;

    gsl_rng_env_setup();
    
    P = gsl_rng_default;
    r = gsl_rng_alloc (P);


  FILE *in;
  FILE *ifp;
  FILE *out;
    
    if(argc!=4)
    {
        printf("./whiten_data H/L Tobs trig_time\n");
        return 1;
    }
    
    Oflag = atoi(argv[1]);  // 0 for H1, 1 for L1
    Tobs = atof(argv[2]);  // duration
    ttrig = atof(argv[3]);  // trigger time
    tstart = ttrig - Tobs/2.0;
    

    
    starttime = tstart;
    endtime = tstart + Tobs;
    
    tti = (int)(ttrig);
    
    /* Here we read in the data, set up some arrays, Tukey window and FFT */
    
    sprintf(command, "frame_%d_%d_%d.dat", (int)(Tobs), tti, Oflag);
    in = fopen(command,"r");
    
    
    N = -1;
    while(!feof(in))
    {
        fscanf(in,"%lf%lf", &x, &y);
        N++;
    }
    rewind(in);
    
    dt = Tobs/(double)(N); // cadence
    df = 1.0/Tobs;  // frequency resolution
    fny = 1.0/(2.0*dt);  // Nyquist
    
    times = (double*)malloc(sizeof(double)*(N));
    D = (double*)malloc(sizeof(double)*(N));
    
    for (i = 0; i < N; ++i) fscanf(in,"%lf%lf", &times[i], &D[i]);
    fclose(in);
    
    // Tukey window parameter. Flat for (1-alpha) of data
    t_rise = 0.5; // Standard LAL setting
    alpha = (2.0*t_rise/Tobs);
    
    
    // Tukey window
    tukey(D, alpha, N);
   
    printf("Fourier Transforming Data\n");
    
    // FFT
    gsl_fft_real_radix2_transform(D, 1, N);
    
    Sn = (double*)malloc(sizeof(double)*(N/2));
    specD = (double*)malloc(sizeof(double)*(N/2));
    sspecD = (double*)malloc(sizeof(double)*(N/2));
    
    printf("Computing spectral estimate\n");
    
    // Form spectral model for whitening data (lines plus a smooth component)
    spectrum(D, Sn, specD, sspecD, df, r, N);
    
    // whiten data
    whiten(D, specD, N);
    
    fac = Tobs/((double)(N)*(double)(N));
    
     printf("Saving Results\n");
    
    sprintf(command, "PSD_%d_%d_%d.dat", Oflag, (int)(Tobs), tti);
    out = fopen(command,"w");
    for (i = 1; i < N/2; ++i)
    {
        fprintf(out,"%.15e %.15e %.15e %.15e\n", (double)(i)/Tobs, Sn[i]*fac, specD[i]*fac, sspecD[i]*fac);
    }
    fclose(out);
    
    sprintf(command, "freq_%d_%d_%d.dat", Oflag, (int)(Tobs), tti);
    out = fopen(command,"w");
    for (i = 1; i < N/2; ++i)
    {
        fprintf(out,"%.15e %.15e %.15e\n", (double)(i)/Tobs, D[i], D[N-i]);
    }
    fclose(out);
    
    gsl_fft_halfcomplex_radix2_inverse(D, 1, N);
    
    av = 0.0;
    var = 0.0;
    
    for (i = N/4; i < (N-N/4); ++i)
    {
        av += D[i];
        var += D[i]*D[i];
    }
    
    av /= (double)(N/2);
    var /= (double)(N/2);
    
    var = sqrt(var -av*av);
    
    //printf("Mean %e 1/Deviation %e  %e\n", av, 1.0/var, sqrt((double)(2*N)));
    
    
    fac = 1.0/var;
    
    
    sprintf(command, "time_%d_%d_%d.dat", Oflag, (int)(Tobs), tti);
    out = fopen(command,"w");
    for (i = 0; i < N; ++i)
    {
        fprintf(out,"%.15e %.15e\n", (double)(i)*dt, fac*D[i]);
    }
    fclose(out);
    
    free(D);
    free(times);
    free(Sn);
    free(specD);
    free(sspecD);
    
    
    
    return 0;

}


void spectrum(double *data, double *S, double *Sn, double *Smooth, double df, gsl_rng * r, int N)
{
    double Df, Dfmax, x, y;
    double Df1, Df2;
    int mw, k, i, j;
    int mm, kk;
    int end1, end2, end3;
    double med;
    double *chunk;
    

    // log(2) is median/2 of chi-squared with 2 dof
    
    
    for(i=1; i< N/2; i++) S[i] = 2.0*(data[i]*data[i]+data[N-i]*data[N-i]);
    S[0] = S[1];
    
    
    Dfmax = 16.0; // is the  width of smoothing window in Hz
    
    // Smaller windows used initially where the spectrum is steep
    Df2 = Dfmax/2.0;
    Df1 = Dfmax/4.0;
    
    // defines the ends of the segments where smaller windows are used
    end1 = (int)(16.0/df);
    end2 = 2*end1;
    
    mw = (int)(Dfmax/df)+1;  // size of median window
    //printf("numer of bins in smoothing window %d\n", mw);
    k = (mw+1)/2;
    chunk = (double*)malloc(sizeof(double)*(mw));
    
    end3 = N/2-k;  // end of final chunk
    
    // Fill the array so the ends are not empty - just to be safe
    for(i=0;i< N/2;i++)
    {
        Sn[i] = S[i];
        Smooth[i] = S[i];
    }
    
    mw = (int)(Df1/df)+1;  // size of median window
    k = (mw+1)/2;
    
    for(i=4; i< k; i++)
    {
        mm = i/2;
        kk = (mm+1)/2;
        
        for(j=0;j< mm;j++)
        {
            chunk[j] = S[i-kk+j];
        }
        
        Sn[i] = gsl_stats_median(chunk, 1, mm)/LN2;  // chi-squared with two degrees of freedom
        Smooth[i] = Sn[i];
        
    }
    
    
    i = k;
    do
    {
        for(j=0;j< mw;j++)
        {
            chunk[j] = S[i-k+j];
        }
        
        Sn[i] = gsl_stats_median(chunk, 1, mw)/LN2;  // chi-squared with two degrees of freedom
        Smooth[i] = Sn[i];
        
        i++;
        
    }while(i < end1);
    
    
    
    mw = (int)(Df2/df)+1;  // size of median window
    k = (mw+1)/2;
    
    
    do
    {
        for(j=0;j< mw;j++)
        {
            chunk[j] = S[i-k+j];
        }
        
        Sn[i] = gsl_stats_median(chunk, 1, mw)/LN2;  // chi-squared with two degrees of freedom
        Smooth[i] = Sn[i];
        
        i++;
        
    }while(i < end2);
    
    mw = (int)(Dfmax/df)+1;  // size of median window
    k = (mw+1)/2;
    
    do
    {
        for(j=0;j< mw;j++)
        {
            chunk[j] = S[i-k+j];
        }
        
        Sn[i] = gsl_stats_median(chunk, 1, mw)/LN2;  // chi-squared with two degrees of freedom
        Smooth[i] = Sn[i];
        
        i++;
        
    }while(i < end3);
    
    
    for(i=end3; i< N/2-4; i++)
    {
        mm = (N/2-i)/2;
        kk = (mm+1)/2;
        
        for(j=0;j< mm;j++)
        {
            chunk[j] = S[i-kk+j];
        }
        
        Sn[i] = gsl_stats_median(chunk, 1, mm)/LN2;  // chi-squared with two degrees of freedom
        Smooth[i] = Sn[i];
        
    }
    
    
    free(chunk);

    
    // zap the lines.
    for(i=1;i< N/2;i++)
    {
        x = S[i]/Sn[i];
        if(x > 10.0)
        {
            y = gsl_ran_exponential (r, 2.0);
            Sn[i] *= (x/y);
        }
    }
    
    
    
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


