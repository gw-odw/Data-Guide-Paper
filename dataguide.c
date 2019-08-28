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

#include "DataHeader.h"
#include "Constants.h"

// gcc -o dataguide dataguide.c -lm -lgsl

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
  double *DHraw, *DLraw;
  double *DH, *DL;
  double *Dds;
  double *SnH, *SnL;
  double *SwH, *SwL;
  double *sdata;
  double *intime, *sqf;
  double sigmean, sigmedian;
  int subscale, octaves;
  int mmax;
  double SNRsq, SNRold, pmax;
    double SNRH, SNRL, pw, alpha;
   double t_rise, s1, s2, ascale, fac, Tpad;
    double av, var, Tfull;
  double ttrig, tstart, tstart_clean, Tclean, starttime, endtime, Dfmax;
    
    int Oflag;
    
  int modelprint;
    
   double *linef, *linew, *lineh, *lineQ;
    
    double *DHfull, *DLfull;
    double *DHcopy, *DLcopy;
	
  char filename[1024];
  char command[1024];
  char Dname[1024];
   

  int n;
    
  FILE *in;
  FILE *ifp;
  FILE *out;
    
    Oflag = 0;

    Tobs = 4.0;
    starttime = 1126257414.0;  // start of data file
    ttrig = 1126259462.0;  // trigger time
    
    // Use multiple segements of length 4 seconds to get the spectral estimate
    Nf = 256;
    Tfull = (double)(Nf)*Tobs;
    
    tti = (int)(ttrig);

        dt= 1.0/4096.0;  // cadence
        df = 1.0/Tobs;  // frequency resolution
        N = (int)(Tobs/dt);  // Number of time samples in segment
        M = (int)(Tfull/dt);  // Allowing for padding
        fny = 1.0/(2.0*dt);  // Nyquist
    
        // printf("%d %d\n", N, M);
    
        DH = (double*)malloc(sizeof(double)*(N));
        DL = (double*)malloc(sizeof(double)*(N));
        DHfull = (double*)malloc(sizeof(double)*(M));
        DLfull = (double*)malloc(sizeof(double)*(M));
        DHcopy = (double*)malloc(sizeof(double)*(M));
        DLcopy = (double*)malloc(sizeof(double)*(M));
    
    /* # Gravitational wave strain for GW150914 for L1 & H1
     # This file has 4096 samples per second
     # starting GPS 1126257414 duration 4096 */
    
    k = (4096*4096);  // length of data file
    j = (int)((ttrig-starttime-Tfull/2.0)/dt); // where we start using the data
    
    printf("Reading in the Hanford data\n");
    if ((in = fopen("H-H1_LOSC_4_V2-1126257414-4096.txt","r")) == NULL)
    {
        printf("Error! opening file");
        exit(1);
    }
    // strip the header
    fgets(command, 1024, in);
    fgets(command, 1024, in);
    fgets(command, 1024, in);
    //skip the first j samples
    for (i = 0; i < j; ++i) fscanf(in,"%lf", &x);
    //read in the data
    for (i = 0; i < M; ++i) fscanf(in,"%lf", &DHfull[i]);
    fclose(in);
    
    printf("Reading in the Livingston data\n");
    if ((in = fopen("L-L1_LOSC_4_V2-1126257414-4096.txt","r")) == NULL)
    {
        printf("Error! opening file");
        exit(1);
    }
    // strip the header
    fgets(command, 1024, in);
    fgets(command, 1024, in);
    fgets(command, 1024, in);
    //skip the first j samples
    for (i = 0; i < j; ++i) fscanf(in,"%lf", &x);
    //read in the data
    for (i = 0; i < M; ++i) fscanf(in,"%lf", &DLfull[i]);
    fclose(in);
    
    printf("Data read complete\n");
    
    // save a copy of the original data
    for (i = 0; i < M; ++i) DHcopy[i] = DHfull[i];
    for (i = 0; i < M; ++i) DLcopy[i] = DLfull[i];
    
    // select chunk of duration Tobs
    // Tukey window parameter. Flat for (1-alpha) of data
    
    t_rise = 0.5; // Standard LAL setting
    alpha = (2.0*t_rise/Tobs);
    
    tukey_scale(&s1, &s2, alpha, N);
    
    fac = Tobs/((double)(N)*(double)(N));
    
    SnH = (double*)malloc(sizeof(double)*(N/2));
    SwH = (double*)malloc(sizeof(double)*(N/2));
    SnL = (double*)malloc(sizeof(double)*(N/2));
    SwL = (double*)malloc(sizeof(double)*(N/2));
    
    for (i = 0; i < N/2; ++i) SwH[i] = 0.0;
    for (i = 0; i < N/2; ++i) SwL[i] = 0.0;
    
    k = 0;
    ii = 0;
    m = 0;
    for (j = 0; j < (2*Nf-1); ++j)
    {
        for (i = 0; i < N; ++i) DH[i] = DHfull[i+k];
        for (i = 0; i < N; ++i) DL[i] = DLfull[i+k];
        k += N/2;   // Move over a half block
        if(k !=  M/2-N/2) // skipping the block containing the signal - offsource estimate
        {
        tukey(DH, alpha, N);
        tukey(DL, alpha, N);
        gsl_fft_real_radix2_transform(DH, 1, N);
        gsl_fft_real_radix2_transform(DL, 1, N);
        for(i=1; i< N/2; i++) SnH[i] = 2.0*(DH[i]*DH[i]+DH[N-i]*DH[N-i]);
        for(i=1; i< N/2; i++) SnL[i] = 2.0*(DL[i]*DL[i]+DL[N-i]*DL[N-i]);
        for(i=1; i< N/2; i++) SwH[i] += SnH[i];
        for(i=1; i< N/2; i++) SwL[i] += SnL[i];
        ii++;
        }
        m++;
    }
    
    for(i=1; i< N/2; i++) SwH[i] /= ((double)(ii));
    for(i=1; i< N/2; i++) SwL[i] /= ((double)(ii));
    
    printf("Saving Spectra\n");
    
    sprintf(command, "Welch_%d_%d_%d.dat", Oflag, (int)(Tobs), tti);
    out = fopen(command,"w");
    for (i = 1; i < N/2; ++i)
    {
        fprintf(out,"%.15e %.15e %.15e\n", (double)(i)/Tobs, fac*SwH[i]/s2, fac*SwL[i]/s2);
    }
    fclose(out);
    
    
    // select central block
    k = M/2-N/2;
    for (i = 0; i < N; ++i) DH[i] = DHfull[i+k];
    for (i = 0; i < N; ++i) DL[i] = DLfull[i+k];
    
    // Tukey window
    tukey(DH, alpha, N);
    tukey(DL, alpha, N);
    
     // FFT
    gsl_fft_real_radix2_transform(DH, 1, N);
    gsl_fft_real_radix2_transform(DL, 1, N);
    
    sprintf(command, "pspec_%d_%d_%d.dat", Oflag, (int)(Tobs), tti);
    out = fopen(command,"w");
    for (i = 0; i < N/2; ++i)
    {
        fprintf(out,"%.15e %.15e %.15e\n", (double)(i)/Tobs, fac*2.0*(DH[i]*DH[i]+DH[N-i]*DH[N-i])/s2, fac*2.0*(DL[i]*DL[i]+DL[N-i]*DL[N-i])/s2);
    }
    fclose(out);
    
    
    
    // select central block
    k = M/2-N/2;
    for (i = 0; i < N; ++i) DH[i] = DHcopy[i+k];
    for (i = 0; i < N; ++i) DL[i] = DLcopy[i+k];
    
    // FFT
    gsl_fft_real_radix2_transform(DH, 1, N);
    gsl_fft_real_radix2_transform(DL, 1, N);
    
    sprintf(command, "freqnowindow_%d_%d_%d.dat", Oflag, (int)(Tobs), tti);
    out = fopen(command,"w");
    for (i = 0; i < N/2; ++i)
    {
        x = atan2(DH[N-i],DH[i]);
        if(x < 0.0) x += TPI;
        y = atan2(DL[N-i],DL[i]);
        if(y < 0.0) y += TPI;
        fprintf(out,"%.15e %.15e %.15e %.15e %.15e\n", (double)(i)/Tobs, fac*2.0*(DH[i]*DH[i]+DH[N-i]*DH[N-i]), x, fac*2.0*(DL[i]*DL[i]+DL[N-i]*DL[N-i]), y);
    }
    fclose(out);
    
    
    // select central block
    k = M/2-N/2;
    for (i = 0; i < N; ++i) DH[i] = DHcopy[i+k];
    for (i = 0; i < N; ++i) DL[i] = DLcopy[i+k];
    
    sprintf(command, "raw_%d_%d_%d.dat", Oflag, (int)(Tobs), tti);
    out = fopen(command,"w");
    for (i = 0; i < N; ++i)
    {
        fprintf(out,"%.15e %.15e %.15e\n", (double)(i-N/2)*dt, DH[i], DL[i]);
    }
    fclose(out);
    
    // Tukey window
    tukey(DH, alpha, N);
    tukey(DL, alpha, N);
    
    sprintf(command, "windowed_%d_%d_%d.dat", Oflag, (int)(Tobs), tti);
    out = fopen(command,"w");
    for (i = 0; i < N; ++i)
    {
        fprintf(out,"%.15e %.15e %.15e\n", (double)(i-N/2)*dt, DH[i], DL[i]);
    }
    fclose(out);
    
    
    printf("Fourier Transforming Data\n");
    
    // FFT
    gsl_fft_real_radix2_transform(DH, 1, N);
    gsl_fft_real_radix2_transform(DL, 1, N);
    
    // whiten data
    whiten(DH, SwH, N);
    whiten(DL, SwL, N);
    
    sprintf(command, "freq_%d_%d_%d.dat", Oflag, (int)(Tobs), tti);
    out = fopen(command,"w");
    for (i = 0; i < N/2; ++i)
    {
        x = atan2(DH[N-i],DH[i]);
        if(x < 0.0) x += TPI;
        y = atan2(DL[N-i],DL[i]);
        if(y < 0.0) y += TPI;
        fprintf(out,"%.15e %.15e %.15e %.15e %.15e\n", (double)(i)/Tobs, fac*2.0*(DH[i]*DH[i]+DH[N-i]*DH[N-i]), x, fac*2.0*(DL[i]*DL[i]+DL[N-i]*DL[N-i]), y);
    }
    fclose(out);
    
    gsl_fft_halfcomplex_radix2_inverse(DH, 1, N);
    gsl_fft_halfcomplex_radix2_inverse(DL, 1, N);
    
    av = 0.0;
    var = 0.0;
    
    for (i = N/4; i < (N-N/4); ++i)
    {
        av += DH[i]+DL[i];
        var += DH[i]*DH[i]+DL[i]*DL[i];
    }
    
    av /= (double)(N);
    var /= (double)(N);
    
    var = sqrt(var-av*av);
    
    
    sprintf(command, "time_%d_%d_%d.dat", Oflag, (int)(Tobs), tti);
    out = fopen(command,"w");
    for (i = 0; i < N; ++i)
    {
        DH[i] /= var;
        DL[i] /= var;
        fprintf(out,"%.15e %.15e %.15e\n", (double)(i-N/2)*dt, DH[i], DL[i]);
    }
    fclose(out);
    
    
    // Now bandspass to use same analysis range as the PRL
    
    fmax = 350.0;
    fmin = 35.0;
    
    bwbpf(DH, DH, 1, N, 8, 1.0/dt, fmax, fmin);
    bwbpf(DH, DH, -1, N, 8, 1.0/dt, fmax, fmin);
    bwbpf(DL, DL, 1, N, 8, 1.0/dt, fmax, fmin);
    bwbpf(DL, DL, -1, N, 8, 1.0/dt, fmax, fmin);
    
    sprintf(command, "bandpassed_%d_%d_%d.dat", Oflag, (int)(Tobs), tti);
    out = fopen(command,"w");
    for (i = 0; i < N; ++i)
    {
        fprintf(out,"%.15e %.15e %.15e\n", (double)(i-N/2)*dt, DH[i], DL[i]);
    }
    fclose(out);

    
    
    // The filter decays in time with an efold scale of 1/fmin. Need to throw out several/fmin
    // zero-phase Butterworth band-pass filter (zero phase by going forward/back)
    bwbpf(DHfull, DHfull, 1, M, 8, 1.0/dt, fmax, fmin);
    bwbpf(DHfull, DHfull, -1, M, 8, 1.0/dt, fmax, fmin);
    bwbpf(DLfull, DLfull, 1, M, 8, 1.0/dt, fmax, fmin);
    bwbpf(DLfull, DLfull, -1, M, 8, 1.0/dt, fmax, fmin);
    
    
    // Recompute variance for whitened bandpassed data
    // select central block
    k = M/2-N/2;
    for (i = 0; i < N; ++i) DH[i] = DHfull[i+k];
    for (i = 0; i < N; ++i) DL[i] = DLfull[i+k];
    // Tukey window
    tukey(DH, alpha, N);
    tukey(DL, alpha, N);
    // FFT
    gsl_fft_real_radix2_transform(DH, 1, N);
    gsl_fft_real_radix2_transform(DL, 1, N);
    // whiten data
    whiten(DH, SwH, N);
    whiten(DL, SwL, N);
    gsl_fft_halfcomplex_radix2_inverse(DH, 1, N);
    gsl_fft_halfcomplex_radix2_inverse(DL, 1, N);
    
    av = 0.0;
    var = 0.0;
    
    for (i = N/4; i < (N-N/4); ++i)
    {
        av += DH[i]+DL[i];
        var += DH[i]*DH[i]+DL[i]*DL[i];
    }
    
    av /= (double)(N);
    var /= (double)(N);
    var = sqrt(var-av*av);
    
    sprintf(command, "timebp_%d_%d_%d.dat", Oflag, (int)(Tobs), tti);
    out = fopen(command,"w");
    for (i = 0; i < N; ++i)
    {
        DH[i] /= var;
        DL[i] /= var;
        fprintf(out,"%.15e %.15e %.15e\n", (double)(i-N/2)*dt, DH[i], DL[i]);
    }
    fclose(out);
    
    
    double *tmplH, *tmplL, *tmplT, *tmplX;
    int NT;
    double dtt;
    
    // NR template starts at 0.25 seconds and is sampled at 16384 Hz, have to downsample to 4096 Hz and
    // zero pad
    
    NT = 860;
    
    tmplT = (double*)malloc(sizeof(double)*(N));
    tmplH = (double*)malloc(sizeof(double)*(N));
    tmplL = (double*)malloc(sizeof(double)*(N));
    tmplX = (double*)malloc(sizeof(double)*(N));
    
    for(i=0; i< N; i++)
    {
        tmplT[i] = 0.0;
    }
    
    // https://www.gw-openscience.org/GW150914data/P150914/fig2-unfiltered-waveform-H.txt
    if ((in = fopen("fig2-unfiltered-waveform-H.txt","r")) == NULL)
    {
        printf("Error! opening file");
        exit(1);
    }
    // strip the header
    fgets(command, 1024, in);
    
    j = 0;
    for(i=0; i< NT*4; i++)
    {
        fscanf(in,"%lf%lf", &x, &y);
        if(i%4==0)
        {
        y *= 1.0e-21;
        z = (double)(j);
        x = 1.0;
        if(j < 100)  // window the abrupt start
        {
            x = sin((z/100.0)*PI/2.0);
        }
        x = x*x;
        tmplT[j+N/2+1024] = x*y;
        j++;
        }
    }
    fclose(in);

    tukey(tmplT, alpha, N);
    
    for (i = 0; i < N; ++i)
    {
        tmplL[i] = tmplT[i];
        tmplH[i] = tmplT[i];
    }
    
    out = fopen("ref_template.dat","w");
    for (i = 0; i < N; ++i)
    {
        fprintf(out,"%.15e %.15e\n", (double)(i-N/2)*dt, tmplH[i]);
    }
    
    gsl_fft_real_radix2_transform(tmplH, 1, N);
    gsl_fft_real_radix2_transform(tmplL, 1, N);
    whiten(tmplH, SwH, N);
    whiten(tmplL, SwL, N);

    
    // select central block
    k = M/2-N/2;
    for (i = 0; i < N; ++i) DH[i] = DHfull[i+k];
    for (i = 0; i < N; ++i) DL[i] = DLfull[i+k];
    tukey(DH, alpha, N);
    tukey(DL, alpha, N);
    gsl_fft_real_radix2_transform(DH, 1, N);
    whiten(DH, SwH, N);
    gsl_fft_real_radix2_transform(DL, 1, N);
    whiten(DL, SwL, N);
    
    
    double tH, pH, nH;
    double tL, pL, nL;
    
    Align(DH, tmplH, N, Tobs, &tH, &pH, &nH);
    Align(DL, tmplL, N, Tobs, &tL, &pL, &nL);
    
    printf("\n");
    
    printf("H1: time %e phase %e amplitude %e\n", tH, pH, nH);
    printf("L1: time %e phase %e amplitude %e\n", tL, pL, nL);
    printf("H1/L1: time %e phase %e amplitude  %e\n", tH-tL, pH-pL, nH/nL);
    printf("\n");
    
    shift(tmplH, Tobs, tH, pH, nH, N);
    shift(tmplL, Tobs, tL, pL, nL, N);
    
    x = 0.0;
    y = 0.0;
    for (i = 1; i < N/2; ++i)
    {
        x+= tmplH[i]*tmplH[i]+tmplH[N-i]*tmplH[N-i];
        y+= tmplL[i]*tmplL[i]+tmplL[N-i]*tmplL[N-i];
    }
    
    x *= 4.0;
    y *= 4.0;
    
    printf("SNR H1 %f SNR L1 %f Network %f\n", sqrt(x), sqrt(y), sqrt(x+y));
    
    
    // make a Pi/2 shifted copy of H template
    tmplX[0] = 0.0;
    tmplX[N/2] = 0.0;
    for (i = 1; i < N/2; ++i)
    {
        tmplX[i] = tmplH[N-i];
        tmplX[N-i] = -tmplH[i];
    }
    
    gsl_fft_halfcomplex_radix2_inverse(tmplX, 1, N);
    gsl_fft_halfcomplex_radix2_inverse(DH, 1, N);
    gsl_fft_halfcomplex_radix2_inverse(tmplH, 1, N);
    gsl_fft_halfcomplex_radix2_inverse(DL, 1, N);
    gsl_fft_halfcomplex_radix2_inverse(tmplL, 1, N);
    
    y = 0.0;
    out = fopen("Henv.dat","w");
    for (i = 0; i < N; ++i)
    {
        x = sqrt(tmplH[i]*tmplH[i]+tmplX[i]*tmplX[i])/var;
        fprintf(out,"%.15e %.15e\n", (double)(i-N/2)*dt, x);
        if(x > y)
        {
            y = x;
            t = (double)(i-N/2)*dt;
        }
    }
    fclose(out);
    
    printf("H max at %f\n", t);
    
    out = fopen("res.dat","w");
    for (i = 0; i < N; ++i)
    {
        fprintf(out,"%.15e %.15e %.15e %.15e %.15e\n", (double)(i-N/2)*dt, DH[i]/var, tmplH[i]/var, DL[i]/var, tmplL[i]/var);
    }
    fclose(out);
    
    
    out = fopen("data_corr_full.dat","w");
    cross_correlate(DL, DH, 2.25, 2.46, Tobs, N, out);
    fclose(out);
    
    out = fopen("data_corr_first.dat","w");
    cross_correlate(DL, DH, 2.25, 2.355, Tobs, N, out);
    fclose(out);
    
    out = fopen("data_corr_last.dat","w");
    cross_correlate(DL, DH, 2.355, 2.46, Tobs, N, out);
    fclose(out);
    
    out = fopen("data_corr_small.dat","w");
    cross_correlate(DL, DH, 2.39, 2.43, Tobs, N, out);
    fclose(out);
    
    
    for (i = 0; i < N; ++i)
    {
        DH[i] -= tmplH[i];
        DL[i] -= tmplL[i];
    }
    
    out = fopen("res_corr_full.dat","w");
    cross_correlate(DL, DH, 2.25, 2.46, Tobs, N, out);
    fclose(out);
    
    out = fopen("res_corr_first.dat","w");
    cross_correlate(DL, DH, 2.25, 2.355, Tobs, N, out);
    fclose(out);
    
    out = fopen("res_corr_last.dat","w");
    cross_correlate(DL, DH, 2.355, 2.46, Tobs, N, out);
    fclose(out);
    
    out = fopen("res_corr_small.dat","w");
    cross_correlate(DL, DH, 2.39, 2.43, Tobs, N, out);
    fclose(out);
    
    int q, p, l;
    double *cr;
    double *hst;
    
    dx = (0.43-0.39);
    
    q = (int)(dx/dt);
    p = (int)(4.0/dx);
    
    cr = (double*)malloc(sizeof(double)*(q));
    hst = (double*)malloc(sizeof(double)*(100));
    
    for (j = 0; j < 100; ++j) hst[j] = 0.0;
    
    k = 0;
    for (i = 0; i < p; ++i)
    {
        t = dx*(double)(i+1);
        if(t > 0.5 && t < 3.5)
        {
        cc(DL, DH, t-dx, t, Tobs, N, cr);
         for (j = 0; j < q; ++j)
         {
             k++;
             l = (int)(0.5*(1.0+cr[j])*100.0);
             hst[l] += 1.0;
         }
        }
    }
    
    for (j = 0; j < 100; ++j) hst[j] /= (double)(k);
    
    out = fopen("corr_dist.dat","w");
    for (j = 0; j < 100; ++j)
    {
        fprintf(out,"%f %f\n", ((double)(j)+0.5)/50.0-1.0, hst[j]);
    }
    fclose(out);
    
    
    for (i = 0; i < N; ++i)
    {
        tmplL[i] = tmplT[i];
        tmplH[i] = tmplT[i];
    }

    
    gsl_fft_real_radix2_transform(tmplH, 1, N);
    gsl_fft_real_radix2_transform(tmplL, 1, N);
    
    shift(tmplH, Tobs, tH, pH, nH, N);
    shift(tmplL, Tobs, tL, pL, nL, N);
    
    // make a Pi/2 shifted copy of H template
    tmplX[0] = 0.0;
    tmplX[N/2] = 0.0;
    for (i = 1; i < N/2; ++i)
    {
        tmplX[i] = tmplH[N-i];
        tmplX[N-i] = -tmplH[i];
    }
    
    gsl_fft_halfcomplex_radix2_inverse(tmplX, 1, N);
    gsl_fft_halfcomplex_radix2_inverse(tmplH, 1, N);
    gsl_fft_halfcomplex_radix2_inverse(tmplL, 1, N);
    
    out = fopen("templates.dat","w");
    for (i = 0; i < N; ++i)
    {
        fprintf(out,"%.15e %.15e %.15e\n", (double)(i-N/2)*dt, tmplH[i],  tmplL[i]);
    }
    fclose(out);
    
    
    out = fopen("Hanford.dat","w");
    for (i = 0; i < N; ++i)
    {
        fprintf(out,"%.15e %.15e %.15e %.15e\n", (double)(i-N/2)*dt, tmplH[i], tmplX[i],  sqrt(tmplH[i]*tmplH[i]+tmplX[i]*tmplX[i]));
    }
    fclose(out);

    
    // Now compute the specta and phase for the residuals
    
    k = M/2-N/2;
    for (i = 0; i < N; ++i) DH[i] = DHfull[i+k];
    for (i = 0; i < N; ++i) DL[i] = DLfull[i+k];
    for (i = 0; i < N; ++i) DH[i] -= tmplH[i];
    for (i = 0; i < N; ++i) DL[i] -= tmplL[i];
    
    tukey(DH, alpha, N);
    tukey(DL, alpha, N);
    
    // FFT
    gsl_fft_real_radix2_transform(DH, 1, N);
    gsl_fft_real_radix2_transform(DL, 1, N);
    
    sprintf(command, "freqres_%d_%d_%d.dat", Oflag, (int)(Tobs), tti);
    out = fopen(command,"w");
    for (i = 0; i < N/2; ++i)
    {
        x = atan2(DH[N-i],DH[i]);
        if(x < 0.0) x += TPI;
        y = atan2(DL[N-i],DL[i]);
        if(y < 0.0) y += TPI;
        fprintf(out,"%.15e %.15e %.15e %.15e %.15e\n", (double)(i)/Tobs, fac*2.0*(DH[i]*DH[i]+DH[N-i]*DH[N-i]), x, fac*2.0*(DL[i]*DL[i]+DL[N-i]*DL[N-i]), y);
    }
    fclose(out);
    
    
    // double check that the subtracrtion was done properly
    
    whiten(DH, SwH, N);
    whiten(DL, SwL, N);
    
    gsl_fft_halfcomplex_radix2_inverse(DH, 1, N);
    gsl_fft_halfcomplex_radix2_inverse(DL, 1, N);
    
    out = fopen("check.dat","w");
    for (i = 0; i < N; ++i)
    {
        fprintf(out,"%.15e %.15e %.15e\n", (double)(i-N/2)*dt, DH[i],  DL[i]);
    }
    fclose(out);
    
    
    // Compute spectrum using unfiltered and un-windowed resdiual
    
    k = M/2-N/2;
    for (i = 0; i < N; ++i) DH[i] = DHcopy[i+k];
    for (i = 0; i < N; ++i) DL[i] = DLcopy[i+k];
    // subtract ML template
    for (i = 0; i < N; ++i) DH[i] -= tmplH[i];
    for (i = 0; i < N; ++i) DL[i] -= tmplL[i];
    
    // FFT
    gsl_fft_real_radix2_transform(DH, 1, N);
    gsl_fft_real_radix2_transform(DL, 1, N);
    
    sprintf(command, "freqnowindowres_%d_%d_%d.dat", Oflag, (int)(Tobs), tti);
    out = fopen(command,"w");
    for (i = 0; i < N/2; ++i)
    {
        x = atan2(DH[N-i],DH[i]);
        if(x < 0.0) x += TPI;
        y = atan2(DL[N-i],DL[i]);
        if(y < 0.0) y += TPI;
        fprintf(out,"%.15e %.15e %.15e %.15e %.15e\n", (double)(i)/Tobs, fac*2.0*(DH[i]*DH[i]+DH[N-i]*DH[N-i]), x, fac*2.0*(DL[i]*DL[i]+DL[N-i]*DL[N-i]), y);
    }
    fclose(out);

    
    
    free(DH);
    free(DL);
    free(DHfull);
    free(DLfull);
    free(DHcopy);
    free(DLcopy);
    
    free(SnH);
    free(SwH);
    free(SnL);
    free(SwL);
    
    free(tmplT);
    free(tmplH);
    free(tmplL);
    free(tmplX);
    
    
    return 0;

}

void cc(double *Hf, double *Lf, double t1, double t2, double T, int M,  double *cr)
{
    int i, j, k, q;
    int istart, istop;
    double dt, t;
    double Hnorm, Lnorm, Corr;
    
    
    dt = T/(double)(M);
    
    q = (int)((t2-t1)/(2.0*dt));
    
    istart = (int)(t1/dt);
    istop = (int)(t2/dt);
    
    Hnorm = 0.0;
    Lnorm = 0.0;
    for(i=istart; i<=istop; i++)
    {
        Hnorm += Hf[i]*Hf[i];
        Lnorm += Lf[i]*Lf[i];
    }
    
    
    for(j=-q; j<=q; j++)
    {
        Corr = 0.0;
        for(i=istart; i<=istop; i++)
        {
            // using a periodic interval
            k = i+j;
            if(k < istart) k = istop - (istart - k);
            if(k > istop) k = istart + (k-istop);
            Corr += Hf[i]*Lf[k];
        }
        Corr /= sqrt(Hnorm*Lnorm);
        cr[j+q] =  Corr;
    }
    
    
}



void cross_correlate(double *Hf, double *Lf, double t1, double t2, double T, int M,  FILE *out)
{
    int i, j, k, q;
    int istart, istop;
    double dt, t;
    double Hnorm, Lnorm, Corr;

    
    dt = T/(double)(M);
    
    q = (int)((t2-t1)/(2.0*dt));
    
    istart = (int)(t1/dt);
    istop = (int)(t2/dt);
    
    Hnorm = 0.0;
    Lnorm = 0.0;
    for(i=istart; i<=istop; i++)
    {
        Hnorm += Hf[i]*Hf[i];
        Lnorm += Lf[i]*Lf[i];
    }
    
    
    for(j=-q; j<=q; j++)
    {
        Corr = 0.0;
        for(i=istart; i<=istop; i++)
        {
            // using a periodic interval
            k = i+j;
            if(k < istart) k = istop - (istart - k);
            if(k > istop) k = istart + (k-istop);
            Corr += Hf[i]*Lf[k];
        }
        Corr /= sqrt(Hnorm*Lnorm);
        t = dt*(double)(j);
        fprintf(out,"%.14e %e\n", t, Corr);
    }
    
    
    
    
}


// a is data, b is the template

void Align(double *a, double *b, int n, double T, double *tx, double *px, double *nx)
{
    double max, min;
    int i, j, k;
    int index;
    double tmp, mx, SNRsq, hsq, norm;
    double cn, sn, Phase, freq, df, ts, ps;
    double hpr, hpi, x;
    double *AC, *AF, *Sn;
    double *corr;
    double psh;
    
    int imin, imax;
    
    AC=(double*)malloc(sizeof(double)*(n));
    AF=(double*)malloc(sizeof(double)*(n));
    corr = (double*)malloc(sizeof(double)*(n));
    
    df = 1./(T);
    
    for(i = 0; i < n; i++) corr[i] = 0.0;
    
    phase_blind_time_shift(AC, AF, a, b, n);
    
    for(i = 0; i < n; i++) corr[i] += sqrt(AC[i]*AC[i]+AF[i]*AF[i]);
    
    max_array_element(&max, &index, corr, n);
    
    SNRsq = corr[index];
    
    hsq = 0.0;
    for(i = 0; i < n/2; i++) hsq += 2.0*(b[i]*b[i]+b[n-i]*b[n-i]);
    
    norm = SNRsq/hsq;
    
    ps = atan2(AF[index],AC[index]);
    
    if(index < (n/2)-1)
    ts = ((double) index)/((double) n)*T;
    else if(index >= n/2)
    ts = ((double) (index - n))/((double) n)*T;
    
    *tx = -ts;
    *px = -ps;
    *nx = norm;
    
    free(AC);
    free(AF);
    free(corr);
    
    return;
}

void shift(double *b, double T, double ts, double ps, double norm, int n)
{
    
    int i;
    double Phase, freq, cn, sn, hpr, hpi, df;
    
    df = 1.0/T;
    

    for(i=1; i< n/2; i++)
    {
        
        freq = (double)(i)*df;
        
        Phase = TPI*freq*ts - ps;
        
        
        cn = cos(Phase);
        sn = sin(Phase);
        
        hpr = (b[i]*cn - b[n-i]*sn);
        hpi = (b[i]*sn + b[n-i]*cn);
        
        b[i] = hpr*norm;
        b[n-i] = hpi*norm;
    }
    
    b[0] = 0.0;
    b[n/2-1] = 0.0;
    
    
}




void bwbpf(double *in, double *out, int fwrv, int M, int n, double s, double f1, double f2)
{
    /* Code from http://www.exstrom.com/journal/sigproc/ */
    
    /*
     *                            COPYRIGHT
     *
     *  Copyright (C) 2014 Exstrom Laboratories LLC
     *
     *  This program is free software; you can redistribute it and/or modify
     *  it under the terms of the GNU General Public License as published by
     *  the Free Software Foundation; either version 2 of the License, or
     *  (at your option) any later version.
     *
     *  This program is distributed in the hope that it will be useful,
     *  but WITHOUT ANY WARRANTY; without even the implied warranty of
     *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     *  GNU General Public License for more details.
     *
     *  A copy of the GNU General Public License is available on the internet at:
     *  http://www.gnu.org/copyleft/gpl.html
     *
     *  or you can write to:
     *
     *  The Free Software Foundation, Inc.
     *  675 Mass Ave
     *  Cambridge, MA 02139, USA
     *
     *  Exstrom Laboratories LLC contact:
     *  stefan(AT)exstrom.com
     *
     *  Exstrom Laboratories LLC
     *  Longmont, CO 80503, USA
     *
     */
    
    /* Butterworth bandpass filter
     n = filter order 4,8,12,...
     s = sampling frequency
     f1 = upper half power frequency
     f2 = lower half power frequency  */
    
    if(n % 4){ printf("Order must be 4,8,12,16,...\n"); return;}
    
    int i, j;
    double a = cos(PI*(f1+f2)/s)/cos(PI*(f1-f2)/s);
    double a2 = a*a;
    double b = tan(PI*(f1-f2)/s);
    double b2 = b*b;
    double r;
    
    n = n/4;
    double *A = (double *)malloc(n*sizeof(double));
    double *d1 = (double *)malloc(n*sizeof(double));
    double *d2 = (double *)malloc(n*sizeof(double));
    double *d3 = (double *)malloc(n*sizeof(double));
    double *d4 = (double *)malloc(n*sizeof(double));
    double *w0 = (double *)malloc(n*sizeof(double));
    double *w1 = (double *)malloc(n*sizeof(double));
    double *w2 = (double *)malloc(n*sizeof(double));
    double *w3 = (double *)malloc(n*sizeof(double));
    double *w4 = (double *)malloc(n*sizeof(double));
    double x;
    
    for(i=0; i<n; ++i)
    {
        r = sin(PI*(2.0*(double)i+1.0)/(4.0*(double)n));
        s = b2 + 2.0*b*r + 1.0;
        A[i] = b2/s;
        d1[i] = 4.0*a*(1.0+b*r)/s;
        d2[i] = 2.0*(b2-2.0*a2-1.0)/s;
        d3[i] = 4.0*a*(1.0-b*r)/s;
        d4[i] = -(b2 - 2.0*b*r + 1.0)/s;
        w0[i] = 0.0;
        w1[i] = 0.0;
        w2[i] = 0.0;
        w3[i] = 0.0;
        w4[i] = 0.0;
    }
    
    for(j=0; j< M; ++j)
    {
        if(fwrv == 1) x = in[j];
        if(fwrv == -1) x = in[M-j-1];
        for(i=0; i<n; ++i)
        {
            w0[i] = d1[i]*w1[i] + d2[i]*w2[i]+ d3[i]*w3[i]+ d4[i]*w4[i] + x;
            x = A[i]*(w0[i] - 2.0*w2[i] + w4[i]);
            w4[i] = w3[i];
            w3[i] = w2[i];
            w2[i] = w1[i];
            w1[i] = w0[i];
        }
        if(fwrv == 1) out[j] = x;
        if(fwrv == -1) out[M-j-1] = x;
    }
    
    free(A);
    free(d1);
    free(d2);
    free(d3);
    free(d4);
    free(w0);
    free(w1);
    free(w2);
    free(w3);
    free(w4);
    
    return;
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

void tukey_scale(double *s1, double *s2, double alpha, int N)
{
    int i, imin, imax;
    double x1, x2;
    double filter;
    
    imin = (int)(alpha*(double)(N-1)/2.0);
    imax = (int)((double)(N-1)*(1.0-alpha/2.0));
    
    x1 = 0.0;
    x2 = 0.0;
    for(i=0; i< N; i++)
    {
        filter = 1.0;
        if(i < imin) filter = 0.5*(1.0+cos(PI*( (double)(i)/(double)(imin)-1.0 )));
        if(i>imax) filter = 0.5*(1.0+cos(PI*( (double)(i)/(double)(imin)-2.0/alpha+1.0 )));
        x1 += filter;
        x2 += filter*filter;
    }
    x1 /= (double)(N);
    x2 /= (double)(N);
    
    *s1 = x1;
    *s2 = x2;
    
}





void max_array_element(double *max, int *index, double *array, int n)
{
    int i;
    
    *max = array[0];
    *index = 0;
    
    for(i = 1; i <= n-1; i++)
    {
        if(array[i] > *max)
        {
            *max = array[i];
            *index = i;
        }
    }
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


void phase_blind_time_shift(double *corr, double *corrf, double *data1, double *data2, int n)
{
    int nb2, i, l, k, j;
    double scale;
    
    nb2 = n / 2;
    
    scale = (double)(n);
    
    for (i=1; i < nb2; i++)
    {
        l=i;
        k=n-i;
        corr[l]	= scale*(data1[l]*data2[l] + data1[k]*data2[k]);
        corr[k]	= scale*(data1[k]*data2[l] - data1[l]*data2[k]);
        corrf[l] = corr[k];
        corrf[k] = -corr[l];
    }
    
    corr[0] = 0.0;
    corrf[0] = 0.0;
    corr[nb2] = 0.0;
    corrf[nb2] = 0.0;
    
    gsl_fft_halfcomplex_radix2_inverse(corr, 1, n);
    gsl_fft_halfcomplex_radix2_inverse(corrf, 1, n);
    
    
}


void SineGaussianF(double *hs, double *sigpar, double Tobs, int NMAX)
{
    double f0, t0, Q, sf, sx, Amp;
    double fmax, fmin, fac;
    double phi, f, t, x, y, z, zold;
    double tau;
    int imin, imax;
    
    int i, id, N;
    
    t0 = sigpar[0];
    f0 = sigpar[1];
    Q = sigpar[2];
    Amp = sigpar[3];
    phi = sigpar[4];
    
    tau = Q/(TPI*f0);
    
    fmax = f0 + 3.0/tau;  // no point evaluating waveform past this time (many efolds down)
    fmin = f0 - 3.0/tau;  // no point evaluating waveform before this time (many efolds down)
    
    fac = sqrt(Tobs);
    
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
            sf = (Amp/fac)*RTPI/2.0*tau*exp(-PI*PI*tau*tau*(f-f0)*(f-f0));
            sx = exp(-Q*Q*f/f0);
            hs[i] = sf*(cos(TPI*f*t0-phi)+sx*cos(TPI*f*t0+phi));
            hs[NMAX-i] = sf*(sin(TPI*f*t0-phi)+sx*sin(TPI*f*t0+phi));
        }
    }
    
    
}

void SineGaussianC(double *hs, double *sigpar, double Tobs, int NMAX)
{
    double f0, t0, Q, sf, sx, Amp;
    double fmax, fmin, fac;
    double phi, f, t, x, y, z, zold;
    double tau, dt;
    int imin, imax;
    
    int i, id, N;
    
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

void pshift(double *a, double R, double delt, double pshift, int n, double Tobs)
{
    int i, j, k;
    double  f;
    double ReA, ImA;
    double cx, sx;
    
    for(i=0; i<n/2; i++)
    {
        f = (double)(i)/Tobs;
        j = i; // real
        k = n-i; // imaginary
        cx = cos(-pshift+TPI*delt*f);
        sx = sin(-pshift+TPI*delt*f);
        ReA = a[j]*cx-a[k]*sx;
        ImA = a[k]*cx+a[j]*sx;
        a[j] = R*ReA;
        a[k] = R*ImA;
    }
    
    return;
    
}


double fourier_nwip_shift(double *a, double *b, double delt, double pshift, int n, double Tobs, int imin, int imax)
{
    int i, j, k;
    double arg, product, f;
    double ReA, ReB, ImA, ImB;
    double cx, sx;
    
    // Does f_nwip with a given time and phase shift
    
    arg = 0.0;
    for(i=1; i<n/2; i++)
    {
        if(i > imin && i < imax)
        {
            f = (double)(i)/Tobs;
            j = i; // real
            k = n-i; // imaginary
            ReA = a[j];
            ImA = a[k];
            cx = cos(-pshift+TPI*delt*f);
            sx = sin(-pshift+TPI*delt*f);
            ReB = b[j]*cx-b[k]*sx;
            ImB = b[k]*cx+b[j]*sx;
            product = ReA*ReB + ImA*ImB;
            arg += product;
        }
    }
    
    return(arg);
    
}

