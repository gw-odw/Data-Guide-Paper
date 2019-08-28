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
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_spline.h>

#define PI 3.1415926535897932

// gcc -o fig1_corr fig1_corr.c -lm -lgsl

void cross_correlate(double *Hf, double *Lf, double t1, double t2, double T, int M,  int *q, double *corr);
void bwbpf(double *in, double *out, int fwrv, int M, int n, double s, double f1, double f2);

int main(int argc, char *argv[])
{
    int i, ii, jj, q, j, N;
    double T, x, y, z, dt, t;
    double *H, *L, *time, *corr;
    double fmax, fmin;
    char command[1024];
    
    FILE *in;
    FILE *out;
    
    /* Needs the GWOSC files 
     https://www.gw-openscience.org/GW150914data/P150914/fig1-observed-H.txt
     https://www.gw-openscience.org/GW150914data/P150914/fig1-observed-L.txt
     https://www.gw-openscience.org/GW150914data/P150914/fig1-residual-H.txt
     https://www.gw-openscience.org/GW150914data/P150914/fig1-residual-L.txt */
    
    N = 3441;
    H = (double*)malloc(sizeof(double)*(N));
    L = (double*)malloc(sizeof(double)*(N));
    corr = (double*)malloc(sizeof(double)*(N));
    time = (double*)malloc(sizeof(double)*(N));
    
    in = fopen("fig1-residual-H.txt", "r");
    // strip the header
    fgets(command, 1024, in);
    for (i = 0; i < N; ++i)
    {
        fscanf(in,"%lf%lf", &time[i], &H[i]);
    }
    fclose(in);
    
    in = fopen("fig1-residual-L.txt", "r");
    // strip the header
    fgets(command, 1024, in);
    for (i = 0; i < N; ++i)
    {
        fscanf(in,"%lf%lf", &time[i], &L[i]);
    }
    fclose(in);
    
    dt = time[1]-time[0];
    T = dt*(double)N;
    
    cross_correlate(L, H, 0.0, T, T, N, &q, corr);
    out = fopen("fig1_full_corr.dat","w");
    for(j=-q; j<=q; j++)
    {
        t = dt*(double)(j);
        fprintf(out,"%.14e %e\n", t*1000.0, corr[j+q]);
    }
    fclose(out);
    
    
    out = fopen("fig1_first_corr.dat","w");
    cross_correlate(L, H, 0.0, T/2.0, T, N, &q, corr);
    for(j=-q; j<=q; j++)
    {
        t = dt*(double)(j);
        fprintf(out,"%.14e %e\n", t*1000.0, corr[j+q]);
    }
    fclose(out);
    
    out = fopen("fig1_last_corr.dat","w");
    cross_correlate(L, H, T/2.0, T, T, N, &q, corr);
    for(j=-q; j<=q; j++)
    {
        t = dt*(double)(j);
        fprintf(out,"%.14e %e\n", t*1000.0, corr[j+q]);
    }
    fclose(out);
    
    // record starts at 0.25, so 0.39 to 0.43 is
    
    out = fopen("fig1_small_corr.dat","w");
    cross_correlate(L, H, 0.14, 0.18, T, N, &q, corr);
    for(j=-q; j<=q; j++)
    {
        t = dt*(double)(j);
        fprintf(out,"%.14e %e\n", t*1000.0, corr[j+q]);
        
    }
    fclose(out);
    
    double dx;
    int p, l, k;
    double *hst;
    
    dx = 0.18-0.14;
    p = (int)(T/dx);
    printf("# short segements %d\n", p);
    hst = (double*)malloc(sizeof(double)*100);
    
    for (j = 0; j < 100; ++j) hst[j] = 0.0;
    
    k = 0;
    for (i = 0; i < p; ++i)
    {
            t = dx*(double)(i+1);
            cross_correlate(L, H, t-dx, t, T, N, &q, corr);
            //printf("%d %d\n", i, q);
            for (j = -q; j <= q; ++j)
            {
                k++;
                l = (int)(0.5*(1.0+corr[j+q])*100.0);
                //printf("%d %d\n", i, j);
                hst[l] += 1.0;
            }
    }
    
    for (j = 0; j < 100; ++j) hst[j] /= (double)(k);
    
    out = fopen("corr_dist_fig1.dat","w");
    for (j = 0; j < 100; ++j)
    {
        fprintf(out,"%f %f\n", ((double)(j)+0.5)/50.0-1.0, hst[j]);
    }
    fclose(out);
    
    
    
    in = fopen("fig1-observed-H.txt", "r");
    // strip the header
    fgets(command, 1024, in);
    for (i = 0; i < N; ++i)
    {
        fscanf(in,"%lf%lf", &time[i], &H[i]);
    }
    fclose(in);
    
    in = fopen("fig1-observed-L.txt", "r");
    // strip the header
    fgets(command, 1024, in);
    for (i = 0; i < N; ++i)
    {
        fscanf(in,"%lf%lf", &time[i], &L[i]);
    }
    fclose(in);
    
    dt = time[1]-time[0];
    T = dt*(double)N;
    
    cross_correlate(L, H, 0.0, T, T, N, &q, corr);
    out = fopen("fig1d_full_corr.dat","w");
    for(j=-q; j<=q; j++)
    {
        t = dt*(double)(j);
        fprintf(out,"%.14e %e\n", t*1000.0, corr[j+q]);
    }
    fclose(out);
    
    
    out = fopen("fig1d_first_corr.dat","w");
    cross_correlate(L, H, 0.0, T/2.0, T, N, &q, corr);
    for(j=-q; j<=q; j++)
    {
        t = dt*(double)(j);
        fprintf(out,"%.14e %e\n", t*1000.0, corr[j+q]);
    }
    fclose(out);
    
    out = fopen("fig1d_last_corr.dat","w");
    cross_correlate(L, H, T/2.0, T, T, N, &q, corr);
    for(j=-q; j<=q; j++)
    {
        t = dt*(double)(j);
        fprintf(out,"%.14e %e\n", t*1000.0, corr[j+q]);
    }
    fclose(out);
    
    // record starts at 0.25, so 0.39 to 0.43 is
    
    out = fopen("fig1d_small_corr.dat","w");
    cross_correlate(L, H, 0.14, 0.18, T, N, &q, corr);
    for(j=-q; j<=q; j++)
    {
        t = dt*(double)(j);
        fprintf(out,"%.14e %e\n", t*1000.0, corr[j+q]);
    }
    fclose(out);

    
    
    
    
    const gsl_rng_type * P;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    
    P = gsl_rng_default;
    r = gsl_rng_alloc (P);
    
    
    for (i = 0; i < N; ++i)
    {
        H[i] = gsl_ran_gaussian(r,1.0);
        L[i] = gsl_ran_gaussian(r,1.0);
    }
    
    out = fopen("rand.dat","w");
    for (i = 0; i < N; ++i)
    {
        fprintf(out,"%.14e %e %e\n", dt*(double)(i), H[i], L[i]);
    }
    
    fmax = 350.0;
    fmin = 35.0;
    
    // The filter decays in time with an efold scale of 1/fmin. Need to throw out several/fmin
    // zero-phase Butterworth band-pass filter (zero phase by going forward/back)
    bwbpf(H, H, 1, N, 8, 1.0/dt, fmax, fmin);
    bwbpf(H, H, -1, N, 8, 1.0/dt, fmax, fmin);
    bwbpf(L, L, 1, N, 8, 1.0/dt, fmax, fmin);
    bwbpf(L, L, -1, N, 8, 1.0/dt, fmax, fmin);
    
    out = fopen("rand_dp.dat","w");
    for (i = 0; i < N; ++i)
    {
        fprintf(out,"%.14e %e %e\n", dt*(double)(i), H[i], L[i]);
    }
    
    double av, var;
    
    
    out = fopen("rand_small_corr.dat","w");
    cross_correlate(L, H, 0.14, 0.18, T, N, &q, corr);
    for(j=-q; j<=q; j++)
    {
        t = dt*(double)(j);
        fprintf(out,"%.14e %e\n", t*1000.0, corr[j+q]);
    }
    fclose(out);
    
    double *hist;
    int bin;
    
    
    hist = (double*)malloc(sizeof(double)*(100));
    for (i = 0; i < 100; ++i) hist[i] = 0.0;
    
    
    jj = 0.0;
    av = 0.0;
    var = 0.0;
    
    for (ii = 0; ii < 10000; ++ii)
    {
    
    for (i = 0; i < N; ++i)
    {
        H[i] = gsl_ran_gaussian(r,1.0);
        L[i] = gsl_ran_gaussian(r,1.0);
    }
    bwbpf(H, H, 1, N, 8, 1.0/dt, fmax, fmin);
    bwbpf(H, H, -1, N, 8, 1.0/dt, fmax, fmin);
    bwbpf(L, L, 1, N, 8, 1.0/dt, fmax, fmin);
    bwbpf(L, L, -1, N, 8, 1.0/dt, fmax, fmin);
    cross_correlate(L, H, 0.14, 0.18, T, N, &q, corr);
        
        for(j=-q; j<=q; j++)
        {
            bin = (int)((1.0+corr[j+q])*50.0);
            hist[bin] += 1.0;
            av += corr[j+q];
            var += corr[j+q]*corr[j+q];
            jj++;
        }
        
    }
    
    for (i = 0; i < 100; ++i) hist[i] /= (double)(jj);
    

    
    
    av /= (double)(jj);
    var /= (double)(jj);
    
    printf("Small segment: mean %e variance %e deviation %e\n", av, var-av*av, sqrt(var-av*av));
    
    x = sqrt(1.0/(2.0*PI*(var-av*av)))/50.0;
    y = 2.0*(var-av*av);
    
    out = fopen("small_corr.dat","w");
    for(j=0; j< 100; j++)
    {
        z = -1.0+((double)(j)+0.5)/50.0;
        fprintf(out,"%.14e %e %e\n", z, hist[j], x*exp(-z*z/y));
    }
    fclose(out);
    
    
    jj = 0.0;
    av = 0.0;
    var = 0.0;
    
    for (ii = 0; ii < 100; ++ii)
    {
        
        for (i = 0; i < N; ++i)
        {
            H[i] = gsl_ran_gaussian(r,1.0);
            L[i] = gsl_ran_gaussian(r,1.0);
        }
        bwbpf(H, H, 1, N, 8, 1.0/dt, fmax, fmin);
        bwbpf(H, H, -1, N, 8, 1.0/dt, fmax, fmin);
        bwbpf(L, L, 1, N, 8, 1.0/dt, fmax, fmin);
        bwbpf(L, L, -1, N, 8, 1.0/dt, fmax, fmin);
        cross_correlate(L, H, T/4.0, 3.0*T/4.0, T, N, &q, corr);
        
        for(j=-q; j<=q; j++)
        {
            av += corr[j+q];
            var += corr[j+q]*corr[j+q];
            jj++;
        }
        
    }
    
    av /= (double)(jj);
    var /= (double)(jj);
    
    printf("Half segments: mean %e variance %e deviation %e\n", av, var-av*av, sqrt(var-av*av));

    jj = 0.0;
    av = 0.0;
    var = 0.0;
    
    for (ii = 0; ii < 100; ++ii)
    {
        
        for (i = 0; i < N; ++i)
        {
            H[i] = gsl_ran_gaussian(r,1.0);
            L[i] = gsl_ran_gaussian(r,1.0);
        }
        bwbpf(H, H, 1, N, 8, 1.0/dt, fmax, fmin);
        bwbpf(H, H, -1, N, 8, 1.0/dt, fmax, fmin);
        bwbpf(L, L, 1, N, 8, 1.0/dt, fmax, fmin);
        bwbpf(L, L, -1, N, 8, 1.0/dt, fmax, fmin);
        cross_correlate(L, H, 0.0, T, T, N, &q, corr);
        
        for(j=-q; j<=q; j++)
        {
            av += corr[j+q];
            var += corr[j+q]*corr[j+q];
            jj++;
        }
        
    }
    
    av /= (double)(jj);
    var /= (double)(jj);
    
    printf("full segment: mean %e variance %e deviation %e\n", av, var-av*av, sqrt(var-av*av));
    
    
    free(H);
    free(L);
    free(corr);
    free(time);
    
    
    
    
    return 0;
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



void cross_correlate(double *Hf, double *Lf, double t1, double t2, double T, int M,  int *p, double *cr)
{
    int i, j, k, q;
    int istart, istop;
    double dt, t;
    double Hnorm, Lnorm, Corr;
    
    
    dt = T/(double)(M);
    
    q = (int)((t2-t1)/(2.0*dt));
    
    *p = q;
    
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
        cr[q+j] = Corr;
    }
    
    
}
