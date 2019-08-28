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
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_cdf.h>

#define TPI 6.2831853071795862319959269370884     // 2 Pi
#define SQPI 2.5066282746310002  // sqrt(TPI)

// gcc -o Gauss_Test_Freq Gauss_Test_Freq.c -lm -lgsl


int main(int argc, char *argv[])
{

    char filename[1024];
    
    FILE *dFile;
    FILE *nFile;
    FILE *out;
    
    int N, M, n,nn,i,kx,cnt;
    double f, t, dre,dim,Sn,Sn25,Sn75,Sn05,Sn95;
    double mean, var, std, dt, df, Tx, T, dx, x, y, z, u;
    double lx, ly;
    double fmax, fmin;
    int imax, imin;
    
    int Nsamp, Nsec;
    int ii, jj;
    double *samples=NULL;
    double *samps=NULL;
    double *times=NULL;
    double *freqs=NULL;
    
    double *hist = malloc (400 * sizeof (double));
    
 
    
    if (argc < 2) {
        
        printf("Usage: %s filename\n", argv[0]);
        return 1;
    }

    
    // how many samples in the data?
    N=0;
    
    //get data file
    dFile = fopen(argv[1],"r");
    
    fmin = 32;
    fmax = 500;
    
    Nsamp = 0;
    
    while(!feof(dFile))
    {
        fscanf(dFile,"%lg %lg %lg", &f, &x, &y);
        N++;
        if(f >= fmin && f <= fmax) Nsamp += 2;
    }
    N--;
    
    rewind(dFile);
    
    // get memory to ingest file
    samples=malloc(Nsamp*sizeof(double));
    freqs=malloc(N*sizeof(double));
    

    out=fopen("freq_RI.dat","w");
    i = 0;
    for(n=0; n<N; n++)
    {
        fscanf(dFile,"%lg %lg %lg", &f, &x, &y);
        freqs[n] = f;
        if(f >= fmin && f <= fmax)
        {
        samples[2*i] = 2.0*x;
        samples[2*i+1] = 2.0*y;
        fprintf(out,"%e %e %e\n", f, samples[2*i], samples[2*i+1]);
        i++;
        }
     }
    fclose(dFile);
    fclose(out);
    
    df = freqs[1]-freqs[0];
    T = 1.0/df;
    
    printf("%f\n", T);
    
    printf("%d %d\n", N, Nsamp);
    
    samps=malloc(Nsamp*sizeof(double));
    
    
    mean = 0.0;
    var = 0.0;
    
    for(n=0; n<Nsamp; n++)
    {
        mean += samples[n];
        var += samples[n]*samples[n];
    }
        
    mean /= (double)(Nsamp);
    var /= (double)(Nsamp);
    var -= mean*mean;
    std = sqrt(var);
    
    printf("Mean %f  Variance %f\n", mean, var);
    
    for(n=0; n<Nsamp; n++)
    {
        samps[n] = samples[n];
    }
    
        
        for(i=0; i< 400; i++) hist[i] = 0.0;
        for(n=0; n<Nsamp; n++)
        {
            kx = (int)(((samps[n]+8.0)/16.0)*400.0);
            if(kx > 0 && kx < 400) hist[kx] += 1.0;
        }
    
       printf("%f\n", hist[200]);
    
        for(i=0; i< 400; i++) hist[i] /= (double)(Nsamp);
    
        dx = (16.0/400.0);
    
        y = 1.0/sqrt(TPI);
    
        x = -8.0+(((double)(200)+0.5)/400.0)*16.0;
    
        printf("%f\n", y*exp(-x*x/2.0)*(double)(Nsamp));
    
        z = 0.0;
        u = 0.0;
        out=fopen("hist_freq.dat","w");
        for(i=0; i< 400; i++)
        {
            x = -8.0+(((double)(i)+0.5)/400.0)*16.0;
            fprintf(out,"%e %e %e %e\n", x, hist[i]/dx, y*exp(-x*x/2.0), dx*y*exp(-x*x/2.0)*(double)(Nsamp));
            z += hist[i];
            u += y*exp(-x*x/2.0)*dx;
        }
        fclose(out);
    
     printf("%f %f\n", z, u);
    
        
    
    gsl_sort(samps,1,Nsamp);
    
    
    double S,A;
    
    // Anderson-Darling statistics
    S=0.0;
    for(n=0; n<Nsamp; n++)
    {
 
        x = gsl_cdf_ugaussian_P(samps[n]);
        if(x < 0.999)
        {
            lx = log(x);
            ly = log(1.0-x);
        }
        else
        {
            u = samps[n];
            z = 1.0-1.0/(u*u)+3.0/(u*u*u*u)-15.0/(u*u*u*u*u*u)+105.0/(u*u*u*u*u*u*u*u);
            y = z*exp(-u*u/2.0)/(u*SQPI);
            lx = -y;
            ly = -u*u/2.0 +log(z/(u*SQPI));
            //printf("%f %e %e %e %e\n", u, log(x), log(1.0-x), lx, ly);
        }
        
        
         S += ((2.*(double)(n+1)-1.0)*lx + ((double)(2*Nsamp)-2.*(double)(n+1)+1.0)*ly)/(double)(Nsamp);
    }
    
    
    
    A = (double)(-1*Nsamp) - S;
    
    A *= (1.0+0.75/(double)(Nsamp)+2.25/(double)(Nsamp*Nsamp));
    
    
    printf("%lg %lg\n", S, A);
    
    
    free(samps);
    free(samples);
    free(times);
    free(freqs);


    
    return 0;
    
}



