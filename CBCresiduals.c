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

// gcc -o CBCresiduals CBCresiduals.c -lm -lgsl

/*  prototypes */
double adinf(double z);
double errfix(int n,double x);
double AD(int n,double z);


int main()
{
  int i, Nsamp;
  double SNR;
  double junk;
  double h1t, l1t, f, HR, HI, LR, LI, x, y, dx, te;
  char filename[1024];
  int n;
  double *histH, *histL;
  double *SH, *SL;
  double *quantH, *quantL;
  double *qb, *qn;
  double *H, *L;
  FILE *h1;
  FILE *l1;
  FILE *time;
  FILE *out;
    
    histH = malloc (100 * sizeof (double));
    histL = malloc (100 * sizeof (double));
    quantH = malloc (100 * sizeof (double));
    quantL = malloc (100 * sizeof (double));
    qb = malloc (100 * sizeof (double));
    qn = malloc (100 * sizeof (double));
    
    h1 = fopen("clean_frequency_residual_199.dat.0","r");
    l1 = fopen("clean_frequency_residual_199.dat.1","r");

	out = fopen("CBC_resdiuals.dat","w");
    
    x = sqrt(2.0);
    
    for(n=0; n< 100; n++)
    {
        histH[n] = 0.0;
        histL[n] = 0.0;
        qb[n] = -5.0+10.0*(double)(n)/100.0;  // quantile boundaries
    }
    
    Nsamp = 2*1920;
    
    H=malloc(Nsamp*sizeof(double));
    L=malloc(Nsamp*sizeof(double));
    SH=malloc(Nsamp*sizeof(double));
    SL=malloc(Nsamp*sizeof(double));
  
    for(n=0; n< 1920; n++)
    {
        fscanf(h1,"%lf%lf%lf%lf", &f, &HR, &HI, &SH[n]);
        fscanf(l1,"%lf%lf%lf%lf", &f, &LR, &LI, &SL[n]);
        HR *= x;
        HI *= x;
        LR *= x;
        LI *= x;
        H[2*n] = HR;
        H[2*n+1] = HI;
        L[2*n] = LR;
        L[2*n+1] = LI;
        fprintf(out,"%e %e %e %e %e\n", f, HR, HI, LR, LI);
        i = (int)(((HR+5.0)/10.0)*100.0);
        if(i > -1 && i < 100) histH[i] += 1.0;
        i = (int)(((HI+5.0)/10.0)*100.0);
        if(i > -1 && i < 100) histH[i] += 1.0;
        i = (int)(((LR+5.0)/10.0)*100.0);
        if(i > -1 && i < 100) histL[i] += 1.0;
        i = (int)(((LI+5.0)/10.0)*100.0);
        if(i > -1 && i < 100) histL[i] += 1.0;
    }
    
    fclose(h1);
    fclose(l1);
    fclose(out);
    
    gsl_sort(H,1,Nsamp);
    gsl_sort(L,1,Nsamp);
    
    // using the "every kth" method, otherwise too many points
    out=fopen("PP.dat","w");
    for(n=0; n< Nsamp; n++)
    {
        if(n%10==0) fprintf(out,"%f %f %f\n",  gsl_cdf_ugaussian_Pinv((double)(n)/(double)(Nsamp)),  H[n], L[n]);
    }
    fclose(out);
    
    out=fopen("PPref.dat","w");
    fprintf(out,"%f %f\n",  -4.0, -4.0);
    fprintf(out,"%f %f\n",  4.0, 4.0);
    fclose(out);
    
    
    
    
    
    for(n=0; n< 100; n++)
    {
        histH[n] /= (double)(Nsamp);
        histL[n] /= (double)(Nsamp);
    }
    
    
    dx = (10.0/100.0);
    
    y = 1.0/sqrt(TPI);
    
   
    out=fopen("hist_freq_H.dat","w");
    for(i=0; i< 100; i++)
    {
        x = -5.0+(((double)(i)+0.5)/100.0)*10.0;
        fprintf(out,"%e %e %e %e\n", x, histH[i]/dx, y*exp(-x*x/2.0), dx*y*exp(-x*x/2.0)*(double)(Nsamp));
    }
    fclose(out);
    
    out=fopen("hist_freq_L.dat","w");
    for(i=0; i< 100; i++)
    {
        x = -5.0+(((double)(i)+0.5)/100.0)*10.0;
        fprintf(out,"%e %e %e %e\n", x, histL[i]/dx, y*exp(-x*x/2.0), dx*y*exp(-x*x/2.0)*(double)(Nsamp));
    }
    fclose(out);
    
    double mean, var, std;
    
    mean = 0.0;
    var = 0.0;
    
    for(n=0; n<Nsamp; n++)
    {
        mean += H[n];
        var += H[n]*H[n];
    }
    
    mean /= (double)(Nsamp);
    var /= (double)(Nsamp);
    var -= mean*mean;
    std = sqrt(var);
    
    printf("Hanford Mean %f  Variance %f\n", mean, var);
    

    
    
    mean = 0.0;
    var = 0.0;
    
    for(n=0; n<Nsamp; n++)
    {
        mean += L[n];
        var += L[n]*L[n];
    }
    
    mean /= (double)(Nsamp);
    var /= (double)(Nsamp);
    var -= mean*mean;
    std = sqrt(var);
    
    printf("Livingston Mean %f  Variance %f\n", mean, var);
    
    
    double A;
    double lx, ly, u, z, p;
    
    // Anderson-Darling statistics
    A = -(double)(Nsamp);
    for(n=0; n<Nsamp; n++)
    {
        A -= ( log(gsl_cdf_ugaussian_P(H[n])) + log(1.0 - gsl_cdf_ugaussian_P(H[Nsamp-1-n])) )*(2.*(double)(n+1)-1.0)/(double)(Nsamp);
    }
    
    p = 1.0-AD(Nsamp, A);
    
    
    printf("Hanford %lg %lg\n", A, p);
    
    
    // Anderson-Darling statistics
    A = -(double)(Nsamp);
    for(n=0; n<Nsamp; n++)
    {
        A -= ( log(gsl_cdf_ugaussian_P(L[n])) + log(1.0 - gsl_cdf_ugaussian_P(L[Nsamp-1-n])) )*(2.*(double)(n+1)-1.0)/(double)(Nsamp);
    }
    
    
    p = 1.0-AD(Nsamp, A);
    
    
    printf("Livingston %lg %lg\n", A, p);
    

    
    
    
}



/*
 
 Anderson-Darling test code from
 https://github.com/cran/DescTools/blob/master/src/AnDarl.c
 
 Anderson-Darling test for uniformity.   Given an ordered set
 x_1<x_2<...<x_n
 of purported uniform [0,1) variates,  compute
 a = -n-(1/n)*[ln(x_1*z_1)+3*ln(x_2*z_2+...+(2*n-1)*ln(x_n*z_n)]
 where z_1=1-x_n, z_2=1-x_(n-1)...z_n=1-x_1, then find
 v=adinf(a) and return  p=v+errfix(v), which should be uniform in [0,1),
 that is, the p-value associated with the observed x_1<x_2<...<x_n.
 */


/* Short, practical version of full ADinf(z), z>0.   */
double adinf(double z)
{ if(z<2.) return exp(-1.2337141/z)/sqrt(z)*(2.00012+(.247105-  \
                                                      (.0649821-(.0347962-(.011672-.00168691*z)*z)*z)*z)*z);
    /* max |error| < .000002 for z<2, (p=.90816...) */
    return
    exp(-exp(1.0776-(2.30695-(.43424-(.082433-(.008056 -.0003146*z)*z)*z)*z)*z));
    /* max |error|<.0000008 for 4<z<infinity */
}

/*
 The procedure  errfix(n,x)  corrects the error caused
 by using the asymptotic approximation, x=adinf(z).
 Thus x+errfix(n,x) is uniform in [0,1) for practical purposes;
 accuracy may be off at the 5th, rarely at the 4th, digit.
 */
double errfix(int n, double x)
{
    double c,t;
    if(x>.8) return
        (-130.2137+(745.2337-(1705.091-(1950.646-(1116.360-255.7844*x)*x)*x)*x)*x)/n;
    c=.01265+.1757/n;
    if(x<c){ t=x/c;
        t=sqrt(t)*(1.-t)*(49*t-102);
        return t*(.0037/(n*n)+.00078/n+.00006)/n;
    }
    t=(x-c)/(.8-c);
    t=-.00022633+(6.54034-(14.6538-(14.458-(8.259-1.91864*t)*t)*t)*t)*t;
    return t*(.04213+.01365/n)/n;
}

/* The function AD(n,z) returns Prob(A_n<z) where
 A_n = -n-(1/n)*[ln(x_1*z_1)+3*ln(x_2*z_2+...+(2*n-1)*ln(x_n*z_n)]
 z_1=1-x_n, z_2=1-x_(n-1)...z_n=1-x_1, and
 x_1<x_2<...<x_n is an ordered set of iid uniform [0,1) variates.
 */

double AD(int n,double z){
    double c,v,x;
    x=adinf(z);
    /* now x=adinf(z). Next, get v=errfix(n,x) and return x+v; */
    if(x>.8)
    {v=(-130.2137+(745.2337-(1705.091-(1950.646-(1116.360-255.7844*x)*x)*x)*x)*x)/n;
        return x+v;
    }
    c=.01265+.1757/n;
    if(x<c){ v=x/c;
        v=sqrt(v)*(1.-v)*(49*v-102);
        return x+v*(.0037/(n*n)+.00078/n+.00006)/n;
    }
    v=(x-c)/(.8-c);
    v=-.00022633+(6.54034-(14.6538-(14.458-(8.259-1.91864*v)*v)*v)*v)*v;
    return x+v*(.04213+.01365/n)/n;
}

/* You must give the ADtest(int n, double *x) routine a sorted array
 x[0]<=x[1]<=..<=x[n-1]
 that you are testing for uniformity.
 It will return the p-value associated
 with the Anderson-Darling test, using
 the above adinf() and errfix( ,   )
 Not well-suited for n<7,
 (accuracy could drop to 3 digits).
 */

double ADtest(int n, double *x)
{ int i;
    double t,z=0;
    for(i=0;i<n;i++)   {
        t=x[i]*(1.-x[n-1-i]);
        z=z-(i+i+1)*log(t);}
    return AD(n,-n+z/n);
}




