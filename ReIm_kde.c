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

// gcc -o ReIm_kde ReIm_kde.c -lm -lgsl


int main(int argc, char *argv[])
{

    char filename[1024];
    
    FILE *in;
    FILE *out;
    
    int N, M, i, n, ix, iy;
    double f, t ;
    double x, y, z, u;
    double dx, dy;
    double Re, Im;
    double sigma2;
    double range;
    
    int NG, NT;
    
    NG = 300;  // x, y grid size
    NT = NG*NG;
    
    range = 3.5;   // x,y go from -range to range
    
    dx = 2.0*range/(double)(NG);
    dy = 2.0*range/(double)(NG);
    
    double *hist = malloc (NT*sizeof(double));
    
    for(n=0; n<NT; n++) hist[n] = 0.0;
    
    // kde smoothing scale
    sigma2 = 0.1*0.1;
    u = 1.0/(TPI*sigma2);
    
    
    N = 119809;
    
   // N = 10000;
    

    in=fopen("freq_RI.dat","r");
    for(n=0; n<N; n++)
    {
        fscanf(in,"%lf %lf %lf", &f, &Re, &Im);
        
        for(ix=0; ix < NG; ix++)
        {
            x = -range + (double)(ix)*dx;
            for(iy=0; iy < NG; iy++)
            {
              y = -range + (double)(iy)*dy;
              z = ((x-Re)*(x-Re)+(y-Im)*(y-Im))/(2.0*sigma2);
              hist[ix*NG+iy] += u*exp(-z);
            }
        }
        
     }
    fclose(in);
    
    for(n=0; n<NT; n++) hist[n] /= (double)(N);
    
    
    out = fopen("freq_kde.dat", "w");
    
    for(ix=0; ix < NG; ix++)
    {
        x = -range + (double)(ix)*dx;
        for(iy=0; iy < NG; iy++)
        {
            y = -range + (double)(iy)*dy;
            fprintf(out, "%f %f %f\n", y, x, hist[ix*NG+iy]);
        }
        fprintf(out, "\n");
    }

    fclose(out);

    
    free(hist);
    
    return 0;
    
}



