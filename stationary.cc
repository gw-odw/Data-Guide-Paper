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

// The WAT library is used https://gitlab.com/gwburst/public/library/tree/public/wat
// Wavelet Analysis Tool
//--------------------------------------------------------------------
// Implementation of
// Bi-othogonal wavelet transforms using lifting scheme
// References:
//   A.Cohen, I.Daubechies, J.Feauveau Bases of compactly supported wavelets
//   Comm. Pure. Appl. Math. 45, 485-560, 1992
//   W. Sweldens - Building your own wavelets at home
//--------------------------------------------------------------------



#include <iostream>
#include <cmath>
#include <climits>
#include <fstream>

// WAT headers
#include "Meyer.hh"
#include "wseries.hh"
#include "wavearray.hh"

using namespace std;

double gaussian(double mean, double var);


int main(int argc, char* argv[]) {
    // Meyer wavelet, dyadic decomp, 1024 filter length
    Meyer<double> md(512,0);
    // Meyer wavelet, binary decomp, 1024 filter length
    Meyer<double> mb(512,2);


    
    if (argc < 5) {
        
        std::cerr << "Usage: " << argv[0] << " filename" << " level" << "f1" << "f2" << std::endl;
        return 1;
    }
    
    ifstream inFile;  // object for reading from a file
    ofstream outFile; // object for writing to a file
    int size;
    int cnt, i;
    double t;
    double f;
    double tr;
    
    int number_of_lines = 0;
    std::string line;
    std::ifstream myfile(argv[1]);
    
    while (std::getline(myfile, line))
        ++number_of_lines;
    std::cout << "Number of lines in input file: " << number_of_lines << endl;
    
    size = number_of_lines;

	// Number of decomposition steps
    int ndecomp = atoi(argv[2]);
    
    double f1, f2;
    
    f1 = atof(argv[3]);
    f2 = atof(argv[4]);

   inFile.open(argv[1], ios::in);
   if (!inFile) {
    cerr << "Can't open input file " << argv[1] << endl;
   exit(1);
   }

   // Expects input data to be time then value
    wavearray<double> x;
    x.resize(size);
    wavearray<double> time;
    time.resize(size);

    for(i = 0; i<x.size(); i++){
	inFile >> time.data[i] >> x.data[i];
    }
    inFile.close();
    
    // Intializes random number generator
    srand(955);
    // testing with unit Gaussian data
   // for(i = 0; i<x.size(); i++) x.data[i] = gaussian(0.0, 1.0);
    
    
    
    double dt;
    double Tobs, fmax;
    
    dt = time.data[1]-time.data[0];
    fmax = 1.0/(2.0*dt);
    Tobs = (double)(x.size())*dt;
    
   cout << dt << " " << fmax << " " << Tobs << endl;
    
    cout << "analysis range " << f1 << " " << f2 << "  Hz" << endl;

   outFile.open("Binary.dat", ios::out);
    // Binary decomposition 
    WSeries<double> bin;
    //cout << "Binary decomposition, " << ndecomp << " steps." << endl;
    bin.Forward(x, mb, ndecomp);
    // layer extraction
    
    int j;
    j = bin.maxLayer();
    wavearray<double> layer;
    bin.getLayer(layer, j);
    
    cout << "Frequency layers: " << (j+1) << " " << bin.maxLayer() << " Layer size: " << layer.size() << endl;
    
    double DT, DF;
    
    DT = Tobs/(double)(layer.size());
    DF = fmax/(double)(j+1);
    
    cout << "Frequency spacing (Hz): " << DF << " Time Spacing (s) " << DT << endl;
    
    int tsamp = layer.size();
    double *powt, *pown, *powa;
    
    int k, kk;
    
    kk = 0;
    
    for(i = 0; i<=bin.maxLayer(); i++){
        f = ((double)i)*DF;
        if(f >= f1 && f <= f2) kk++;
    }
    
    cout << "Frequency samples in analysis range " << kk << endl;

    
    powa = (double*) malloc (tsamp*sizeof(double));
    powt = (double*) malloc (tsamp*sizeof(double));
    pown = (double*) malloc (tsamp*sizeof(double));
    
    for(int j = 0; j<tsamp; j++) powt[j] = 0.0;
    
    
    for(i = 0; i<=bin.maxLayer(); i++){
        wavearray<double> layer;
        bin.getLayer(layer, i);
        //cout << "Frequency layer: " << i << endl;
        for(int j = 0; j<layer.size(); j++){
	  //cout << layer.data[j] << " ";
            f = ((double)i)*DF;
            outFile << -Tobs/2.0+(0.5+(double)j)*DT << " " << f+0.5*DF << " " << layer.data[j] << endl;
            if(f >= f1 && f <= f2) powt[j] += layer.data[j]*layer.data[j];
        }
       outFile  << endl ;
    }
   outFile.close();
    
    double sigma;
    sigma = sqrt(2.0/(double)(k));
    
    double mk;
    mk = sqrt((double)(2*k)-1.0);
    
    //for(int j = 0; j<tsamp; j++) pown[j] = sqrt(2.0*powt[j])-mk;
    
    for(int j = 0; j<tsamp; j++) powa[j] = powt[j]/(double)(kk);
    
    
    // Wilson-Hilferty
    for(int j = 0; j<tsamp; j++) pown[j] = (pow(powt[j]/(double)(kk),1.0/3.0)-(1.0-2.0/((double)(9*kk))))/sqrt(2.0/((double)(9*kk)));
    
    // Standard
    //for(int j = 0; j<tsamp; j++) pown[j] = (powt[j]-(double)(kk))/sqrt((double)(2*kk));
    
    
    
    double chisq;
    chisq = 0.0;
    
    double mn, var;
    
    tr = Tobs/2.0-2.0;
    if(4.0*DT > 2.0)
    {
    tr = Tobs/2.0-4.0*DT;  // remove edges to quash filter effects
    }
    
    cout << "Time range " << -tr <<" " << tr << endl;
    
    
    mn = 0.0;
    var = 0.0;
    k = 0;
    for(int j = 0; j<tsamp; j++)
    {
        t = -Tobs/2.0+(0.5+(double)j)*DT;
        if(t >= -tr && t <= tr)
        {
            mn += pown[j];
            var += pown[j]*pown[j];
            k++;
            chisq += pown[j]*pown[j];
        }
    }
    
    mn /= (double)k;
    var /= (double)k;
    
    cout << "Mean " << mn << " Deviation " << sqrt(var-mn*mn) << " Theoretical Deviation " << 1.0 << endl;
    


    chisq = 0.0;
    mn = 0.0;
    var = 0.0;
    k = 0;
    for(int j = 0; j<tsamp; j++)
    {
        t = -Tobs/2.0+(0.5+(double)j)*DT;
        if(t >= -tr && t <= tr)
        {
            mn += powt[j];
            var += powt[j]*powt[j];
            k++;
            chisq += (powt[j])*(powt[j]);
        }
    }
    
    mn /= (double)k;
    var /= (double)k;
    
    cout << "Mean " << mn << " Deviation " << sqrt(var-mn*mn)  << endl;
    
    
    outFile.open("Scales.dat", ios::out);
    outFile <<  kk << " " << k << endl;
    outFile.close();
    
    
    outFile.open("PowerTime.dat", ios::out);
    
    for(int j = 0; j<tsamp; j++){
        t = -Tobs/2.0+(0.5+(double)j)*DT;
        if(t >= -tr && t <= tr)
        {
        outFile << t << " " << pown[j] << " " << powt[j] << " " << powa[j] << endl;
        }
    }
    outFile.close();
    
    
 
    return 0;
}

double gaussian(double mean, double var){
    
    double x1=0.0, x2=0.0, y1=0.0, y2=0.0, w=2.0;
    
    while( w > 1.0 ){
        x1 = 2*(rand() / float(INT_MAX))-1;
        x2 = 2*(rand() / float(INT_MAX))-1;
        w = x1*x1 + x2*x2;
    }
    
    w = sqrt(-2.0*log(w)/w);
    return (x1*w)*sqrt(var)+mean;
}



