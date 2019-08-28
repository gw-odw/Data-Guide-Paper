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


int main(int argc, char* argv[]) {
    
    int filtlen = 512;
    
    // Meyer wavelet, dyadic decomp, 512 filter length
    Meyer<double> md(filtlen,0);
    // Meyer wavelet, binary decomp, 512 filter length
    Meyer<double> mb(filtlen,2);
    
    if (argc < 4) {
        
        std::cerr << "Usage: " << argv[0] << " filename " << "f1 " << "f2 " << std::endl;
        return 1;
    }
    
    double fstart = atof(argv[2]); // start of frequency range;
    double fend =   atof(argv[3]); // end of frequency range;
    
    ifstream inFile;  // object for reading from a file
    ofstream outFile; // object for writing to a file
    ofstream plotFile; // object for writing to a plotting file
    int size;
    int cnt;
    double t;
    double f;
    
    int number_of_lines = 0;
    std::string line;
    std::ifstream myfile(argv[1]);
    
    while (std::getline(myfile, line))
        ++number_of_lines;
    std::cout << "Number of lines in input file: " << number_of_lines << endl;
    
    size = number_of_lines;

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

    for(int i = 0; i<x.size(); i++){
	inFile >> time.data[i] >> x.data[i];
    }
    inFile.close();
    
    
    
    double dt;
    double Tobs, fmax;
    
    dt = time.data[1]-time.data[0];
    fmax = 1.0/(2.0*dt);
    Tobs = (double)(x.size())*dt;
    
   cout << dt << " " << fmax << " " << Tobs << endl;
    
    double DT = 1.0;  // length of time chunks in seconds
 
    
    int Tchunks = (int)(Tobs/DT);
 
   
    // Binary decomposition 
    WSeries<double> bin;
    //cout << "Binary decomposition, " << ndecomp << " steps." << endl;
    bin.Forward(x, mb, 1);
    
    int ndecompmax = bin.getMaxLevel();
    int ndecompmin = ndecompmax-7;
    
    if(ndecompmin < 2) ndecompmin = 2;
    
    std::cout << "Maximum level " << ndecompmax << " Minimum level " << ndecompmin << endl;
    
    for(int res = ndecompmin; res<=ndecompmax; res++){
        
     bin.Forward(x, mb, res);
        
    char buffer[32];
    snprintf(buffer, sizeof(char) * 32, "samples_%i.dat", res);
    outFile.open(buffer, ios::out);
        
    char buff[32];
    snprintf(buff, sizeof(char) * 32, "plot_%i.dat", res);
    plotFile.open(buff, ios::out);
        
    int imin = (int)((fstart/fmax)*(double)bin.maxLayer());
    int imax = (int)((fend/fmax)*(double)bin.maxLayer())+1;
        
    f = fmax*(double)imax/(double)bin.maxLayer();
    if(f < fend) imax++;
        
    f = fmax*(double)imin/(double)bin.maxLayer();
    if(f < fstart) imin++;
        
        
    // layer extraction 
    for(int i = imin; i<=imax; i++){
        
        f = fmax*(double)i/(double)bin.maxLayer();
        
        
        wavearray<double> layer;
        bin.getLayer(layer, i);

        for(int j = 0; j<layer.size(); j++){

            outFile << -Tobs/2.0+Tobs*((double)j)/((double)layer.size()-1.0) << " " << f << " " << layer.data[j] << endl;
            plotFile << -Tobs/2.0+Tobs*((double)j)/((double)layer.size()-1.0) << " " << f << " " << layer.data[j] << endl;
        }
      
      // needed for gnuplot
     plotFile << endl;
        
    }
   outFile.close();
   plotFile.close();

   int j;
   j = bin.maxLayer();
   wavearray<double> layer;
   bin.getLayer(layer, j);

   cout << "Frequency layers: " << (j+1) << " Layer size: " << layer.size() << endl;
        
    }
    

    return 0;
}

