/*----------------------------------------------------------------------------
 
 RHF - Ray Histogram Fusion
 
 Copyright (c) 2013, A. Buades <toni.buades@uib.es>,
 M. Delbracio <mdelbra@gmail.com>, 
 J-M. Morel <morel@cmla.ens-cachan.fr>,
 P. Muse <muse@fing.edu.uy>
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as
 published by the Free Software Foundation, either version 3 of the
 License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU Affero General Public License for more details.
 
 You should have received a copy of the GNU Affero General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.
 
 ----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "io_exr.h"

 
int main(int argc, char **argv) {


    // read input EXR images
    int nx1,ny1,nx2,ny2,nc;
    int nxyc;
    
    float *img1, *img2;
    
    float psnr;

    if(argc < 3)
    {
        printf("not enough arguments...\n");
        printf("usage:\n");
        printf("    %s img1.exr img2.exr [sigma] [imgdiff]\n",argv[0]);
        exit(0);
        
    }
    /*Only works with EXR color images*/
    nc = 3;

    img1 = ReadImageEXR(argv[1], &nx1, &ny1);    
    img2 = ReadImageEXR(argv[2], &nx2, &ny2);
    
    nxyc = nx1*ny1*nc;
    
    double acumdif = 0.0f;

    //Compute difference image and PSNR
    //Compute difference and convert from [-4 sigma, 4 sigma] to [0,1]

    if(argv[3] && argv[4])
    {
        float sigma = atof(argv[3]);
        sigma *= 4.0f;
        
        float *difference = new float[nxyc];
                
        for (int i=0; i < nxyc ;  i++) {
            
            double dif = img1[i] - img2[i];
                
            //If one of the two pixels is not saturated
            if(img1[i]<=1 || img2[i]<=1)
            {
                acumdif += dif*dif; 
            }
            
            dif = (dif + sigma) / (2.0f * sigma);
            difference[i] = dif;
            
            
        }
        
        WriteImageEXR(argv[4], difference, nx1, ny1);
        
    }
    //Do not compute difference image only compute PSNR
    else {
       
        for (int i=0; i < nxyc ;  i++) {
            
            //If one of the two pixels is not saturated
            if(img1[i]<=1 || img2[i]<=1)
            {
                double dif = img1[i] - img2[i];
                acumdif += dif*dif; 
                
            }
        }
        
    }
    
    acumdif = acumdif/nxyc;
    psnr = (float) 10*log10(1/acumdif);
    
    printf("%5.3f\n",psnr);
    
    return 0;

}



