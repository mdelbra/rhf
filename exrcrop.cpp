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

/*VERSION 02.08.13*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#include "io_exr.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


int main(int argc, char *argv[])
{
    
    int nx, ny, ncha, cx, cy;
    int xmin, xmax, ymin, ymax;
    int x0,x1,y0,y1;
    int cxy,nxy;
    float *data, *d;
    
    char flag_mc = 1;
    
    if(argc < 6)
    {
        printf("not enough arguments...\n");
        printf("usage:\n");
        printf("    %s x0 y0 x1 y1 img-exr img-crop-exr\n",argv[0]);
        exit(0);
        
    }
    
    data = readMultiImageEXR(argv[5], &nx, &ny, &ncha);
    
    if (!data) {
        data = ReadImageEXR(argv[5], &nx, &ny);
        ncha = 3;
        flag_mc = 0;
        
        if(!data) {
            printf("error :: %s not a correct exr image \n", argv[5]);
            exit(-1);
        }
        
    }
        
        
    x0 = atoi(argv[1]);
    y0 = atoi(argv[2]);
    x1 = atoi(argv[3]);
    y1 = atoi(argv[4]);
    
    xmin = MAX(0,MIN(x0,x1));
    xmax = MIN(nx-1,MAX(x0,x1));
    ymin = MAX(0,MIN(y0,y1));
    ymax = MIN(ny-1,MAX(y0,y1));
    
    cx = xmax - xmin;
    cy = ymax - ymin;
    
    printf("Cropping: xmin:%d xmax:%d ymin:%d ymax:%d\n",xmin,xmax,ymin,ymax);
    //if cx or cy is negative parameters wrong
    if(cx <= 0 || cy <= 0)
    {
        printf("error :: incorrect crop values\n");
        exit(-1);
    }
     
    nxy = nx * ny;
    cxy = cx * cy;   

    d = (float*)malloc(sizeof(float) * cxy * ncha);

    
    /*Crop Image*/
    for(int y=0;y<cy;y++)
        for(int x=0;x<cx;x++)
            for(int c=0;c<ncha;c++)
                d[x + cx*y + cxy*c] = data[x + xmin + nx*(y+ymin) + nxy*c];
    

    
    
    if(flag_mc)
    {
        writeMultiImageEXR(argv[6], d, cx, cy, ncha);
    }
    else
    {
        WriteImageEXR(argv[6], d, cx, cy);
    }
        
    
    free(data);
    free(d);
    
    return 0;
    
}
