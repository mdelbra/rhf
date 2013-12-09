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
#include <limits.h>
#include <math.h>

#include "io_exr.h"
#include "io_png.h"




int main(int argc, char *argv[])
{
    
    int nrow, ncol;
    float *data;
    int i;
    
    if(argc < 2)
    {
        printf("not enough arguments...\n");
        printf("usage:\n");
        printf("    %s imgexr imgpng\n",argv[0]);
        exit(0);
        
        
    }
        

    data = ReadImageEXR(argv[1], &ncol, &nrow);
    
    /*Rescale image to 0-255*/
    for(i=0;i<nrow*ncol*3;i++)
        data[i] = data[i]*255;
    
    

    io_png_write_f32(argv[2], data, (size_t) ncol, (size_t) nrow, 3);
    
    
    free(data);
    
    
    return 0;
    
}
