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

#ifndef IO_EXR_H_
#define IO_EXR_H_


float  *ReadImageEXR(const char fileName[], int *nx, int *ny);

void WriteImageEXR(const char *name, float *pixels,
                          int xRes, int yRes);

void WriteImageEXR(const char *name, float **pixels,
                   int xRes, int yRes);

void writeMultiImageEXR (const char *fileName,
                         float *zPixels,
                         int width,
                         int height,
                         int nchan);

void writeMultiImageEXR (const char *fileName,
                         float **zPixels,
                         int width,
                         int height,
                         int nchan);

float * readMultiImageEXR(const char fileName[],
                          int *width, int *height, int *ncha);

#endif
