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

#ifndef _LIBAUXILIAR_H_
#define _LIBAUXILIAR_H_


#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cmath>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <cassert>



#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#define MIN(i,j) ( (i)<(j) ? (i):(j) )


#define dTiny 1e-10
#define fTiny 0.00000001f
#define fLarge 100000000.0f


using namespace std;

void fpClear(float *fpI,float fValue, int iLength);


float fiChiSquareNDfFloatDist(int *df, float **u0,float **u1,int i0,int j0,
                              int i1,int j1,int radius,int channels,
                              int width0, int width1);

float fiChiSquareNDfFloatDist(int *df, float *k0,float *k1,float *u0,float *u1,
                              int i0,int j0,int i1,int j1,int radius,
                              int width0, int width1);


void compute_knn_index(int k, float *ivect_dist, int *ovect_ind,  int n);


float** gaussian_sampler( float** in, int nx, int ny, int nch,
                         int* nx_o, int* ny_o,
                         float scale, float sigma_scale );


float** bicubic_interpolation(float** in, int nx, int ny, int nch,
                              int nxs, int nys);

int compute_filter_mask(float** in, int x, int y, int nx, int ny, int nch);

void alpha_mul (float *rgb, float *alpha, int np);

void alpha_div (float *rgb, float *alpha, int np);

#endif
