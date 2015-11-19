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

#include "libauxiliar.h"
#include "math.h"
#include <climits>



void fpClear(float *fpI,float fValue, int iLength) {
    for (int ii=0; ii < iLength; ii++) fpI[ii] = fValue;
}


float fiChiSquareNDfFloatDist(int *df, float *k0, float *k1, float *u0,
                              float *u1,int i0,int j0,int i1,int j1,int radius,
                              int width0, int width1) {
    
    float dist=0.0;
    for (int s=-radius; s<= radius; s++) {
        
        int l = (j0+s)*width0 + (i0-radius);
        float *ptr0 = &u0[l];
        float *ptrK0 = &k0[l];
        
        l = (j1+s)*width1 + (i1-radius);
        float *ptr1 = &u1[l];
        float *ptrK1 = &k1[l];
        
        
        for (int r=-radius; r<=radius; r++,ptr0++, ++ptrK0, ptr1++, ++ptrK1) {
            
            float sum = (*ptr0 + *ptr1);
            
            if(sum>1.00f)//to avoid problems due to little values.
            {
                float dif =  (*ptr0)*(*ptrK1) - (*ptr1)*(*ptrK0);
                dist += (dif*dif)/((*ptrK0)*(*ptrK1)*sum);
                (*df)++;
            }
            
        }
        
    }
    
    return dist;
}



float fiChiSquareNDfFloatDist(int *df, float **u0,float **u1,int i0,int j0,
                              int i1,int j1,int radius,int channels, 
                              int width0, int width1) {
    
    float dif = 0.0f;
    
    for (int ii=0; ii < channels-1; ii++) {
        
        dif += fiChiSquareNDfFloatDist(df, u0[channels-1], u1[channels-1],
                                       u0[ii], u1[ii], i0, j0,i1,j1,radius,
                                       width0, width1);
        
    }
    
    
    return dif;
}


 


/*    
The functions gaussian_kernel, gaussian_kernel_notNorm and gaussian_sampler
are based on those from the LSD algorithm.
*/

/*----------------------------------------------------------------------------
 
 LSD - Line Segment Detector on digital images
 
 Copyright 2007-2010 rafael grompone von gioi (grompone@gmail.com)
 
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

    
/*----------------------------------------------------------------------------*/
/*----------------------------- Gaussian filter ------------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** Compute a Gaussian kernel of length 'kernel->dim',
 standard deviation 'sigma', and centered at value 'mean'.
 
 For example, if mean=0.5, the Gaussian will be centered
 in the middle point between values 'kernel->values[0]'
 and 'kernel->values[1]'.
 */
static void gaussian_kernel(float * kernel, int ksize, float sigma, float mean)
{
    float sum = 0.0;
    float val;
    int i;
    
    /* check parameters */
    
    /* compute Gaussian kernel */
    for(i=0;i<ksize;i++)
    {
        val = ( (float) i - mean ) / sigma;
        kernel[i] = (float)exp( -0.5 * val * val );
        sum += kernel[i];
    }
    
    /* normalization */
    if( sum >= 0.0 ) for(i=0;i<ksize;i++) kernel[i] /= sum;
}




/*----------------------------------------------------------------------------*/
/** Scale the input image 'in' by a factor 'scale' by Gaussian sub-sampling.
 
 For example, scale=0.8 will give a result at 80% of the original size.
 
 The image is convolved with a Gaussian kernel
 @f[
 G(x,y) = \frac{1}{2\pi\sigma^2} e^{-\frac{x^2+y^2}{2\sigma^2}}
 @f]
 before the sub-sampling to prevent aliasing.
 
 The standard deviation sigma given by:
 -  sigma = sigma_scale / scale,   if scale <  1.0
 -  sigma = sigma_scale,           if scale >= 1.0
 
 To be able to sub-sample at non-integer steps, some interpolation
 is needed. In this implementation, the interpolation is done by
 the Gaussian kernel, so both operations (filtering and sampling)
 are done at the same time. The Gaussian kernel is computed
 centered on the coordinates of the required sample. In this way,
 when applied, it gives directly the result of convolving the image
 with the kernel and interpolated to that particular position.
 
 A fast algorithm is done using the separability of the Gaussian
 kernel. Applying the 2D Gaussian kernel is equivalent to applying
 first a horizontal 1D Gaussian kernel and then a vertical 1D
 Gaussian kernel (or the other way round). The reason is that
 @f[
 G(x,y) = G(x) * G(y)
 @f]
 where
 @f[
 G(x) = \frac{1}{\sqrt{2\pi}\sigma} e^{-\frac{x^2}{2\sigma^2}}.
 @f]
 The algorithm first applies a combined Gaussian kernel and sampling
 in the x axis, and then the combined Gaussian kernel and sampling
 in the y axis.
 */
float** gaussian_sampler( float** in, int nx, int ny, int nch,
                         int* nx_o, int* ny_o,
                        float scale, float sigma_scale )
{
    
    float *kernel, **aux, **out;
    
    
    int N,M,h,n,x,y,i,ii;
    int xc,yc,j,double_x_size,double_y_size;
    float sigma,xx,yy,sum,prec;
    
    /* compute new image size and get memory for images */
    N = (int) ceil( nx * scale );
    M = (int) ceil( ny * scale );
    
    aux = (float**)malloc(sizeof(float*)*nch);
    out = (float**)malloc(sizeof(float*)*nch);;
    
    for (int ii=0; ii < nch; ii++) {
        aux[ii] = (float*)malloc(sizeof(float)*N*ny);
        out[ii] = (float*)malloc(sizeof(float)*N*M);
    }
    
    
    /* sigma, kernel size and memory for the kernel */
    //sigma = scale < 1.0 ? sigma_scale / scale : sigma_scale;
    
    sigma = scale <1.0 ? sigma_scale * sqrt(1/(scale*scale) - 1) : sigma_scale;
    
    /*The size of the kernel is selected to guarantee that the
     the first discarded term is at least 10^prec times smaller
     than the central value. For that, h should be larger than x, with
     e^(-x^2/2sigma^2) = 1/10^prec.
     Then,
     x = sigma * sqrt( 2 * prec * ln(10) ).
     */
    prec = 2.0;
    h = (int) ceil( sigma * sqrt( 2.0 * prec * log(10.0) ) );
    n = 1+2*h; /* kernel size */
    
    kernel = (float*) malloc(n*sizeof(float));
    
    /* auxiliary double image size variables */
    double_x_size = (int) (2 * nx);
    double_y_size = (int) (2 * ny);
    
    /* First subsampling: x axis */
    for(x=0;x<N;x++)
    {
        /*
         x   is the coordinate in the new image.
         xx  is the corresponding x-value in the original size image.
         xc  is the integer value, the pixel coordinate of xx.
         */
        xx = ((float) x +0.5f) / scale;
        /* coordinate (0.0,0.0) is in the center of pixel (0,0),
         so the pixel with xc=0 get the values of xx from -0.5 to 0.5 */
        xc = (int) floor( xx );
        
        /*If the image is an histogram image do not normalize the 
         interpolation kernel (calculate the histogram union)*/

       // if(nch>3)
       //     gaussian_kernel_notNorm( kernel, n, sigma, 
       //                             (float) h + xx - (float) xc -0.5f);
       // else
            gaussian_kernel(kernel, n, sigma, (float)h +xx -(float) xc -0.5f);
        
        /* the kernel must be computed for each x because the fine
         offset xx-xc is different in each case */
        
        for(y=0;y<ny;y++)
        {
            for(ii=0;ii<nch;ii++)
            {
                sum = 0.0;
                for(i=0;i<n;i++)
                {
                    
                    j = xc - h + i ;

                    /* symmetry boundary condition */
                    while( j < 0 ) j += double_x_size;
                    while( j >= double_x_size ) j -= double_x_size;
                    if( j >= (int) nx ) j = double_x_size-1-j;
                
                    sum += in[ii][ j + y * nx] * kernel[i];
                }
                aux[ii][x + y * N] = sum;
            }
        }
    }
    
    /* Second subsampling: y axis */
    for(y=0;y<M;y++)
    {
        /*
         y   is the coordinate in the new image.
         yy  is the corresponding x-value in the original size image.
         yc  is the integer value, the pixel coordinate of xx.
         */
        yy = ((float) y + 0.5f) / scale;

        /* coordinate (0.0,0.0) is in the center of pixel (0,0),
         so the pixel with yc=0 get the values of yy from -0.5 to 0.5 */
        //yc = (int) floor( yy + 0.5 );
        yc = (int) floor( yy  );

        
        gaussian_kernel( kernel, n, sigma, (float)h +yy - (float) yc -0.5f);

        /* the kernel must be computed for each y because the fine
         offset yy-yc is different in each case */
        
        for(x=0;x<N;x++)
        {
            for(ii=0;ii<nch;ii++)
            {
                sum = 0.0;
                for(i=0;i<n;i++)
                {
                    j = yc - h + i;
                
                    /* symmetry boundary condition */
                    while( j < 0 ) j += double_y_size;
                    while( j >= double_y_size ) j -= double_y_size;
                    if( j >= (int) ny ) j = double_y_size-1-j;
                    
                    sum += aux[ii][ x + j * N] * kernel[i];
                }
                out[ii][ x + y * N ] = sum;
            }
        }
    }
    
    *nx_o = N;
    *ny_o = M;
    
    /* free memory */
    for (int ii=0; ii < nch; ii++) {
        free(aux[ii]);
    }
 
    free(kernel);
    free(aux);
  
    
    return out;
}


 


/* extract image value (even outside image domain) */
float v(float *in, int x, int y, int nx, int ny, float bg)
{
    if (x<0 || x>=nx || y<0 || y>=ny)
        return bg; else return in[y*nx+x];
}

/* extract image value (even outside image domain) */
float v(float *in, int x, int y, int nx, int ny)
{
 
    x = MIN(x,2*nx-x-1);
    y = MIN(y,2*ny-y-1);
    
    x = MAX(x,-x-1);
    y = MAX(y,-y-1);
    
    return in[y*nx+x];
}

/* c[] = values of interpolation function at ...,t-2,t-1,t,t+1,... */

/* coefficients for cubic interpolant (Keys' function) */
void keys(float *c,float t,float a)
{
    float t2,at;
    
    t2 = t*t;
    at = a*t;
    c[0] = a*t2*(1.f-t);
    c[1] = (2.f*a+3.f - (a+2.f)*t)*t2 - at;
    c[2] = ((a+2.f)*t - a-3.f)*t2 + 1.f;
    c[3] = a*(t-2.f)*t2 + at;
}

/*------------------------ MAIN MODULE ---------------------------------*/

float** bicubic_interpolation(float** in, int nx, int ny, int nch,
                              int nsx, int nsy)
{
    
    int    x,y,xi,yi,d,k;
    float  zx,zy,res,xp,yp,u,c[12];
    float **ref, **tmp, **out;
    
    
    /* COMPUTE OUTPUT DIMENSIONS AND ZOOM FACTORS */
    zx = ((float)nsx) / ((float)nx);
    zy = ((float)nsy) / ((float)ny);
   
    
    tmp = (float**) malloc(sizeof(float*) * nch);
    out = (float**) malloc(sizeof(float*) * nch);

    for(k=0;k<nch;k++)
    {
        tmp[k] = (float*)malloc(sizeof(float)*ny*nsx);
        out[k] = (float*)malloc(sizeof(float)*nsy*nsx);
    }
    
    ref = in;
    /********** FIRST LOOP (x) **********/
    
    for (x=0;x<nsx;x++) {
        
        xp = ( (float)x + 0.5f )/zx;
        
        if (xp<0. || xp>(float)nx) 
            for (y=0;y<ny;y++) 
                for(k=0;k<nch;k++) 
                    tmp[k][y*nsx+x]=0.0f; 
                
                else {
                    
                    xp -= 0.5;
                    xi = (int)floor((float)xp); 
                    u = xp-(float)xi;
                    
                    keys(c,u,-0.5);
                   
                    /* this test saves computation time */
                    if (xi-1>=0 && xi+2<nx) {
                        for (y=0;y<ny;y++) {
                            for (k=0;k<nch;k++) {
                                for (d=-1,res=0.;d<=2;d++) 
                                    res += c[2-d]*ref[k][y*nx+xi+d];
                                tmp[k][y*nsx+x] = res;
                            }
                        }
                        
                    } else 
                        for (y=0;y<ny;y++) {
                            for(k=0;k<nch;k++) {
                                for (d=-1,res=0.;d<=2;d++) 
                                    res += c[2-d]*v(ref[k],xi+d,y,nx,ny);
                                tmp[k][y*nsx+x] = res;
                            }
                        }
                }
        
    }
    
    ref = tmp;
    
    /********** SECOND LOOP (y) **********/
        
        for (y=0;y<nsy;y++) {
            
            yp = ( (float)y + 0.5f )/zy;
            
            if (yp<0. || yp>(float)ny) 
                for (x=0;x<nsx;x++) 
                     for(k=0;k<nch;k++) 
                         out[k][y*nsx+x]=0.0f; 
                    
                    else {
                        
                        yp -= 0.5;
                        
                        yi = (int)floor((float)yp); 
                        u = yp-(float)yi;
                        keys(c,u,-0.5); 
                        
                        /* this test saves computation time */
                        if (yi-1>=0 && yi+2<ny) {
                            for (x=0;x<nsx;x++) {
                                for(k=0;k<nch;k++) {
                                    for (d=-1,res=0.;d<=2;d++) 
                                        res += c[2-d]*ref[k][(yi+d)*nsx+x];
                                    out[k][y*nsx+x] = res;
                                }
                            }
                        } else 
                            for (x=0;x<nsx;x++) {
                                for(k=0;k<nch;k++) {
                                    for (d=-1,res=0.;d<=2;d++) 
                                        res += c[2-d]*v(ref[k],x,yi+d,nsx,ny);
                                    out[k][y*nsx+x] = res;
                                }
                            }
                    }
        }
    
    
    /*Clean*/
    for(k=0;k<nch;k++)
        free(tmp[k]);
    
    free(tmp);
    
    return out;
}

int compute_filter_mask(float** in, int x, int y, int nx, int ny, int nch)
{
    
    //Assume RGB Channels!
    int nc = (nch-1)/3;
    int mask = 0;
    
    int i=0;
    float a1=0;
    float a2=0;
    float a3=0;
    float b1=0;
    float b2=0;
    float b3=0;
    
    float nr = in[nch-1][x + y*nx];
    while (i<nc && (mask<3))
    {
        a1=in[i][x + y*nx];
        
        if(abs(a1+b1-nr) < 1e-3)
            mask++;
        
        b1=a1;
        
        a2=in[i+nc][x + y*nx];
        if(abs(a2+b2-nr) < 1e-3)
            mask++;
        
        b2=a2;
        
        a3=in[i+2*nc][x + y*nx];
        if(abs(a3+b3-nr) < 1e-3)
            mask++;
        
        b3=a3;
        i++;
        
    }
    
    return mask;
    
}

void compute_knn_index(int k, float *ivect_dist, int *ovect_ind,  int n)
{
    
    int *ind = new int[n];
    
    for (int i=0;i<n;i++)
        ind[i] = i;
    
    
    float minv;
    int minind;
    
    /*Outer Loop*/
    for(int i=0; i<k; i++)
    {
        
        minv = ivect_dist[i]; /*Big Number*/
        minind = i;
        
        /*inner loop: find minimum value*/
        for(int j=i+1; j<n; j++)
        {
            if (ivect_dist[j] < minv) 
            {
                minv = ivect_dist[j];
                minind = j;
            }
            
        }
        
        /*Swap index*/
        int ind_aux = ind[i];
        ind[i] = ind[minind];
        ind[minind] = ind_aux;
        
        /*Swap values*/
        float val_aux = ivect_dist[i];
        ivect_dist[i] = ivect_dist[minind];
        ivect_dist[minind] = val_aux;
        
    }
    //Return all the indices also the part no sorted
    for(int i=0;i < n; i++)
        ovect_ind[i] = ind[i];
    
    delete[] ind;
    
}

