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

#include "libdenoising.h"
#include "io_exr.h"

//small value
#define EPSILON 0.0001

// Gaussian subsampling blur
#define SIGMASCALE 0.55

void rhf_multiscale(int iDWin,       // Half size of patch
                    int iDBloc,      // Half size of research window
                    float fDistance, // Max-Distance parameter
                    int knn,         // Minimum Number of Neighbors
                    int iNscales,    // Number of Scales
                    float **fhI,     // Histogram
                    float **fpI,     // Input
                    float **fpO,     // Output
                    int iChannels,   // Number of channels
                    int iWidth,      // Image width
                    int iHeight,     // Image height
                    int iBins)       // Number of bins Histogram image

{
    
    //0 Scales:  1x1
    //1 Scales:  1x1, 2x2
    //2 Scales:  1x1, 2x2, 4x4
    //3 Scales:  1x1, 2x2, 4x4, 16x16
    
    //Need a buffer....
    float **fpOs_old = (float**) malloc(sizeof(float*)*iChannels);
    int nxSold,nySold;
    nxSold=nySold=-1;
    
    for (int ii=0; ii< iChannels;ii++)
        fpOs_old[ii] = (float*) malloc(sizeof(float)*iWidth * iHeight);
    
    double sigma_scale = SIGMASCALE;
    
    //Total number of samples in the whole image
    double dtotal = 0.0f;
    
    for(int ii=0;ii < iWidth*iHeight; ii++)
    {
        dtotal += fhI[iBins-1][ii];
    }
    
    for(int s=iNscales-1; s>=0; s--)
    {
        
        // Generate Image and Histogram at desired scale
        float **fpIs;// = (float**) malloc(sizeof(float*)*iChannels);
        float **fpOs = (float**) malloc(sizeof(float*)*iChannels);
        float **fhIs;// = (float**) malloc(sizeof(float*)*iBins);
        
        int nxS, nyS, nyM, nxM, nxSS, nySS;
        
        double scale = pow(0.5f,s);
        
        if(s>0) //If it is not the last scale...
        {
            fhIs = gaussian_sampler(fhI,
                                    iWidth, iHeight, iBins,
                                    &nxS, &nyS,
                                    (float)scale, (float)sigma_scale);
            
            //Renormalize weights to keep total number of samples
            double dtotalS = 0;
            for(int ii=0;ii < nxS*nyS; ii++)
            {
                dtotalS += fhIs[iBins-1][ii];
            }
            float samples_factor = (float)(dtotal/dtotalS);
            
            for(int ii=0; ii < nxS*nyS; ii++)
            {
                for(int jj=0; jj< iBins; jj++)
                    fhIs[jj][ii] = samples_factor*fhIs[jj][ii];
            }
            
            fpIs = gaussian_sampler(fpI,
                                    iWidth, iHeight, iChannels,
                                    &nxS, &nyS,
                                    (float)scale, (float)sigma_scale);
        }
        else
        {
            fpIs = fpI;
            fhIs = fhI;
            nxS = iWidth;
            nyS = iHeight;
        }
        
        for (int ii=0; ii < iChannels; ii++)
        {
            fpOs[ii] = (float*) malloc(sizeof(float)*nxS*nyS);
        }
        
        //Filter Scale
        printf("-->Filtering Scale %d BEGIN\n",s);
        
        //In the scale 0 force a minimum number of neighbors
        int knnS = (s>0)?0:knn;
        
        if(knnS>0)
        {
            rhf_knn(iDWin,       // Half size of patch
                    iDBloc,      // Half size of research window
                    fDistance,   // Max-Distance parameter
                    knn,         // Minimum number of neighbors
                    fhIs,        // Histogram
                    fpIs,        // Input
                    fpOs,        // Output
                    iChannels, nxS, nyS, iBins);
        }
        else
        {
            rhf(iDWin,      // Half size of patch
                iDBloc,     // Half size of research window
                fDistance,  // Max-Distance parameter
                fhIs,       // Histogram
                fpIs,       // Input
                fpOs,       // Output
                iChannels, nxS, nyS, iBins);
        }
        
        printf("-->Filtering Scale %d END\n",s);
        
        if(s < iNscales - 1) //This is not the last Scale
        {
            
            float **fpOs_PosTerm; //= (float**) malloc(sizeof(float*)*iChannels);
            
            fpOs_PosTerm = bicubic_interpolation(fpOs_old, nxSold, nySold,
                                                 iChannels, nxS, nyS);
            nxM = nxS;
            nyM = nyS;
            
            float **fpOs_NegTermD; //= (float**) malloc(sizeof(float*)*iChannels);
            float **fpOs_NegTerm; //= (float**) malloc(sizeof(float*)*iChannels);
            
            fpOs_NegTermD = gaussian_sampler(fpOs,
                                             nxS, nyS, iChannels,
                                             &nxM, &nyM,
                                             0.5f,    (float)sigma_scale);
            
            fpOs_NegTerm = bicubic_interpolation(fpOs_NegTermD, nxM, nyM,
                                                 iChannels, nxS, nyS);
            nxSS = nxS;
            nySS = nyS;
            
            nxS = MIN(nxS, nxSS);
            nyS = MIN(nyS, nySS);
            
            for(int ii=0; ii < iChannels; ii++)
                for(int x=0;x<nxS*nyS;x++)
                {
                    fpOs[ii][x] += fpOs_PosTerm[ii][x] - fpOs_NegTerm[ii][x];
                }
            
            //Cleaning
            for (int ii=0; ii < iChannels; ii++) {
                free(fpOs_NegTermD[ii]);
                free(fpOs_NegTerm[ii]);
                free(fpOs_PosTerm[ii]);
            }
            
            free(fpOs_NegTermD);
            free(fpOs_NegTerm);
            free(fpOs_PosTerm);
        }
        
        //Clean last fpOs_old and create the new one
        nxSold = nxS;
        nySold = nyS;
        
        for (int ii=0; ii < iChannels; ii++) {
            for(int x=0;x<nxS*nyS;x++)
                fpOs_old[ii][x] = fpOs[ii][x];
        }
        
        //If scale==0 copy output
        if(s==0)
        {
            for(int ii=0; ii < iChannels; ii++)
                for(int x=0;x<iWidth*iHeight;x++)
                    fpO[ii][x] = fpOs[ii][x];
            
            /*Clean*/
            for (int ii=0; ii < iChannels; ii++) {
                free(fpOs_old[ii]);
            }
            
            free(fpOs_old);
        }
        
        //Clean the Step
        if(s>0)
        {
            for (int ii=0; ii < iBins; ii++) {
                free(fhIs[ii]);
            }
            
            for (int ii=0; ii < iChannels; ii++) {
                free(fpIs[ii]);
            }
            
            free(fpIs);
            free(fhIs);
        }
        
        for (int ii=0; ii < iChannels; ii++) {
            free(fpOs[ii]);
        }
        
        free(fpOs);
        
    } //for Scales
}


void rhf_knn(int iDWin,       // Half size of patch
             int iDBloc,      // Half size of search window
             float fDistance, // Max-Distance parameter
             int knn,         // Minimum k-nearest neighbours
             float **fhI,     // Histogram Image
             float **fpI,     // Input image
             float **fpO,     // Output image
             int iChannels,   // Number of channels
             int iWidth,      // Image width
             int iHeight,     // Image height
             int iBins)       // Number of bins Histogram image
{
    
    printf("---->rhf_knn: dmax = %f, k = %d\n", fDistance,knn);
    
    //k nearest neighbors + the current patch
    int knnT = knn + 1;
    
    // length of each channel
    int iwxh = iWidth * iHeight;
    
    //  length of comparison window
    int ihwl = (2*iDWin+1);
    int iwl = (2*iDWin+1) * (2*iDWin+1);
    
    // auxiliary variable
    // number of denoised values per pixel
    float *fpCount = new float[iwxh];
    fpClear(fpCount, 0.0f,iwxh);
    
    // clear output
    for (int ii=0; ii < iChannels; ii++) fpClear(fpO[ii], 0.0f, iwxh);
    
    // PROCESS STARTS
    // for each pixel (x,y)
#pragma omp parallel shared(fpI, fpO)
    {
        
#pragma omp for schedule(dynamic) nowait
        
        for (int y=0; y < iHeight ; y++) {
            
            // auxiliary variable
            // denoised patch centered at a certain pixel
            float **fpODenoised = new float*[iChannels];
            for (int ii=0; ii < iChannels; ii++)
                fpODenoised[ii] = new float[iwl];
            
            for (int x=0 ; x < iWidth;  x++) {
                
                //Reduce the size of comparison window  near the boundary
                int iDWin0 = MIN(iDWin,MIN(iWidth-1-x,
                                           MIN(iHeight-1-y,MIN(x,y))));
                
                //Research zone depending on the boundary/size of the window
                int imin=MAX(x-iDBloc,iDWin0);
                int jmin=MAX(y-iDBloc,iDWin0);
                
                int imax=MIN(x+iDBloc,iWidth-1-iDWin0);
                int jmax=MIN(y+iDBloc,iHeight-1-iDWin0);
                
                
                //  clear current denoised patch
                for (int ii=0; ii < iChannels; ii++)
                    fpClear(fpODenoised[ii], 0.0f, iwl);
                
                /*Check if we need to denoise this pixel!!*/
                // sum of weights
                float fTotalWeight = 0.0f;
                
                // weights
                float fWeight = 1.0f;
                
                int dj = jmax-jmin+1;
                int di = imax-imin+1;
                
                int *ovect_ind = new int[dj*di];
                float *fDif_all = new float[dj*di];
                
                for (int j=jmin; j <= jmax; j++)
                    for (int i=imin ; i <= imax; i++)
                    {
                        
                        int df = 0;
                        float fDifHist = fiChiSquareNDfFloatDist(&df,fhI,fhI,
                                                                 x,y,
                                                                 i,j,iDWin0,
                                                                 iBins,
                                                                 iWidth,
                                                                 iWidth);
                        
                        fDif_all[(j-jmin)+(i-imin)*dj] = fDifHist/(df+(float)EPSILON);
                        
                    }
                
                compute_knn_index(knnT, fDif_all, ovect_ind, dj*di);
                
                //ALWAYS: select at least KNN similar patchs.
                int kk;
                for(kk=0;kk<knnT;kk++)
                {
                    
                    fTotalWeight += fWeight;
                    
                    //Reconvert index
                    int i = ovect_ind[kk]/dj + imin;
                    int j = ovect_ind[kk]%dj + jmin;
                    
                    for (int is=-iDWin0; is <=iDWin0; is++) {
                        int aiindex = (iDWin+is) * ihwl + iDWin;
                        int ail = (j+is)*iWidth+i;
                        
                        for (int ir=-iDWin0; ir <= iDWin0; ir++) {
                            
                            int iindex = aiindex + ir;
                            int il= ail +ir;
                            
                            for (int ii=0; ii < iChannels; ii++)
                                fpODenoised[ii][iindex] +=  fWeight*fpI[ii][il];
                        }
                    }
                }
                
                /*SOMETIMES: select those patchs at distance < dmax*/
                for (kk=knnT;kk<dj*di;kk++)
                {
                    if (fDif_all[kk] < fDistance)
                    {
                        
                        fTotalWeight += fWeight;
                        
                        int i = ovect_ind[kk]/dj + imin;
                        int j = ovect_ind[kk]%dj + jmin;
                        
                        for (int is=-iDWin0; is <=iDWin0; is++) {
                            int aiindex = (iDWin+is) * ihwl + iDWin;
                            int ail = (j+is)*iWidth+i;
                            
                            for (int ir=-iDWin0; ir <= iDWin0; ir++) {
                                
                                int iindex = aiindex + ir;
                                int il= ail +ir;
                                
                                for (int ii=0; ii < iChannels; ii++)
                                    fpODenoised[ii][iindex] +=
                                    fWeight*fpI[ii][il];
                            }
                        }
                    }
                }
                
                //Normalize average value when fTotalweight is not near zero
                
                for (int is=-iDWin0; is <=iDWin0; is++) {
                    int aiindex = (iDWin+is) * ihwl + iDWin;
                    int ail=(y+is)*iWidth+x;
                    
                    for (int ir=-iDWin0; ir <= iDWin0; ir++) {
                        int iindex = aiindex + ir;
                        int il=ail+ ir;
                        
#pragma omp atomic
                        fpCount[il]++;
                        
                        for (int ii=0; ii < iChannels; ii++) {
#pragma omp atomic
                            
                            fpO[ii][il] += fpODenoised[ii][iindex]/fTotalWeight;
                            
                        }
                    }
                }
                
                delete[] ovect_ind;
                delete[] fDif_all;
                
            }//x loop end
            
            
            for (int ii=0; ii < iChannels; ii++) delete[] fpODenoised[ii];
            delete[] fpODenoised;
            
        }//yloop end
    }//pragma end
    
    for (int ii=0; ii < iwxh; ii++)
        if (fpCount[ii]>0.0) {
            for (int jj=0; jj < iChannels; jj++)  fpO[jj][ii] /= fpCount[ii];
            
        }       else {
            
            for (int jj=0; jj < iChannels; jj++)  fpO[jj][ii] = fpI[jj][ii];
        }
    
    // delete memory
    delete[] fpCount;
    
}


void rhf(int iDWin,        // Half size of patch
         int iDBloc,       // Half size of research window
         float fDistance,  // Max-Distance parameter
         float **fhI,      // Histogram
         float **fpI,      // Input
         float **fpO,      // Output
         int iChannels,    // Number of channels
         int iWidth,       // Image width
         int iHeight,      // Image height
         int iBins)        // Number of bins Histogram image
{
    
    printf("---->rhf: dmax = %f\n", fDistance);
    
    // length of each channel
    int iwxh = iWidth * iHeight;
    
    //  length of comparison window
    int ihwl = (2*iDWin+1);
    int iwl = (2*iDWin+1) * (2*iDWin+1);
    
    // auxiliary variable
    // number of denoised values per pixel
    float *fpCount = new float[iwxh];
    fpClear(fpCount, 0.0f,iwxh);
    
    // clear output
    for (int ii=0; ii < iChannels; ii++) fpClear(fpO[ii], 0.0f, iwxh);
    
    // PROCESS STARTS
    // for each pixel (x,y)
#pragma omp parallel shared(fpI, fpO)
    {
        
#pragma omp for schedule(dynamic) nowait
        
        for (int y=0; y < iHeight ; y++) {
            
            // auxiliary variable
            // denoised patch centered at a certain pixel
            float **fpODenoised = new float*[iChannels];
            for (int ii=0; ii < iChannels; ii++)
                fpODenoised[ii] = new float[iwl];
            
            for (int x=0 ; x < iWidth;  x++)
            {
                /*Check if we need to denoise this pixel!!*/
                // sum of weights
                float fTotalWeight = 0.0f;
                
                // weights
                float fWeight = 1.0f;
                
                //Reduce the size of comparison window near the boundary
                int iDWin0 = MIN(iDWin,MIN(iWidth-1-x,
                                           MIN(iHeight-1-y,MIN(x,y))));
                
                //Clear current denoised patch
                for (int ii=0; ii < iChannels; ii++)
                    fpClear(fpODenoised[ii], 0.0f, iwl);
                
                //Research zone depending on the boundary/size of the window
                int imin=MAX(x-iDBloc,iDWin0);
                int jmin=MAX(y-iDBloc,iDWin0);
                
                int imax=MIN(x+iDBloc,iWidth-1-iDWin0);
                int jmax=MIN(y+iDBloc,iHeight-1-iDWin0);
                
                for (int j=jmin; j <= jmax; j++)
                    for (int i=imin ; i <= imax; i++)
                        if (i!=x || j!=y)
                        {
                            
                            int df=0;
                            float fDifHist = fiChiSquareNDfFloatDist(&df,fhI,
                                                                     fhI,x,y,
                                                                     i,j,iDWin0,
                                                                     iBins,
                                                                     iWidth,
                                                                     iWidth);
                            if(fDifHist <  fDistance*df)
                            {
                                fTotalWeight += fWeight;
                                
                                for (int is=-iDWin0; is <=iDWin0; is++) {
                                    int aiindex = (iDWin+is) * ihwl + iDWin;
                                    int ail = (j+is)*iWidth+i;
                                    
                                    for (int ir=-iDWin0; ir <= iDWin0; ir++) {
                                        
                                        int iindex = aiindex + ir;
                                        int il= ail +ir;
                                        
                                        for (int ii=0; ii < iChannels; ii++)
                                            fpODenoised[ii][iindex] +=
                                            fWeight * fpI[ii][il];
                                    }
                                }
                            }
                        }
                
                // current patch with fMaxWeight
                for (int is=-iDWin0; is <=iDWin0; is++) {
                    int aiindex = (iDWin+is) * ihwl + iDWin;
                    int ail=(y+is)*iWidth+x;
                    
                    for (int ir=-iDWin0; ir <= iDWin0; ir++) {
                        
                        int iindex = aiindex + ir;
                        int il=ail+ir;
                        
                        for (int ii=0; ii < iChannels; ii++)
                            fpODenoised[ii][iindex] += fWeight * fpI[ii][il];
                    }
                }
                
                fTotalWeight += fWeight;
                
                // normalize average value when fTotalweight is not near zero
                if (fTotalWeight > fTiny) {
                    
                    for (int is=-iDWin0; is <=iDWin0; is++) {
                        int aiindex = (iDWin+is) * ihwl + iDWin;
                        int ail=(y+is)*iWidth+x;
                        
                        for (int ir=-iDWin0; ir <= iDWin0; ir++) {
                            int iindex = aiindex + ir;
                            int il=ail+ ir;
                            
#pragma omp atomic
                            fpCount[il]++;
                            
                            for (int ii=0; ii < iChannels; ii++) {
#pragma omp atomic
                                fpO[ii][il] +=
                                fpODenoised[ii][iindex]/fTotalWeight;
                            }
                        }
                    }
                }//end if Tiny
            }
            
            for (int ii=0; ii < iChannels; ii++) delete[] fpODenoised[ii];
            delete[] fpODenoised;
        }
    }
    
    for (int ii=0; ii < iwxh; ii++)
        if (fpCount[ii]>0.0) {
            for (int jj=0; jj < iChannels; jj++)  fpO[jj][ii] /= fpCount[ii];
            
        }       else {
            
            for (int jj=0; jj < iChannels; jj++)  fpO[jj][ii] = fpI[jj][ii];
        }
    
    // delete memory
    delete[] fpCount;
    
}
