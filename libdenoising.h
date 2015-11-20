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

#ifndef _LIBDENOISING_H_
#define _LIBDENOISING_H_

#include "libauxiliar.h"

/** @brief RHF filter function for the case the user sets a
 *         minimum number of similar patches (k), single scale
 *
 *  @param iDWin half size of patch
 *  @param iDBloc half size of research window
 *  @param fDistance max-distance threshold
 *  @param knn minimum number of similar patchs
 *  @param fhI histogram
 *  @param fpI input image
 *  @param fpO output (filtered) image
 *  @param iChannels number of color channels (should be 3)
 *  @param iWidth  image width
 *  @param iHeight image height
 *  @param iBins number of bins per pixel in histogram
 *  @return void
 */
void rhf_knn(int iDWin,            // Half size of patch
             int iDBloc,           // Half size of research window
             float fDistance,      // Max-Distance parameter
             int knn,
             float **fhI,          // Histogram
             float **fpI,          // Input
             float **fpO,          // Output
             float *alpha,         // Alpha channel, can be NULL
             int iChannels, int iWidth,int iHeight, int iBins);


/** @brief RHF filter function for the case the user does not set a
 *         minimum number of similar patches, single scale
 *
 *  @param iDWin half size of patch
 *  @param iDBloc half size of research window
 *  @param fDistance max-distance threshold
 *  @param fhI histogram
 *  @param fpI input image
 *  @param fpO output (filtered) image
 *  @param iChannels number of color channels (should be 3)
 *  @param iWidth  image width
 *  @param iHeight image height
 *  @param iBins number of bins per pixel in histogram
 *  @return void
 */

void rhf(int iDWin,            // Half size of patch
         int iDBloc,           // Half size of research window
         float fDistance,      // Max-Distance parameter
         float **fhI,          // Histogram
         float **fpI,          // Input
         float **fpO,          // Output
         float *alpha,         // Alpha channel, can be NULL
         int iChannels, int iWidth,int iHeight, int iBins);


/** @brief Multi-scale RHF filter function for the case the user sets a
 *         minimum number of similar patches (k)
 *
 *  @param iDWin half size of patch
 *  @param iDBloc half size of research window
 *  @param fDistance max-distance threshold
 *  @param knn minimum number of similar patchs
 *  @param iNscales number of scales
 *  @param fhI histogram
 *  @param fpI input image
 *  @param ffpO output (filtered) image
 *  @param iChannels number of color channels (should be 3)
 *  @param iWidth  image width
 *  @param iHeight image height
 *  @param iBins number of bins per pixel in histogram
 *  @return void
 */
void rhf_multiscale(int iDWin,            // Half size of patch
                    int iDBloc,           // Half size of research window
                    float fDistance,      // Max-Distance parameter
                    int knn,
                    int iNscales,         // Number of Scales
                    float **fhI,          // Histogram
                    float **fpI,          // Input
                    float **fpO,          // Output
                    float *alpha,         // Alpha channel, can be NULL
                    int iChannels, int iWidth,int iHeight, int iBins);



#endif
