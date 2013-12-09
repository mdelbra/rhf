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




void rhf_knn(int iDWin,            // Half size of patch
                         int iDBloc,           // Half size of research window
                         float fDistance,      // Max-Distance parameter
                         int knn,
                         float **fhI,          // Histogram
                         float **fpI,          // Input
                         float **fpO,          // Output
                         int iChannels, int iWidth,int iHeight, int iBins);    


void rhf(int iDWin,            // Half size of patch
                             int iDBloc,           // Half size of research window
                             float fDistance,      // Max-Distance parameter
                             float **fhI,          // Histogram
                             float **fpI,          // Input
                             float **fpO,          // Output
                             int iChannels, int iWidth,int iHeight, int iBins);    




void rhf_multiscale(int iDWin,            // Half size of patch
                            int iDBloc,           // Half size of research window
                            float fDistance,      // Max-Distance parameter
                            int knn,    
                            int iNscales,         // Number of Scales
                            float **fhI,          // Histogram
                            float **fpI,          // Input
                            float **fpO,          // Output
                            int iChannels, int iWidth,int iHeight, int iBins);

    

#endif
