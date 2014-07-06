Ray Histogram Fusion
======================================================================
Version 1.1 - June 09, 2014

by    Mauricio Delbracio <mdelbra@gmail.com>	
      Pablo Muse
      Antoni Buades
      Jean-Michel Morel

Introduction
-----------
RHF is a multi-scale filter that accelerates Monte Carlo renderers.  
Each pixel in the image is characterized by the colors of the rays
that reach its surface. The RHF algorithm uses a statistical 
distance to compare with each other the ray color distributions 
associated with different pixels, at each scale. Based on this 
distance, it decides whether two pixels can share their rays or not.


Usage
-----

The modified PBRT-v2 code implements a film instance that
computes the samples color distribution. This is done by 
specifying the necessary variables (s,M,gamma,nbins)

Film "histo"
     "string filename" ["cornell-path_00256.exr"]
     "integer xresolution" [250] "integer yresolution" [250]
     "string histfilename" ["cornell-path_00256_hist.exr"]
     "integer nbins" [20]
     "float gamma" [2.2]
     "float M" [2.5]
     "float s" [2.0]


Copyright and License
---------------------
 
 RHF - Ray Histogram Fusion
 
 Copyright (c) 2014, A. Buades <toni.buades@uib.es>,
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
 
