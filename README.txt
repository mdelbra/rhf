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

Files
-----

COPYING
Makefile
README.txt
VERSION
exrcrop.cpp
exrdiff.cpp
exrtopng.cpp
io_exr.cpp
io_exr.h
io_png.c
io_png.h
libauxiliar.cpp
libauxiliar.h
libdenoising.cpp
libdenoising.h
rhf.cpp
extras/pbrt-v2-rhf (A modified version of PBRT-v2)


Requirements
------------
- OpenEXR for reading EXR images

Compilation
-----------
Simply use the provided makefile, with the command `make`. You need to set
the directory where the libraries: openEXR have the respective
header and libraries files.

Running
-------

Usage: ./rhf [options] <input file> <output file>
Only EXR images are supported.

Options:
   -h <hist>   The filename with the histogram
   -d <float>  Max-distance between patchs
   -k <int>    Minimum number of similar patchs (default: 2)
   -b <int>    Half the block size  (default: 6)
   -w <int>    Half the windows size (default: 1)
   -s <int>    Number of Scales - Multi-Scale (default: 2)


Please report bugs in RHF to <mdelbra@gmail.com>.


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
 