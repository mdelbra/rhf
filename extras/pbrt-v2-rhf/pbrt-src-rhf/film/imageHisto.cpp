
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// film/image.cpp*
#include "stdafx.h"
#include "film/imageHisto.h"
#include "spectrum.h"
#include "parallel.h"
#include "imageio.h"

// ImageFilm Method Definitions

/*BEGIN ImageFilm Method Definitions by mdelbra 09.08.12 (histo)*/
ImageFilmHisto::ImageFilmHisto(int xres, int yres, Filter *filt, const float crop[4],
                     const string &fn, bool openWindow,
                     const string &histFn, int nB, float gam, float fMval, float fsval)
: Film(xres, yres) {
    filter = filt;
    memcpy(cropWindow, crop, 4 * sizeof(float));
    filename = fn;
    
    //Histogram parameters mdelbra 09.08.12
    histFilename = histFn;
    nBins = nB;
    gamma = gam;
    Mval = fMval;
    sval = fsval;
    
    // Compute film image extent
    xPixelStart = Ceil2Int(xResolution * cropWindow[0]);
    xPixelCount = max(1, Ceil2Int(xResolution * cropWindow[1]) - xPixelStart);
    yPixelStart = Ceil2Int(yResolution * cropWindow[2]);
    yPixelCount = max(1, Ceil2Int(yResolution * cropWindow[3]) - yPixelStart);
    
    // Allocate film image storage
    pixels = new BlockedArray<Pixel>(xPixelCount, yPixelCount);
    
    // Precompute filter weight table
#define FILTER_TABLE_SIZE 16
    filterTable = new float[FILTER_TABLE_SIZE * FILTER_TABLE_SIZE];
    float *ftp = filterTable;
    for (int y = 0; y < FILTER_TABLE_SIZE; ++y) {
        float fy = ((float)y + .5f) *
        filter->yWidth / FILTER_TABLE_SIZE;
        for (int x = 0; x < FILTER_TABLE_SIZE; ++x) {
            float fx = ((float)x + .5f) *
            filter->xWidth / FILTER_TABLE_SIZE;
            *ftp++ = filter->Evaluate(fx, fy);
        }
    }
    
    // Possibly open window for image display
    if (openWindow || PbrtOptions.openWindow) {
        Warning("Support for opening image display window not available in this build.");
    }
}
/*END ImageFilm Method Definitions by mdelbra 09.08.12 (histo)*/


/*END AddSampleHisto mdelbra 04.03.2013*/
void ImageFilmHisto::AddSample(const CameraSample &sample,
                          const Spectrum &L) {
    // Compute sample's raster extent
    float dimageX = sample.imageX - 0.5f;
    float dimageY = sample.imageY - 0.5f;
    int x0 = Ceil2Int (dimageX - filter->xWidth);
    int x1 = Floor2Int(dimageX + filter->xWidth);
    int y0 = Ceil2Int (dimageY - filter->yWidth);
    int y1 = Floor2Int(dimageY + filter->yWidth);
    x0 = max(x0, xPixelStart);
    x1 = min(x1, xPixelStart + xPixelCount - 1);
    y0 = max(y0, yPixelStart);
    y1 = min(y1, yPixelStart + yPixelCount - 1);
    if ((x1-x0) < 0 || (y1-y0) < 0)
    {
        PBRT_SAMPLE_OUTSIDE_IMAGE_EXTENT(const_cast<CameraSample *>(&sample));
        return;
    }
    
    // Loop over filter support and add sample to pixel arrays
    float xyz[3];
    L.ToXYZ(xyz);
    
    float rgb[3];
    L.ToRGB(rgb); //mdelbra 09.08.12
    
    // Precompute $x$ and $y$ filter table offsets
    int *ifx = ALLOCA(int, x1 - x0 + 1);
    for (int x = x0; x <= x1; ++x) {
        float fx = fabsf((x - dimageX) *
                         filter->invXWidth * FILTER_TABLE_SIZE);
        ifx[x-x0] = min(Floor2Int(fx), FILTER_TABLE_SIZE-1);
    }
    int *ify = ALLOCA(int, y1 - y0 + 1);
    for (int y = y0; y <= y1; ++y) {
        float fy = fabsf((y - dimageY) *
                         filter->invYWidth * FILTER_TABLE_SIZE);
        ify[y-y0] = min(Floor2Int(fy), FILTER_TABLE_SIZE-1);
    }
    bool syncNeeded = (filter->xWidth > 0.5f || filter->yWidth > 0.5f);
    for (int y = y0; y <= y1; ++y) {
        for (int x = x0; x <= x1; ++x) {
            // Evaluate filter value at $(x,y)$ pixel
            int offset = ify[y-y0]*FILTER_TABLE_SIZE + ifx[x-x0];
            float filterWt = filterTable[offset];
            
            
            // Update pixel values with filtered sample contribution
            Pixel &pixel = (*pixels)(x - xPixelStart, y - yPixelStart);
            if (!syncNeeded) {
                pixel.Lxyz[0] += filterWt * xyz[0];
                pixel.Lxyz[1] += filterWt * xyz[1];
                pixel.Lxyz[2] += filterWt * xyz[2];
                pixel.weightSum += filterWt;
                
                //Add to histogram contribution mdelbra 09.08.12
                pixel.nsamples += filterWt;
                
                int ibL;
                int ibH;
                float wbL;
                float wbH;
                
                
                for(int j=0;j<3;j++)//Loop XYZ or RGB
                {
                    
                    float v = filterWt*rgb[j];
                    v = v>0?v:0;
                    
                    /*Compress dynamical range*/
                    if(gamma>1) v = pow(v,1/gamma);
                    
                    /*normalize to max_Vale*/
                    if(Mval>0) v  = v/Mval;
                    
                    /*Truncate to SATURE_LEVEL*/
                    v = v> sval ? sval : v;
                    
                    float fbin = v * (nBins-2);
                    
                    int ibinL = (int) (fbin); //Low bin 
                    
                    /*Check out of bounds when ibinL > nbins-1; 
                     only equality is possible?
                     */
                    
                    if(ibinL < nBins-2) //inbounds
                    {
                        //High bin weight, Low bin weight 1-wH
                        
                        float wH = fbin - ibinL;
                        ibL = ibinL;
                        wbL = 1.0f - wH;
                        ibH = ibinL+1;
                        wbH = wH;
                    }
                    else { //out of bounds... v >= 1
                        
                        float wH = (v - 1.0f)/(sval - 1);
                        
                        ibL = nBins-2;
                        wbL = 1.0f - wH;
                        ibH = nBins-1;
                        wbH = wH;
                    }
                    
                    pixel.histLxyz[j][ibL] += wbL;
                    pixel.histLxyz[j][ibH] += wbH;
                    
                    
                }
                
                
                
            }
            
            else {
                // Safely update _Lxyz_ and _weightSum_ even with concurrency
                AtomicAdd(&pixel.Lxyz[0], filterWt * xyz[0]);
                AtomicAdd(&pixel.Lxyz[1], filterWt * xyz[1]);
                AtomicAdd(&pixel.Lxyz[2], filterWt * xyz[2]);
                AtomicAdd(&pixel.weightSum, filterWt);
                
                //Add to histogram contribution
                
                AtomicAdd(&pixel.nsamples, filterWt);
                //AtomicAdd(&pixel.nsamples, 1);
                
                int ibL;
                int ibH;
                float wbL;
                float wbH;
                
                
                for(int j=0;j<3;j++)
                {
                    
                    float v = filterWt*rgb[j];
                    v = v>0?v:0;
                    
                    /*Compress dynamical range*/
                    if(gamma>1) v = pow(v,1/gamma);
                    
                    /*normalize to max_Vale*/
                    if(Mval>0) v  = v/Mval;
                    
                    /*Truncate to SATURE_LEVEL*/
                    v = v>sval ? sval : v;
                    
                    float fbin = v * (nBins-2);
                    
                    int ibinL = (int) (fbin); //Low bin 
                    
                    /*Check out of bounds when ibinL > nbins-1; 
                     only equality is possible?
                     */
                    
                    if(ibinL < nBins-2) //inbounds
                    {
                        //High bin weight, Low bin weight 1-wH
                        
                        float wH = fbin - ibinL;
                        ibL = ibinL;
                        wbL = 1.0f - wH;
                        ibH = ibinL+1;
                        wbH = wH;
                    }
                    else { //out of bounds... v >= 1
                        
                        float wH = (v - 1.0f)/(sval - 1);
                        
                        ibL = nBins-2;
                        wbL = 1.0f - wH;
                        ibH = nBins-1;
                        wbH = wH;
                    }
                    
                    AtomicAdd(&pixel.histLxyz[j][ibL], wbL);
                    AtomicAdd(&pixel.histLxyz[j][ibH], wbH);
                    
                }
                
            }
            
        }
    }
}
/*END AddSampleHisto mdelbra 04.03.2013*/

void ImageFilmHisto::Splat(const CameraSample &sample, const Spectrum &L) {
    if (L.HasNaNs()) {
        Warning("ImageFilm ignoring splatted spectrum with NaN values");
        return;
    }
    float xyz[3];
    L.ToXYZ(xyz);
    int x = Floor2Int(sample.imageX), y = Floor2Int(sample.imageY);
    if (x < xPixelStart || x - xPixelStart >= xPixelCount ||
        y < yPixelStart || y - yPixelStart >= yPixelCount) return;
    Pixel &pixel = (*pixels)(x - xPixelStart, y - yPixelStart);
    AtomicAdd(&pixel.splatXYZ[0], xyz[0]);
    AtomicAdd(&pixel.splatXYZ[1], xyz[1]);
    AtomicAdd(&pixel.splatXYZ[2], xyz[2]);
}


void ImageFilmHisto::GetSampleExtent(int *xstart, int *xend,
                                int *ystart, int *yend) const {
    *xstart = Floor2Int(xPixelStart + 0.5f - filter->xWidth);
    *xend   = Ceil2Int(xPixelStart + 0.5f + xPixelCount +
                       filter->xWidth);

    *ystart = Floor2Int(yPixelStart + 0.5f - filter->yWidth);
    *yend   = Ceil2Int(yPixelStart + 0.5f + yPixelCount +
                       filter->yWidth);
}


void ImageFilmHisto::GetPixelExtent(int *xstart, int *xend,
                               int *ystart, int *yend) const {
    *xstart = xPixelStart;
    *xend   = xPixelStart + xPixelCount;
    *ystart = yPixelStart;
    *yend   = yPixelStart + yPixelCount;
}



/*BEGIN WriteImageHisto mdelbra 05.03.2013*/
void ImageFilmHisto::WriteImage(float splatScale) {
    // Convert image to RGB and compute final pixel values
    int nPix = xPixelCount * yPixelCount;
    float *rgb = new float[3*nPix];
    
    //3*nbins + nsamples per pixel
    float *hist = new float[(3*nBins+1)*nPix];
    
    int offset = 0;
    for (int y = 0; y < yPixelCount; ++y) {
        for (int x = 0; x < xPixelCount; ++x) {
            // Convert pixel XYZ color to RGB
            XYZToRGB((*pixels)(x, y).Lxyz, &rgb[3*offset]);
            
            // Normalize pixel with weight sum
            float weightSum = (*pixels)(x, y).weightSum;
            if (weightSum != 0.f) {
                /*if (weightSum > 0.001) {*/
                float invWt = 1.f / weightSum;
                //if(rgb[3*offset  ]<0 || rgb[3*offset+1]<0 || rgb[3*offset+2 ]<0) 
                //    printf("%d %d\n", x,y);
                rgb[3*offset  ] = max(0.f, rgb[3*offset  ] * invWt);
                rgb[3*offset+1] = max(0.f, rgb[3*offset+1] * invWt);
                rgb[3*offset+2] = max(0.f, rgb[3*offset+2] * invWt);
            }
            
            
            float splatRGB[3];
            XYZToRGB((*pixels)(x, y).splatXYZ, splatRGB);
            rgb[3*offset  ] += splatScale * splatRGB[0];
            rgb[3*offset+1] += splatScale * splatRGB[1];
            rgb[3*offset+2] += splatScale * splatRGB[2];
            
            //Save histogram mdelbra 09.08.12
            for(int i=0;i<nBins;i++)
                for(int j=0;j<3;j++)
                    hist[j*nPix*nBins + i*nPix + offset] = (*pixels)(x,y).histLxyz[j][i];
            
            hist[3*nPix*nBins + offset] = (*pixels)(x,y).nsamples;
            
            ++offset;
        }
    }
    
    // Write RGB image
    ::WriteImage(filename, rgb, NULL, xPixelCount, yPixelCount,
                 xResolution, yResolution, xPixelStart, yPixelStart);
    
    // Write Histogram image
    ::WriteMultiChannelImage(histFilename, hist, NULL, xPixelCount, yPixelCount,
                             xResolution, yResolution, xPixelStart, yPixelStart, 3*nBins + 1);
    
    
    
    // Release temporary image memory
    delete[] rgb;
}
/*END WriteImageHisto mdelbra 05.03.2013*/


void ImageFilmHisto::UpdateDisplay(int x0, int y0, int x1, int y1,
    float splatScale) {
}

/*BEGIN CreateImagefilm mdelbra 05.04.2013*/
ImageFilmHisto *CreateImageFilmHisto(const ParamSet &params, Filter *filter) {
    string filename = params.FindOneString("filename", "");
        
    /*Read histogram Parameters mdelbra 05.04.2013*/
    string histfilename = params.FindOneString("histfilename", "");
    int    iBins   = params.FindOneInt("nbins", 20);
    float  fgamma  = params.FindOneFloat("gamma", 2.2f);
    float  fMval = params.FindOneFloat("M", 2.5f);
    float  fsval = params.FindOneFloat("s", 2.0f);

    
    if (PbrtOptions.imageFile != "") {
        if (filename != "") {
            Warning("Output filename supplied on command line, \"%s\", ignored "
                    "due to filename provided in scene description file, \"%s\".",
                    PbrtOptions.imageFile.c_str(), filename.c_str());
        }
        else
            filename = PbrtOptions.imageFile;
    }
    if (filename == "")
#ifdef PBRT_HAS_OPENEXR
        filename = "pbrt.exr";
#else
        filename = "pbrt.tga";
#endif

    if (histfilename == "")
#ifdef PBRT_HAS_OPENEXR
        histfilename = "pbrtHist.exr";
#else
        Error("EXR is needed to work in Histogram Mode\n");
#endif  
    
    
    int xres = params.FindOneInt("xresolution", 640);
    int yres = params.FindOneInt("yresolution", 480);
    if (PbrtOptions.quickRender) xres = max(1, xres / 4);
    if (PbrtOptions.quickRender) yres = max(1, yres / 4);
    bool openwin = params.FindOneBool("display", false);
    float crop[4] = { 0, 1, 0, 1 };
    int cwi;
    const float *cr = params.FindFloat("cropwindow", &cwi);
    if (cr && cwi == 4) {
        crop[0] = Clamp(min(cr[0], cr[1]), 0., 1.);
        crop[1] = Clamp(max(cr[0], cr[1]), 0., 1.);
        crop[2] = Clamp(min(cr[2], cr[3]), 0., 1.);
        crop[3] = Clamp(max(cr[2], cr[3]), 0., 1.);
    }

    //Mode HISTO? mdelbra 11.06.12
    return new ImageFilmHisto(xres, yres, filter, crop, filename, openwin,
                             histfilename, iBins, fgamma, fMval, fsval);
        
}


