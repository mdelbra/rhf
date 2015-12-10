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

#include <ImfInputFile.h>
#include <ImfOutputFile.h>
#include <ImfRgbaFile.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>
#include <half.h>
#include <stdlib.h>
#include <stdio.h>


using namespace Imf;
using namespace Imath;

// EXR Function Definitions
struct sRGB
{
    float r, g, b;
};


float  *ReadImageEXR(const char fileName[], Box2i &dataWindow, Box2i &displayWindow, float *&alpha)
{
    
    
    try {
        // Open EXR file
        
        InputFile file (fileName);
        
        // Get image dimensions.
        dataWindow = file.header().dataWindow();
        displayWindow = file.header().displayWindow();
        
        int width = dataWindow.max.x - dataWindow.min.x + 1;
        int height = dataWindow.max.y - dataWindow.min.y + 1;
        
        const bool hasAlpha = file.header().channels().findChannel("A") != NULL;
        
        // Allocate memory to read image bits. We will only try to read R, G and B
        // here, but additional channels like A (alpha) could also be added...
        float *pixelsR = new float[width * height];
        float *pixelsG = new float[width * height];
        float *pixelsB = new float[width * height];
        alpha = hasAlpha ? new float[width * height] : NULL;
        
        // Now create the frame buffer to feed the image reader with. We will use
        // the Slice method flexibility to directly read R, G and B data in an
        // interlaced manner, using appropriate x & y stride values.
        FrameBuffer frameBuffer;
        
        frameBuffer.insert("R", Slice(FLOAT, (char*)(&pixelsR[0] -
                                                     dataWindow.min.x - dataWindow.min.y*width),
                                      sizeof(float), 
                                      width * sizeof(float), 
                                      1, 1, // x/y sampling
                                      0.0));
        frameBuffer.insert("G", Slice(FLOAT, (char*)(&pixelsG[0] -
                                                     dataWindow.min.x -dataWindow.min.y*width),
                                      sizeof(float), 
                                      width * sizeof(float),
                                      1, 1, // x/y sampling
                                      0.0));
        frameBuffer.insert("B", Slice(FLOAT, (char*)(&pixelsB[0] -
                                                     dataWindow.min.x -dataWindow.min.y*width),
                                      sizeof(float), 
                                      width * sizeof(float),
                                      1, 1, // x/y sampling
                                      0.0));
        if (alpha)
            frameBuffer.insert("A", Slice(FLOAT, (char*)(alpha -
                                                     dataWindow.min.x -dataWindow.min.y*width),
                                      sizeof(float), 
                                      width * sizeof(float),
                                      1, 1, // x/y sampling
                                      0.0));
        
        file.setFrameBuffer(frameBuffer);
        
        // Saving images
        float *data = new float[3*width * height];
        
        
        try {
            file.readPixels (dataWindow.min.y, dataWindow.max.y);
            
        }
        catch (const std::exception &) {
            data = NULL;
        }
        
        
        for (int i=0; i < width * height; i++) 
        {
            data[i] =  pixelsR[i];
            data[i+width*height] =  pixelsG[i];
            data[i+2*width*height] =  pixelsB[i];
        }
        
        
        delete[] pixelsR;
        delete[] pixelsG;
        delete[] pixelsB;
        
        return data;
    }
    
    catch (const std::exception &) {
        printf("Error reading file: %s\n",fileName);
        exit(-1);
    }
    
}

float  *ReadImageEXR(const char fileName[], int *nx, int *ny)
{
    Box2i dataWindow, displayWindow;
    float *alpha;
    float *rgb = ReadImageEXR (fileName, dataWindow, displayWindow, alpha);
    if (alpha)
        delete[] alpha;
    *nx = dataWindow.max.x - dataWindow.min.x + 1;
    *ny = dataWindow.max.y - dataWindow.min.y + 1;
    return rgb;
}

void WriteImageEXR(const char *name, float *pixels, float *alpha,
                    const Imath::Box2i &dataWindow,
                    const Imath::Box2i &displayWindow, int channelStride)
{
    const int xRes = dataWindow.max.x - dataWindow.min.x + 1;
    const int yRes = dataWindow.max.y - dataWindow.min.y + 1;
    
    Rgba *hrgba = new Rgba[xRes * yRes];
    for (int i = 0; i < xRes * yRes; ++i)
    {

        hrgba[i] = Rgba(pixels[i], pixels[i+channelStride],pixels[i+2*channelStride],
                        alpha ? alpha[i]: 1.f);
                        
    }
    
    try {
        RgbaOutputFile file(name, displayWindow, dataWindow, WRITE_RGBA);
        file.setFrameBuffer(hrgba - dataWindow.min.x - dataWindow.min.y * xRes, 1, xRes);
        file.writePixels(yRes);
    }
    catch (const std::exception &e) {
        printf("Unable to write image file \"%s\": %s", name, e.what());
        exit(-1);
    }
    
    delete[] hrgba;
}

void WriteImageEXR(const char *name, float **pixels,
                   int xRes, int yRes)
{
    
    //this can be a parameter to gerneralize the function
    int xOffset = 0;
    int yOffset = 0;
    int totalXRes = xRes;
    int totalYRes = yRes;
    
    float *alpha = NULL;
    
    Rgba *hrgba = new Rgba[xRes * yRes];
    for (int i = 0; i < xRes * yRes; ++i)
    {

        hrgba[i] = Rgba(pixels[0][i], pixels[1][i],pixels[2][i],
                        alpha ? alpha[i]: 1.f);
        
    }
    
    Box2i displayWindow(V2i(0,0), V2i(totalXRes-1, totalYRes-1));
    Box2i dataWindow(V2i(xOffset, yOffset), 
                     V2i(xOffset + xRes - 1, yOffset + yRes - 1));
    
    try {
        RgbaOutputFile file(name, displayWindow, dataWindow, WRITE_RGBA);
        file.setFrameBuffer(hrgba - xOffset - yOffset * xRes, 1, xRes);
        file.writePixels(yRes);
    }
    catch (const std::exception &e) {
        printf("Unable to write image file \"%s\": %s", name, e.what());
        exit(-1);
    }
    
    delete[] hrgba;
}

void WriteImageEXR(const char *name, float *pixels, int nx, int ny)
{
    Imath::Box2i dw;
    dw.min.x = dw.min.y = 0;
    dw.max.x = nx-1;
    dw.max.y = ny-1;
    WriteImageEXR(name, pixels, NULL, dw, dw, nx*ny);
}

void writeMultiImageEXR (const char *fileName,
          float *zPixels,
          int width,
          int height,
          int nchan)
{
    
    int nhnc = height*width;
    
    Header header (width, height); 
    
    for(int i =0; i<nchan;i++)
    {
        char ch_name[10];
        sprintf(ch_name,"Bin_%04d",i);
        header.channels().insert (ch_name, Channel (FLOAT)); 
    
    }
    
    OutputFile file(fileName, header); 
    
    FrameBuffer frameBuffer;

    for(int i =0; i<nchan;i++)
    {
        char ch_name[10];
        sprintf(ch_name,"Bin_%04d",i);
        frameBuffer.insert (ch_name, // name 
                            Slice (FLOAT, // type 
                                   (char *) &zPixels[i*nhnc], // base 
                                   sizeof (*zPixels) * 1, // xStride
                                   sizeof (*zPixels) * width)); // yStride
        
        
        
    }
    
       
    file.setFrameBuffer (frameBuffer); 
    file.writePixels (height); 

}

void writeMultiImageEXR (const char *fileName,
                         float **zPixels,
                         int width,
                         int height,
                         int nchan)
{
        
    Header header (width, height); 
    
    for(int i =0; i<nchan;i++)
    {
        char ch_name[10];
        sprintf(ch_name,"Bin_%04d",i);
        header.channels().insert (ch_name, Channel (FLOAT)); 
        
    }
    
    OutputFile file(fileName, header); 
    
    FrameBuffer frameBuffer;
    
    for(int i =0; i<nchan;i++)
    {
        char ch_name[10];
        sprintf(ch_name,"Bin_%04d",i);
        frameBuffer.insert (ch_name, // name 
                            Slice (FLOAT, // type 
                                   (char *) zPixels[i], // base 
                                   sizeof (*zPixels) * 1, // xStride
                                   sizeof (*zPixels) * width)); // yStride
        
        
        
    }
    
    
    file.setFrameBuffer (frameBuffer); 
    file.writePixels (height); 
    
}

float * readMultiImageEXR(const char fileName[],
                          int *width, int *height, int *nbins)
{
    
    float *data;
    
    try {
        InputFile file (fileName);
        
        Box2i dw = file.header().dataWindow();
        
        *width = dw.max.x - dw.min.x + 1;
        *height = dw.max.y - dw.min.y + 1;
        int nhnc = (*width) * (*height);
        
        const ChannelList &channelList = file.header().channels();
        
        FrameBuffer frameBuffer;
        
        char ch_name[10];
        int nbin =0;
        sprintf(ch_name,"Bin_%04d",nbin);

        const Channel *channelPtr = channelList.findChannel(ch_name);
        
        
        while(channelPtr)
        {
            nbin++;
            sprintf(ch_name,"Bin_%04d",nbin);
            // printf("%s\n",ch_name);
            channelPtr = channelList.findChannel(ch_name);
        }
        
        
        //reserve memory for the whole array
        data = new float[nhnc*nbin];
        
        
        for(int i=0;i<nbin;i++)
        {        
            sprintf(ch_name,"Bin_%04d",i);
            
            frameBuffer.insert(ch_name, 
                               Slice(FLOAT, 
                                     (char*)(&data[i*nhnc] -
                                     dw.min.x -dw.min.y**width), 
                                     sizeof(float), 
                                     (*width) * sizeof(float)));
            
        }
        
        
        file.setFrameBuffer (frameBuffer);
        
        try {
            file.readPixels (dw.min.y, dw.max.y);
            
        }
        catch (const std::exception &) {
            data = NULL;
        }
        
        
        *nbins = nbin;
        return data;
        
    }
    catch (const std::exception &) {
        printf("Error reading file: %s\n",fileName);
        exit(-1);
    }
    
}
