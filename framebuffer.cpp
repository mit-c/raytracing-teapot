/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#include <stdlib.h>
#include <float.h>
#include <iostream>
#include <fstream>
using namespace std;

#include "framebuffer.h"
#include "colour.h"
#include <algorithm>
#include <list>

FrameBuffer::FrameBuffer(int w, int h)
{
  if ((w<0)||(w>2048)||(h<0)||(h>2048))
  {
    exit(-1);
  }

  this->framebuffer = new Pixel[w * h];

  if (!this->framebuffer)
  {
    exit(-2);
  }

  this->width = w;
  this->height = h;

  for (int i = 0; i < w*h; i += 1)
  {
    this->framebuffer[i].red   = 0.0f;
    this->framebuffer[i].green = 0.0f;
    this->framebuffer[i].blue  = 0.0f;
    this->framebuffer[i].depth = 0.0f;
  }
  
  return;
}

int FrameBuffer::plotPixel(int x, int y, float red, float green, float blue)
{
  if ((x<0)||(x>=this->width)||(y<0)||(y>=this->height))
  {
    return -1;
  }

  //  if ((red > 1.0f) || (red < 0.0f)) cerr<<"out of range\n";
  
  this->framebuffer[y * this->width + x].red = red;
  this->framebuffer[y * this->width + x].green = green;
  this->framebuffer[y * this->width + x].blue = blue;

  return 0;
}

int FrameBuffer::plotDepth(int x, int y, float depth)
{
  if ((x<0)||(x>=this->width)||(y<0)||(y>=this->height))
  {
    return -1;
  }

  this->framebuffer[y * this->width + x].depth = depth;

  return 0;
}

int FrameBuffer::getDepth(int x, int y, float &depth)
{
  if ((x<0)||(x>=this->width)||(y<0)||(y>=this->height))
  {
    return -1;
  }
  
  depth = this->framebuffer[y * this->width + x].depth;

  return 0;
}

int FrameBuffer::getPixel(int x, int y, float &red, float &green, float &blue)
{
  if ((x<0)||(x>=this->width)||(y<0)||(y>=this->height))
  {
    return -1;
  }

  red = this->framebuffer[y * this->width + x].red;
  green = this->framebuffer[y * this->width + x].green;
  blue = this->framebuffer[y * this->width + x].blue;

  return 0;
}


Colour FrameBuffer::getAverage(int x, int y, int window_size)
{
    // Given a valid x,y point in a framebuffer calculate average.
    if(window_size % 2 == 0)
    {
      std::cout << "window size should be even" << endl;
      return Colour();
    }
    // in case where window_size is 3
    // first pixel to get is x - window_size / 2. Note 3/2 = 1 for int div
   
    
    float r_sum=0;
    float g_sum=0;
    float b_sum=0;
    float r,g,b;
    int window_half = window_size/2;
    
    for(int i = x - window_half; i<x+window_half+1; i++)
    {
      for(int j = y - window_half; j<y+window_half+1; j++)
      {
        
        this->getPixel(i,j,r,g,b);
      
        r_sum += r;
        g_sum += g;
        b_sum += b;
      }
    }
    int n = window_size*window_size;
    float r_out = std::min(r_sum/(float)n,10.0f);
    float g_out = std::min(g_sum/(float)n,10.0f);
    float b_out = std::min(b_sum/(float)n,10.0f);
    
    return Colour(r_out,g_out,b_out,0);

}

int FrameBuffer::antiAlias(FrameBuffer &fb_out, int window_size)
{
 // This function outputs an anti-aliased version of the local framebuffer.
  int window_half = window_size/2; 
  int x_start = window_half;
  int y_start = x_start;
  int x_end = width - window_half;
  int y_end = height - window_half;
  int x_pixel,y_pixel;
  std::cout << "aliasing starts here:"  << endl;

  for(int x =x_start;x<x_end;x++)
  {
    for(int y = y_start;y<y_end;y++)
    {
      Colour colour_aa = this->getAverage(x,y,window_size);
      x_pixel = x - x_start;
      y_pixel = y - y_start;
      fb_out.plotPixel(x_pixel,y_pixel,colour_aa.r,colour_aa.g,colour_aa.b);

    }
  }
  return 0;
}

int FrameBuffer::writeRGBFile(char *filename)
{
  float min = 0.0f;
  float max = 0.0f;

  ofstream outfile;

  outfile.open(filename, ofstream::binary);

  if (!outfile.is_open())
  {
    return -1;
  }

  for (int i = 0; i < this->width*this->height; i += 1)
  {
    if (this->framebuffer[i].red > max) max = this->framebuffer[i].red;
    if (this->framebuffer[i].green > max) max = this->framebuffer[i].green;
    if (this->framebuffer[i].blue > max) max = this->framebuffer[i].blue;
    if (this->framebuffer[i].red < min) min = this->framebuffer[i].red;
    if (this->framebuffer[i].green < min) min = this->framebuffer[i].green;
    if (this->framebuffer[i].blue < min) min = this->framebuffer[i].blue;
  }
  std::cout <<"maximum intensity: " <<  max << endl;
  float diff = max - min;
  if (diff == 0.0f) diff = 1.0f;

  outfile << "P6\n";
  outfile << this->width << " " << this->height << "\n255\n";

  for (int j = 0; j<  this->width*this->height; j += 1)
  {
    outfile << (unsigned char)(((this->framebuffer[j].red-min)/diff)*255.0);
    outfile << (unsigned char)(((this->framebuffer[j].green)/diff)*255.0);
    outfile << (unsigned char)(((this->framebuffer[j].blue)/diff)*255.0);
  }
  
  outfile.close();
  return 0;
}



int FrameBuffer::writeDepthFile(char *filename)
{
  float max = 0;
  float min = FLT_MAX;

  ofstream outfile;

  outfile.open(filename, ofstream::binary);

  if (!outfile.is_open())
  {
    return -1;
  }

  for (int i = 0; i < this->width*this->height; i += 1)
  {
    if (this->framebuffer[i].depth > max) max = this->framebuffer[i].depth;
    if (this->framebuffer[i].depth < min) min = this->framebuffer[i].depth;
  }

  float diff = max - min;
  if (diff == 0.0f) diff = 1.0f;

  //  cerr << "Min/max/diff" << min << "/" << max << "/" << diff << endl;

  outfile << "P6\n";
  outfile << this->width << " " << this->height << "\n255\n";

  for (int j = 0; j<  this->width*this->height; j += 1)
  {
    unsigned char pd = (unsigned char)(((this->framebuffer[j].depth-min)/diff)*255.0);
    outfile << pd << pd << pd;
  }
  


  outfile.close();
  return 0;
}


