#include "utils.h"

#include "integral.h"

IplImage *Integral(IplImage *source)
{

  IplImage *img = getGray(source);
  IplImage *int_img = cvCreateImage(cvGetSize(img), IPL_DEPTH_32F, 1);

  int height = img->height;
  int width = img->width;
  int step = img->widthStep/sizeof(float);
  float *data   = (float *) img->imageData;  
  float *i_data = (float *) int_img->imageData;  

  float rs = 0.0f;
  for(int j=0; j<width; j++) 
  {
    rs += data[j]; 
    i_data[j] = rs;
  }
  for(int i=1; i<height; ++i) 
  {
    rs = 0.0f;
    for(int j=0; j<width; ++j) 
    {
      rs += data[i*step+j]; 
      i_data[i*step+j] = rs + i_data[(i-1)*step+j];
    }
  }
  cvReleaseImage(&img);
  return int_img;
}

