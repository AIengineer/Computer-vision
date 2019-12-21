
#ifndef SURFLIB_H
#define SURFLIB_H

#include <cv.h>
#include <highgui.h>

#include "integral.h"
#include "fasthessian.h"
#include "surf.h"
#include "ipoint.h"
#include "utils.h"

inline void surfDetDes(IplImage *img,  
                       std::vector<Ipoint> &ipts, 
                       bool upright = false, 
                       int octaves = OCTAVES, 
                       int intervals = INTERVALS, 
                       int init_sample = INIT_SAMPLE,
                       float thres = THRES )
{

  IplImage *int_img = Integral(img);

  FastHessian fh(int_img, ipts, octaves, intervals, init_sample, thres);

  fh.getIpoints();

  Surf des(int_img, ipts);

  des.getDescriptors(upright);
  cvReleaseImage(&int_img);
}

inline void surfDet(IplImage *img,  
                    std::vector<Ipoint> &ipts,
                    int octaves = OCTAVES, 
                    int intervals = INTERVALS,
                    int init_sample = INIT_SAMPLE,
                    float thres = THRES )
{

  IplImage *int_img = Integral(img);

  FastHessian fh(int_img, ipts, octaves, intervals, init_sample, thres);

  fh.getIpoints();

  cvReleaseImage(&int_img);
}


inline void surfDes(IplImage *img, 
                    std::vector<Ipoint> &ipts, 
                    bool upright = false)
{ 

  IplImage *int_img = Integral(img);

  Surf des(int_img, ipts);

  des.getDescriptors(upright);

  cvReleaseImage(&int_img);
}


#endif
