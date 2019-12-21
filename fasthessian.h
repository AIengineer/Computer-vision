#ifndef FASTHESSIAN_H
#define FASTHESSIAN_H

#include <cv.h>
#include "ipoint.h"

#include <vector>

class ResponseLayer;
static const int OCTAVES = 5;
static const int INTERVALS = 4;
static const float THRES = 0.0004f;
static const int INIT_SAMPLE = 2;


class FastHessian {
  
  public:
   
    FastHessian(std::vector<Ipoint> &ipts, 
                const int octaves = OCTAVES, 
                const int intervals = INTERVALS, 
                const int init_sample = INIT_SAMPLE, 
                const float thres = THRES);

    FastHessian(IplImage *img, 
                std::vector<Ipoint> &ipts, 
                const int octaves = OCTAVES, 
                const int intervals = INTERVALS, 
                const int init_sample = INIT_SAMPLE, 
                const float thres = THRES);

    ~FastHessian();

    void saveParameters(const int octaves, 
                        const int intervals,
                        const int init_sample, 
                        const float thres);

    void setIntImage(IplImage *img);

    void getIpoints();
    
  private:

    void buildResponseMap();
    void buildResponseLayer(ResponseLayer *r);

    int isExtremum(int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b);    
    
    void interpolateExtremum(int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b);
    void interpolateStep(int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b,
                          double* xi, double* xr, double* xc );
    CvMat* deriv3D(int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b);
    CvMat* hessian3D(int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b);

    IplImage *img;
    int i_width, i_height;

    std::vector<Ipoint> &ipts;

    std::vector<ResponseLayer *> responseMap;

    int octaves;

    int intervals;

    int init_sample;

    float thresh;
};


#endif
