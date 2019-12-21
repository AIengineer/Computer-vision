#ifndef SURF_H
#define SURF_H

#include <cv.h>
#include "ipoint.h"
#include "integral.h"

#include <vector>

class Surf {
  
  public:
    Surf(IplImage *img, std::vector<Ipoint> &ipts);

    void getDescriptors(bool bUpright = false);
  
  private:


    void getOrientation();
    

    void getDescriptor(bool bUpright = false);

    inline float gaussian(int x, int y, float sig);
    inline float gaussian(float x, float y, float sig);


    inline float haarX(int row, int column, int size);
    inline float haarY(int row, int column, int size);

    float getAngle(float X, float Y);

    IplImage *img;

    IpVec &ipts;

    int index;
};


#endif
