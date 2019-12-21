#ifndef IPOINT_H
#define IPOINT_H

#include <vector>
#include <math.h>

class Ipoint; 
typedef std::vector<Ipoint> IpVec;
typedef std::vector<std::pair<Ipoint, Ipoint> > IpPairVec;

void getMatches(IpVec &ipts1, IpVec &ipts2, IpPairVec &matches);
int translateCorners(IpPairVec &matches, const CvPoint src_corners[4], CvPoint dst_corners[4]);

class Ipoint {

public:

  ~Ipoint() {};

  Ipoint() : orientation(0) {};
  float operator-(const Ipoint &rhs)
  {
    float sum=0.f;
    for(int i=0; i < 64; ++i)
      sum += (this->descriptor[i] - rhs.descriptor[i])*(this->descriptor[i] - rhs.descriptor[i]);
    return sqrt(sum);
  };

  float x, y;

  float scale;

  float orientation;

  int laplacian;

  float descriptor[64];

  float dx, dy;

  int clusterIndex;
};
#endif
