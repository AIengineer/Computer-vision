
#ifndef UTILS_H
#define UTILS_H

#include <cv.h>
#include "ipoint.h"

#include <vector>
void error(const char *msg);
void showImage(const IplImage *img);
void showImage(char *title,const IplImage *img);
IplImage* getGray(const IplImage *img);
void drawIpoint(IplImage *img, Ipoint &ipt, int tailSize = 0);
void drawIpoints(IplImage *img, std::vector<Ipoint> &ipts, int tailSize = 0);

void drawWindows(IplImage *img, std::vector<Ipoint> &ipts);
void drawFPS(IplImage *img);
void drawPoint(IplImage *img, Ipoint &ipt);
void drawPoints(IplImage *img, std::vector<Ipoint> &ipts);
void saveSurf(char *filename, std::vector<Ipoint> &ipts);
void loadSurf(char *filename, std::vector<Ipoint> &ipts);
inline int fRound(float flt)
{
  return (int) floor(flt+0.5f);
}

#endif
