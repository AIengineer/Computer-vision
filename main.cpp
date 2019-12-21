#include "surflib.h"
#include "kmeans.h"
#include <ctime>
#include <iostream>
int mainImage(void)
{
  // Declare Ipoints and other stuff
  IpVec ipts;
  IplImage *img=cvLoadImage("imgs/sf.jpg");

  // Detect and describe interest points in the image
  clock_t start = clock();
  surfDetDes(img, ipts, false, 5, 4, 2, 0.0004f); 
  clock_t end = clock();

  std::cout<< "OpenSURF found: " << ipts.size() << " interest points" << std::endl;
  std::cout<< "OpenSURF took: " << float(end - start) / CLOCKS_PER_SEC  << " seconds" << std::endl;

  // Draw the detected points
  drawIpoints(img, ipts);
  
  // Display the result
  showImage(img);

  return 0;
}
int main(void) 
{
	mainImage();
}
