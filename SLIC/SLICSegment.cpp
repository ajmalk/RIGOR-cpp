/*
 * SLICSegment.cpp
 *
 *  Created on: Apr 16, 2013
 *      Author: abhijit
 */

#include "SLICSegment.h"
#include "SLIC.h"
#include <iostream>

using namespace cv;
using namespace std;

// BGR + A  --- >  ARGB
Mat makeARGBfromBGR(const Mat& bgr) {
  Mat argb (bgr.rows, bgr.cols, CV_8UC4);
  Mat alpha = Mat::zeros( bgr.rows, bgr.cols, CV_8UC1);
  Mat in[] = {bgr,  alpha};


  int from_to[] = { 0,3, 1,2, 2,1, 3,0 };
  mixChannels( in, 2, &argb, 1, from_to, 4 );
  return argb;
}

// ARGB  --- >  BGR + A
Mat makeBGRfromARGB(const Mat& argb) {
  Mat bgr (argb.rows, argb.cols, CV_8UC3);
  Mat alpha = Mat::zeros( argb.rows, argb.cols, CV_8UC1);
  Mat out[] = {bgr,  alpha};


  int from_to[] = { 0,3, 1,2, 2,1, 3,0 };
  mixChannels( &argb, 1, out, 2, from_to, 4 );
  return bgr;
}

SLICSegment::SLICSegment(int req_num_of_super_pixels, double compactness)
    : req_num_of_super_pixels_(req_num_of_super_pixels),
      compactness_(compactness) {
}

int SLICSegment::performSegmentation(const cv::Mat& input, cv::Mat& result, const ResultType type) const {
  if (!input.data) {
      cout << "ERROR: Input Image is Empty\n";
      return 0;
    }

  Mat in_argb  = makeARGBfromBGR(input);

  unsigned int* pbuff = (unsigned int*) in_argb.data;

  SLIC segment;


  int* klabels = NULL;
  int numlabels(0);
  segment.DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(pbuff, in_argb.cols,
      in_argb.rows, klabels, numlabels, req_num_of_super_pixels_, compactness_);

  // Draw boundaries around segments
  segment.DrawContoursAroundSegments(pbuff, klabels, in_argb.cols, in_argb.rows, 0xff0000);

  if (type == LABELS)
    result = Mat(input.rows, input.cols, CV_32SC1, klabels);
  else
    result = makeBGRfromARGB(in_argb);

  return numlabels;
}
