/*
 * SLICSegment.h
 *
 *  Created on: Apr 16, 2013
 *      Author: abhijit
 */

#ifndef SLIC_SEGMENT_H_
#define SLIC_SEGMENT_H_

#include <opencv2/core/core.hpp>

class SLICSegment {
public:
  SLICSegment(
      int req_num_of_super_pixels = 1000 /*Desired number of superpixels*/,
      double compactness = 20 /*Compactness factor. use a value ranging from 10 to 40 */);

  enum ResultType {
    SUPER_PIXEL_BOUNDARY,
    LABELS
  };

  /** performSegmentation over input image
   * @param input image
   * @param result image
   * @param type ResultType (Default SUPER_PIXEL_BOUNDARY)
   * @return num of labels
   */
  int performSegmentation(const cv::Mat & input, cv::Mat & result, const ResultType type =
      SUPER_PIXEL_BOUNDARY) const;

private:
  int req_num_of_super_pixels_;
  double compactness_;
};

#endif // SLIC_SEGMENT_H_
