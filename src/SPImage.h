//
//  SPImage.h
//  Rigor
//
//  Created by Ajmal Kunnummal on 2/20/15.
//  Copyright (c) 2015 ajmal. All rights reserved.
//

#ifndef __Rigor__SPImage__
#define __Rigor__SPImage__

#include <opencv2/core/core.hpp>
#include "AbstractGraph.h"

class SPImage {
private:
    
    cv::Mat *image, *sp_im, *edges;
    int sx, sy;
    int num_sp;
    
    std::vector<PWParams> pw_params = { { 1, 3.5, 1e-3, 1000}, { 1, 3.5, 1e-3, 1000} };

    uchar *im_data;
    int *sp_data;

    void build_pairwise_edges();

public:
    std::vector<spixel> spixels;
    std::vector< PWEdges > pairwise;

    SPImage(cv::Mat &image, cv::Mat &sp_im, cv::Mat &edges, int num_sp);

    static int slic_segmentation(cv::Mat &image, cv::Mat &sp_im, int num_sp, double compactness);
    static std::vector< std::set<int> > generate_seeds(cv::Mat sp_img, int num_seeds, int radius = 1);
    cv::Mat get_sp_image();
    cv::Mat get_color_sp_im();
    cv::Mat seeds_to_sp_im(std::vector< std::set<int> > &seeds);
    cv::Mat cut_to_image(GraphCut& cut);
    cv::Mat cut_to_image(std::vector<bool>::iterator cut);
};

#endif /* defined(__Rigor__SPImage__) */
