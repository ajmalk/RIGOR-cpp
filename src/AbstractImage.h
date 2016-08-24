//
//  AbstractImage.h
//  Rigor
//
//  Created by Ajmal Kunnummal on 2/20/15.
//  Copyright (c) 2015 ajmal. All rights reserved.
//

#ifndef __Rigor__AbstractImage__
#define __Rigor__AbstractImage__

#include <opencv2/core/core.hpp>
#include <set>

typedef double edgew;

typedef cv::Vec3f color_t;

const edgew infinity = 21475000000;

typedef edgew (*dist_fun)(color_t, color_t);

struct pw_edge {
    int a, b;
    edgew w;
};

struct PWParams{
    float Pc;
    float sig_s;
    float Pw;
    float Po;
};

typedef std::vector<pw_edge> PWEdges;
typedef std::set<int> Seed;

struct spixel{
    int size;
    color_t color;
    bool ext;
};

typedef std::vector<spixel> SPixels;

class GraphCut {
public:
    virtual inline void run() = 0;
    virtual inline bool in_source_seg(int sp) = 0;
};

class AbstractImage {
protected:
    
    cv::Mat *image, *sp_im, *edges;
    int sx, sy;

    std::vector<PWParams> pw_params = { { 1, 3.5, 1e-3, 1000}, { 1, 3.5, 1e-3, 1000} };

    uchar *im_data;
    int *sp_data;

    void build_pairwise_edges();

public:

    virtual inline int getNumPixels() = 0;
    virtual inline int getNumPairwise() = 0;
    virtual inline pw_edge get_pairwise(int set, int index) = 0;
    virtual inline int get_size(int p) = 0;
    virtual inline color_t get_color(int p) = 0;
    virtual inline bool get_ext(int p) = 0;
};

#endif /* defined(__Rigor__AbstractImage__) */
