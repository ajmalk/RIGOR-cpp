//
//  PixelImage.h
//  Rigor
//
//  Created by Ajmal Kunnummal on 2/20/15.
//  Copyright (c) 2015 ajmal. All rights reserved.
//

#ifndef __Rigor__PixelImage__
#define __Rigor__PixelImage__

#include <opencv2/core/core.hpp>
#include "AbstractGraph.h"
#include <cmath>

using namespace cv;
using namespace std;

pair<int, int> neib[] = { {0, 1}, {1, 1}, {1, 0}, {-1, 1} };

class PixelImage : public AbstractImage {
private:

    cv::Mat *image, *edges;

    int size, num_pw;

    std::vector<PWParams> pw_params = { { 1, 3.5, 1e-3, 1000}, { 1, 3.5, 1e-3, 1000} };

    Vec3b *im_data;
    uchar *edge_data;

public:

    const int sx, sy;

    PixelImage(cv::Mat &image, cv::Mat &edges)
            : image(&image), edges(&edges),
              sy(image.size[0]), sx(image.size[1]) {
        size = sx * sy;
        num_pw = (sy - 1) * sx * 2 +
                 sy * (sx - 1) * 2 +
                 (sy - 1) * (sx - 1) * 4;
        im_data = (Vec3b *) image.data;
        edge_data = edges.data;

        cout << sx << ":" << sy << endl;
    }

    static std::vector< std::set<int> > generate_seeds(cv::Mat sp_img, int num_seeds, int radius = 1);
    cv::Mat seeds_to_sp_im(std::vector< std::set<int> > &seeds);
    cv::Mat cut_to_image(GraphCut& cut);

    inline cv::Mat cut_to_image(std::vector<bool>::iterator cut) {
        cv::Mat cut_im = cv::Mat(sy, sx, CV_8U);
        auto im_it = cut_im.begin<uchar>(), im_it_end = cut_im.end<uchar>();
        auto p = 0;
        for(; im_it != im_it_end; ++im_it, ++cut) {
            *im_it = *cut * 255;
        }
        return cut_im;
    }

    inline int getNumPixels() {
        return size;
    }

    inline int getNumPairwise() {
        return num_pw;
    }

    inline pw_edge get_pairwise(int set, int index) {
        int par = index % 2;
        index /= 2;
        int x1 = -1, y1 = -1, x2 = -1, y2 = -1;
        pair<int, int> n = {0, 0};
        bool diag = false;

        if( index < 4 * (sx - 2) * (sy - 1) ) {
            int p = index / 4;
            n = neib[index % 4];
            diag = index % 2;
            x1 = p % (sx - 2) + 1;
            y1 = p / (sx - 2);
            goto ret_pw;
        }

        index -= 4 * (sx - 2) * (sy - 1);

        if( index < 3 * (sy - 1) ) {
            int p = index / 3;
            n = neib[index % 3];
            diag = (index % 3) == 1;
            x1 = 0;
            y1 = p;
            goto ret_pw;
        }

        index -= 3 * (sy - 1);

        if( index < 2 * (sy - 1) ) {
            int p = index / 2;
            n = neib[index % 2];
            diag = index % 2;
            n.first = -n.first;
            x1 = sx - 1;
            y1 = p;
            goto ret_pw;
        }

        index -= 2 * (sy - 1);

        x1 = index;
        y1 = sy - 1;
        n = neib[2];

        ret_pw:

        x2 = x1 + n.first;
        y2 = y1 + n.second;

        int p1 = x1 + y1 * sx;
        int p2 = x2 + y2 * sx;
        edgew w = (*(edge_data + p1) + *(edge_data + p2)) / 2.0;
        if(diag) w = 0;

        auto &params = pw_params[set];
        edgew x = double(w) / 2.0 / params.sig_s / params.sig_s;
//        x = double(w);
        edgew inv = exp( - x );
        edgew ev = (params.Pc * inv + params.Pw + 0.007) * params.Po;

        if(par) {
//            if(x1 < 1 || y1 < 1)
//                cout << x1 << ":" << y1 << " " << x2 << ":" << y2 << " " << w << " " << ev << endl;
//            cout << params.Pc << " " << params.Pw << " " << params.Po << " " << params.sig_s << endl;
            return { p1, p2, ev };
        } else {
            return { p2, p1, ev };
        }
    }

    inline int get_size(int p) {
        return 1;
    }

    inline color_t get_color(int p) {
        return im_data[p];
    }

    inline bool get_ext(int p) {
        return !( (p % sx) && (p % sx - 1) && (p % sy) && (p % sy - 1) );
    }
};

#endif /* defined(__Rigor__PixelImage__) */
