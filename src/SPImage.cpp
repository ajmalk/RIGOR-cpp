#include "SPImage.h"

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/contrib/contrib.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <iostream>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iterator>
#include <climits>
#include <numeric>

#include "SLICSegment.h"

namespace fs = boost::filesystem;
namespace po = boost::program_options;

using namespace std;

int SPImage::slic_segmentation(cv::Mat &image, cv::Mat &sp_im, int num_sp, double compactness){
    
    SLICSegment::ResultType out_sp_type = SLICSegment::LABELS;
    SLICSegment slic(num_sp, compactness);
    num_sp = slic.performSegmentation(image, sp_im, out_sp_type);
    
    return num_sp;
}

template <typename T>
void print_vector (vector<T> vec) {
    cout << "[ ";
    std::copy(vec.begin(), vec.end(),
              std::ostream_iterator<T>(std::cout, " "));
    cout << "]" << endl;
}

void print_set(set<int> vec){
    cout << "{ ";
    std::copy(vec.begin(), vec.end(),
              std::ostream_iterator<int>(std::cout, " "));
    cout << "}" << endl;
}

void SPImage::build_pairwise_edges(){
    auto edge_hash = [this](const pair<int, int> &key) { return key.first * this->spixels.size() + key.second; };
    std::unordered_map< pair<int, int>, vector<edgew>, decltype(edge_hash) > b_pixels(0, edge_hash);
    b_pixels.reserve(spixels.size() * 10);
    
    auto sp_it = sp_im->begin<int>();
    auto edges_it = edges->begin<edgew>();
    for(int y = 0; y < sy - 1; y++, sp_it++, edges_it++) {
        for (int x = 0; x < sx - 1; x++, sp_it++, edges_it++) {
            int sp = *sp_it, right = *(sp_it + 1), below = *(sp_it + sx);
            edgew cur_e = *edges_it, right_e = *(edges_it + 1), below_e = *(edges_it + sx);
            if (sp != right) {
                auto r = b_pixels.emplace( make_pair(sp, right), vector<edgew>() );
                auto l = b_pixels.emplace( make_pair(right, sp), vector<edgew>() );
                auto w = (cur_e + right_e) / 2;
                r.first->second.emplace_back(w);
                l.first->second.emplace_back(w);
            }
            if (sp != below) {
                auto r = b_pixels.emplace( make_pair(sp, below), vector<edgew>() );
                auto l = b_pixels.emplace( make_pair(below, sp), vector<edgew>() );
                auto w = (cur_e + below_e) / 2;
                r.first->second.emplace_back(w);
                l.first->second.emplace_back(w);
            }
        }
        // last column
        int sp = *sp_it, below = *(sp_it + sx);
        double cur_e = *edges_it, below_e = *(edges_it + sx);
        if (sp != below) {
            auto r = b_pixels.emplace( make_pair(sp, below), vector<edgew>() );
            auto l = b_pixels.emplace( make_pair(below, sp), vector<edgew>() );
            auto w = (cur_e + below_e) / 2;
            r.first->second.emplace_back(w);
            l.first->second.emplace_back(w);
        }
    }
    // last row
    for (int x = 0; x < sx - 1; x++, sp_it++, edges_it++) {
        int sp = *sp_it, right = *(sp_it + 1);
        edgew cur_e = *edges_it, right_e = *(edges_it + 1);
        if (sp != right) {
            auto r = b_pixels.emplace( make_pair(sp, right), vector<edgew>() );
            auto l = b_pixels.emplace( make_pair(right, sp), vector<edgew>() );
            auto w = (cur_e + right_e) / 2;
            
            r.first->second.emplace_back(w);
            l.first->second.emplace_back(w);
        }
    }
    
    for ( int i = 0; i < pw_params.size(); i++ ){
        this->pairwise.emplace_back();
        this->pairwise[i].reserve(b_pixels.size());
    }
    
    for(auto const & edge: b_pixels){
        for (int i = 0; i < pw_params.size(); i++) {
            PWEdges &edges = this->pairwise[i];
            auto &params = pw_params[i];
            edgew avg = accumulate(edge.second.begin(), edge.second.end(), 0.0) / edge.second.size();
            edgew inv = exp( - double(avg) / 2.0 / params.sig_s / params.sig_s );
            edgew ev = (params.Pc * inv + params.Pw + 0.007) * params.Po;
//            ev2 = 100;
            edges.emplace_back( pw_edge {edge.first.first, edge.first.second, ev} );
//            cout << edge.first.first << " " << edge.first.second << " " << avg << endl;
//            print_vector<double>(edge.second);
        }
    }
    
//    for (int i = 0; i < spixels.size(); i++) {
//        cout<< pairwise[i].a << " " << pairwise[i].b << ": " << pairwise[i].w << endl;
//    }
    
}

/*
 Takes in a superpixel map and instatiates an SPImage by generating an array of spixels
 */
SPImage::SPImage(cv::Mat &image, cv::Mat &sp_im, cv::Mat &edges, int num_sp){
    
    this->image = &image;
    this->num_sp = num_sp;
    this->sp_im = &sp_im;
    this->edges = &edges;
    this->pw_params = pw_params;
    
    sy = image.size[0];
    sx = image.size[1];
    
    im_data = (uchar *) image.data;
    sp_data = (int *) sp_im.data;
    spixel z = {0, 0};
    spixels = vector<spixel>(num_sp, z);
    
    auto im_it = image.begin<cv::Vec3b>(), im_it_end = image.end<cv::Vec3b>();
    auto sp_it = sp_im.begin<int>();
    for(int i = 0; im_it != im_it_end; ++im_it, ++sp_it, i++) {
        spixel *p = &spixels[*sp_it];
        p->size++;
        p->color += (*im_it);
        p->ext |= i < sx || i > sx * (sy - 1) || i % sx == 0 || i % sx == sx - 1;
    }
    
    for (spixel &p: spixels) {
        p.color /= p.size;
//        cout<< p.color << endl;
    }
    
    build_pairwise_edges();
}

cv::Mat SPImage::get_color_sp_im(){
    cv::Mat color_sp = cv::Mat(sy, sx, CV_8UC3);
    auto im_it = color_sp.begin< cv::Vec<uchar, 3> >();
    auto sp_it = sp_im->begin<int>(), sp_it_end = sp_im->end<int>();
    for(; sp_it != sp_it_end; ++sp_it) {
        spixel *p = &spixels[*sp_it];
        *(im_it++) = p->color;
    }
    return color_sp;
}

cv::Mat SPImage::cut_to_image(GraphCut& cut){
    cv::Mat cut_im = cv::Mat(sy, sx, CV_8U);
    auto im_it = cut_im.begin<uchar>(), im_it_end = cut_im.end<uchar>();
    auto sp_it = sp_im->begin<int>();
    for(; im_it != im_it_end; ++im_it, ++sp_it) {
        *im_it = cut.in_source_seg(*sp_it) * 255;
    }
    return cut_im;
}

cv::Mat SPImage::cut_to_image(vector<bool>::iterator cut){
    
    cv::Mat cut_im = cv::Mat(sy, sx, CV_8U);
    auto im_it = cut_im.begin<uchar>();
    auto sp_it = sp_im->begin<int>();
    for(; im_it != cut_im.end<uchar>(); ++im_it, ++sp_it) {
        *im_it = *(cut + *sp_it) * 255;
    }
    return cut_im;
}

cv::Mat SPImage::seeds_to_sp_im(std::vector< std::set<int> > &seeds){
    cv::Mat img = cv::Mat(sy, sx, CV_8U);
    auto im_it = img.begin<uchar>(), im_it_end = img.end<uchar>();
    auto sp_it = sp_im->begin<int>();
    for(; im_it != im_it_end; ++im_it, ++sp_it) {
        *im_it = 0;
        for (auto &seed: seeds){
            *im_it |= seed.count(*sp_it) * 255;
        }
    }
    return img;
}

vector< set<int> > SPImage::generate_seeds(cv::Mat sp_img, int num_seeds, int radius) {
    
    const int   sy = sp_img.size[0],
                sx = sp_img.size[1],
                rows = ceil(sqrt(sx / sy * num_seeds));
    int * data = (int *) sp_img.data;
    
    vector<int> row_cols = vector<int>(rows);
    for (int i = 0, rem_seeds = num_seeds, rem_rows = rows; i < float(rows) / 2; i++) {
        row_cols[i] = rem_seeds / rem_rows--;
        rem_seeds -= row_cols[i];
        if (rem_rows) {
            row_cols[rows - i - 1] = rem_seeds / rem_rows--;
            rem_seeds -= row_cols[rows - i - 1];
        }
    }
    
    vector< set<int> > seeds = vector< set<int> >(num_seeds);
    for (int row = 0, seed = 0; row < rows; row++) {
        int y = sy / (rows + 1) * (row + 1);
        for (int col = 0; col < row_cols[row]; col++, seed++) {
            int x = sx / (row_cols[row] + 1) * (col + 1);
            for (int i = 0; i < radius; i++) {
                for (int j = 0; j < radius; j++) {
                    seeds[seed].insert( data[x + i + (y + j) * sx] );
                    seeds[seed].insert( data[x - i + (y - j) * sx] );
                    seeds[seed].insert( data[x + i + (y - j) * sx] );
                    seeds[seed].insert( data[x - i + (y + j) * sx] );
                }
            }
        }
    }
    return seeds;
}

