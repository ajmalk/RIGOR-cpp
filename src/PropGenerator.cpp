//
//  PropGenerator.cpp
//  Rigor
//
//  Created by Ajmal Kunnummal on 2/20/15.
//  Copyright (c) 2015 ajmal. All rights reserved.
//

#include <opencv2/core/core.hpp>
#include <opencv2/contrib/contrib.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "bk_dynamicgraphs.h"
#include "SLICSegment.h"

#include "PropGenerator.h"
#include "Parameters.h"
#include "SPImage.h"
#include "PixelImage.h"
#include "FuxinGraph.h"
#include "SimpleBK.h"

#include <iostream>
#include <memory>
#include <vector>


using namespace std;

cv::Mat image;
cv::Mat sp_im;

PropGenerator::PropGenerator(Parameters &params){
    this->params = params;
}

void print_pos_val(int event, int x, int y, int flags, void* image){
    
    if  ( event == cv::EVENT_LBUTTONDOWN ) {
        cv::Mat * imp = (cv::Mat *) image;
        cv::Mat im = *imp;
        cout << "click: " << ((cv::Mat *) image)->at<int>(y, x)  << endl;
    }
}

void show_in_color(string window, cv::Mat &image){
    cv::namedWindow(window, cv::WINDOW_AUTOSIZE);
    double min, max;
    cv::minMaxIdx(image, &min, &max);
    cv::Mat adj_sp_map;
    // expand your range to 0..255. Similar to histEq();
    image.convertTo(adj_sp_map, CV_8UC1, 255 / (max-min), -min);
    // convert to a colormaped image
    cv::applyColorMap(adj_sp_map, adj_sp_map, cv::COLORMAP_JET);
    cv::imshow(window, adj_sp_map);
    cv::setMouseCallback(window, print_pos_val, &image);
}

inline cv::Mat iterator_to_image(AbstractGraph::UnaryIterator it, int sx, int sy) {
    cv::Mat cut_im = cv::Mat(sy, sx, CV_8U);
    auto im_it = cut_im.begin<uchar>(), im_it_end = cut_im.end<uchar>();
    auto p = 0;
    for(; im_it != im_it_end; ++im_it, ++it) {
        auto ed = (uchar) (*it).first;
        *im_it = ed;
        if(0 == (int) (*it).second) {
            cout << (int) (uchar) (*it).first << endl;
        }
    }
    return cut_im;
}

cv::Mat gradient(cv::Mat image){
    
    int scale = 1;
    int delta = 0;
    int ddepth = CV_16S;
    
    GaussianBlur( image, image, cv::Size(3,3), 0, 0, cv::BORDER_DEFAULT );
    
    /// Convert it to gray
    cv::Mat image_gray;
    cvtColor( image, image_gray, CV_RGB2GRAY );
    
    /// Generate grad_x and grad_y
    cv::Mat grad_x, grad_y;
    cv::Mat abs_grad_x, abs_grad_y;
    
    /// Gradient X
    //Scharr( src_gray, grad_x, ddepth, 1, 0, scale, delta, BORDER_DEFAULT );
    Sobel( image_gray, grad_x, ddepth, 1, 0, 3, scale, delta, cv::BORDER_DEFAULT );
    convertScaleAbs( grad_x, abs_grad_x );
    
    /// Gradient Y
    //Scharr( src_gray, grad_y, ddepth, 0, 1, scale, delta, BORDER_DEFAULT );
    Sobel( image_gray, grad_y, ddepth, 0, 1, 3, scale, delta, cv::BORDER_DEFAULT );
    convertScaleAbs( grad_y, abs_grad_y );
    
    /// Total Gradient (approximate)
    cv::Mat grad;
    addWeighted( abs_grad_x, 0.5, abs_grad_y, 0.5, 0, grad );
    
    cv::Mat gradf = cv::Mat(image.size[0], image.size[1], CV_64F);
    auto i_it = grad.begin<uchar>();
    for(auto f_it = gradf.begin<double>() ; f_it != gradf.end<double>(); ++i_it, ++f_it) {
        *f_it = float(*i_it);
    }
    return gradf;
}

cv::Mat struct_edges(cv::Mat image){
    
    /// Convert it to gray
    cv::Mat image_gray;
    cvtColor( image, image_gray, CV_RGB2GRAY );
    
    cv::Mat grayf = cv::Mat(image.size[0], image.size[1], CV_64F);
    auto i_it = image_gray.begin<uchar>();
    for(auto f_it = grayf.begin<double>() ; f_it != grayf.end<double>(); ++i_it, ++f_it) {
        *f_it = float(*i_it);
    }
    return grayf;
}

void PropGenerator::generate(fs::path filename){
    image = cv::imread(filename.c_str());
    
    if (image.empty()) {
        // Check for invalid input
        std::cout <<  "Could not open or find the image" << std::endl;
        return;
    }
    
    if(params.debug){
        cout << "Read image" << endl;
        cv::namedWindow("Input Image", cv::WINDOW_AUTOSIZE);
        cv::imshow("Input Image", image);
    }
    
    cv::Mat edges = gradient(image);
    //cv::Mat edges = cv::imread("/Users/ajmalkunnummal/Pictures/peppers_str_edges_fat.png");
    //edges = struct_edges(edges);
    
    if(params.debug){
        cout << "Ran edge detection" << endl;
        cv::namedWindow( "Edge Detector", CV_WINDOW_AUTOSIZE );
        imshow( "Edge Detector", edges );
    }
    
    int num_sp = SPImage::slic_segmentation(image, sp_im, params.num_sp, params.compactness);
    
    
    if(params.debug) {
        cout << "Segmented image. Num of superpixels: " << num_sp << endl;
        show_in_color("Superpixels", sp_im);
    }
    
    SPImage spixels = SPImage(image, sp_im, edges, num_sp);
    
    if(params.debug) {
        cout << "Setup graph and pairwise edges" << endl;
        cv::namedWindow("Superpixel Colors", cv::WINDOW_AUTOSIZE);
        cv::imshow("Superpixel Colors", spixels.get_color_sp_im());
    }
    
    vector< set<int> > seeds = SPImage::generate_seeds(sp_im, params.seeds, params.seed_radius = 5);

    if(params.debug){
        cout << "Generated " << seeds.size() << " seeds" << endl;
        cv::imshow("Seeds", spixels.seeds_to_sp_im(seeds));
    }
    
    vector< unique_ptr<AbstractGraph> > graphs;
    
    for ( auto& seed: seeds) {
        for (auto graph_type: params.graph_types) {
            graphs.push_back( graph_type(seed, &spixels) );
        }
    }

    AbstractGraph::UnaryIterator unaries = graphs[0]->get_unaries(graphs[0]->lambdas[0],
                                                                  graphs[0]->lambdas[graphs[0]->lambdas.size() - 1]);

    vector<edgew> stats;
//    for (int i = 0; i < graphs[0]->image->getNumPixels(); unaries++, i++) {
//        std::pair<edgew, edgew> unary_pair = *unaries;
////        cout << unary_pair.first << " " << unary_pair.second << endl;
//        stats.push_back(unary_pair.first);
//    }

    for(int i = 1; i < spixels.getNumPairwise(); i+=2) {
        auto pw = spixels.get_pairwise(0, i);
        stats.push_back(pw.w);
    }

    sort(stats.begin(), stats.end());

    cout << spixels.getNumPairwise() << endl;
    cout << stats[0] << endl;
    cout << stats[stats.size() / 4] << endl;
    cout << stats[stats.size() / 2] << endl;
    cout << stats[stats.size() / 4 * 3] << endl;
    cout << stats[stats.size() / 8 * 7] << endl;
    cout << stats[stats.size() / 16 * 15] << endl;
    cout << stats[stats.size() - 7] << endl;
    cout << stats[stats.size() - 1] << endl;
    
    if(params.debug){
        cout << "Setup " << graphs.size() << " graphs (1 for each seed and graph type)" << endl;
    }
    
    auto cuts = nodynamic_param_maxflow(graphs);
    
    if(params.debug){
        cout << "Found " << cuts.size() / spixels.spixels.size() << " cuts" << endl;
//        for(int i = 0; i < cuts.size() / spixels.spixels.size(); i++){
//            for(int j = 0; j < spixels.spixels.size(); j++){
//                cout << cuts[ i * spixels.spixels.size() + j];
//            }
//            cout << endl;
//        }
        for (auto cur_cut = cuts.begin(); cur_cut != cuts.end(); cur_cut += spixels.spixels.size() ) {
            cv::imshow("Cut", spixels.cut_to_image(cur_cut));
            cv::waitKey(0);
        }
    }
    
    
    if(params.debug)
        cv::waitKey(0);
}

void PropGenerator::generate_pix(fs::path filename){
    image = cv::imread(filename.c_str());

    if (image.empty()) {
        // Check for invalid input
        std::cout <<  "Could not open or find the image" << std::endl ;
        return;
    }

    if(params.debug){
        cout << "Read image" << endl;
        cv::namedWindow("Input Image", cv::WINDOW_AUTOSIZE);
        cv::imshow("Input Image", image);
    }

    cv::Mat edges = gradient(image) ;
    //cv::Mat edges = cv::imread("/Users/ajmalkunnummal/Pictures/peppers_str_edges_fat.png");
    //edges = struct_edges(edges);

    if(params.debug){
        cout << "Ran edge detection" << endl;
        cv::namedWindow( "Edge Detector", CV_WINDOW_AUTOSIZE );
        imshow( "Edge Detector", edges);
    }

    PixelImage spixels = PixelImage(image, edges);

    cout << "Num of pixels: " << spixels.getNumPixels() << endl;

    cout << "unary color: " << spixels.get_color(0) << endl;
    cout << "unary size:" << spixels.get_size(0) << endl;
    cout << "unary ext:" << spixels.get_ext(0) << endl;

    vector< set<int> > seeds = { { spixels.getNumPixels() / 2 + spixels.sx / 2 } };

    if(params.debug) {
        cout << "Generated " << seeds.size() << " seeds" << endl;
//        cv::imshow("Seeds", spixels.seeds_to_sp_im(seeds));
    }

    vector< unique_ptr<AbstractGraph> > graphs;

    for ( auto& seed: seeds) {
        for (auto graph_type: params.graph_types) {
            graphs.push_back( graph_type(seed, &spixels) );
        }
    }




    auto unaries2 = graphs[0]->get_unaries(graphs[0]->lambdas[0],
                                     graphs[0]->lambdas[graphs[0]->lambdas.size() - 0 - 1]);

    imshow( "Unaries", iterator_to_image(unaries2, spixels.sx, spixels.sy) );


    AbstractGraph::UnaryIterator unaries = graphs[0]->get_unaries(graphs[0]->lambdas[0],
                                                                  graphs[0]->lambdas[graphs[0]->lambdas.size() - 0 - 1]);
    vector<edgew> stats;

//    for (int i = 0; i < graphs[0]->image->getNumPixels(); unaries++, i++) {
//        std::pair<edgew, edgew> unary_pair = *unaries;
////        cout << unary_pair.first << " " << unary_pair.second << endl;
//        stats.push_back(unary_pair.first);
//    }


    for(int i = 1; i < spixels.getNumPairwise(); i+=2) {
        auto pw = spixels.get_pairwise(0, i);
        stats.push_back(-pw.w);
    }

    sort(stats.begin(), stats.end());

    cout << stats[0] << endl;
    cout << stats[stats.size() / 4] << endl;
    cout << stats[stats.size() / 2] << endl;
    cout << stats[stats.size() / 4 * 3] << endl;
    cout << stats[stats.size() / 8 * 7] << endl;
    cout << stats[stats.size() / 16 * 15] << endl;
    cout << stats[stats.size() - 7] << endl;
    cout << stats[stats.size() - 1] << endl;
    cout << "highest" << endl;

    if(params.debug){
        cout << "Setup " << graphs.size() << " graphs (1 for each seed and graph type)" << endl;
    }

    auto cuts = nodynamic_param_maxflow(graphs);

    if(params.debug){
        cout << "Found " << cuts.size() / spixels.getNumPixels() << " cuts" << endl;
        for (auto cur_cut = cuts.begin(); cur_cut != cuts.end(); cur_cut += spixels.getNumPixels() ) {
            cv::imshow("Cut", spixels.cut_to_image(cur_cut));
            cv::waitKey(0);
        }
    }


    if(params.debug)
        cv::waitKey(0);
}