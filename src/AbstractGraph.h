//
//  AbstractGraph.h
//  rigor
//
//  Created by Ajmal Kunnummal on 4/11/15.
//
//

#ifndef __rigor__AbstractGraph__
#define __rigor__AbstractGraph__

#include <vector>
#include <set>
#include <climits>
#include <opencv2/core/core.hpp>
#include <iostream>
#include <cmath>

typedef double edgew;
typedef cv::Vec3f color_t;

const edgew infinity = 21475000000;

struct spixel{
    int size;
    color_t color;
    bool ext;
};

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

typedef std::vector<spixel> SPixels;
typedef std::vector<pw_edge> PWEdges;
typedef std::set<int> Seed;

struct un_edge {
    edgew   lambda_s,
            lambda_t,
            nonlambda_s,
            nonlambda_t;
};

class GraphCut {
public:
    virtual inline void run() = 0;
    virtual inline bool in_source_seg(int sp) = 0;
};

class AbstractGraph {
protected:
    AbstractGraph(Seed &seed, SPixels &spixels, PWEdges *pairwise, std::string type, std::vector< float > &lambdas) :
        spixels(spixels), pairwise(pairwise), seed(seed), type(type), lambdas(lambdas) { }
    
public:
    
    class UnaryIterator {
    public:
        typedef UnaryIterator self_type;
        typedef std::pair<edgew, edgew> value_type;
        typedef value_type& reference;
        typedef value_type* pointer;
        typedef std::forward_iterator_tag iterator_category;
        typedef int difference_type;
        
        self_type operator++() { sp++; return *this; }
        self_type operator++(int junk) { sp++; return *this; }
        value_type operator*() { return graph.get_unary(sp, lambda, lambda_bar); }
        bool operator==(const self_type& rhs) { return sp == rhs.sp; }
        bool operator!=(const self_type& rhs) { return sp != rhs.sp; }
    private:
        friend class AbstractGraph;
        AbstractGraph& graph;
        int sp;
        double lambda, lambda_bar;
        UnaryIterator(AbstractGraph& graph, int sp, double lambda, double lambda_bar) : graph(graph), sp(sp), lambda(lambda), lambda_bar(lambda_bar) { }
    };
    
    SPixels &spixels;
    PWEdges *pairwise;
    Seed &seed;
    std::string type;
    std::vector< float > &lambdas;
    
    virtual inline edgew get_lambda_s(int sp) = 0;
    virtual inline edgew get_lambda_t(int sp) = 0;
    virtual inline edgew get_nonlambda_s(int sp) = 0;
    virtual inline edgew get_nonlambda_t(int sp) = 0;
    
    inline std::pair<edgew, edgew> get_unary(int sp, double lambda, double lambda_bar) {
//        std::cout << get_lambda_s(sp) + get_nonlambda_s(sp) << " " << get_lambda_t(sp) + get_nonlambda_t(sp) << std::endl;
        return std::make_pair(get_lambda_s(sp) * lambda + get_nonlambda_s(sp),
                              get_lambda_t(sp) * lambda_bar + get_nonlambda_t(sp));
    }
    
    inline edgew get_pairwise(int sp1, int sp2);
    
    static inline edgew eucl_distance(color_t color1, color_t color2) {
        color_t w = color1 - color2;
        return sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
    }
    
    inline edgew seed_inv_distance(int sp, dist_fun distance){
        edgew minim = infinity;
        for ( auto s: seed ) {
            minim = std::min(minim, distance(spixels[sp].color, spixels[s].color));
        }
        return exp(-minim / 2.0);
    }
    
    inline UnaryIterator get_unaries(double lambda, double lambda_bar){
        return UnaryIterator(*this, 0, lambda, lambda_bar);
    }
    
};

typedef std::unique_ptr<AbstractGraph> (*create_graph_func)(Seed &, SPixels &, std::vector<PWEdges> &);

#endif