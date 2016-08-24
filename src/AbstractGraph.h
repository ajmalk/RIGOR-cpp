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
#include "AbstractImage.h"


struct un_edge {
    edgew   lambda_s,
            lambda_t,
            nonlambda_s,
            nonlambda_t;
};


class AbstractGraph {
protected:
    AbstractGraph(Seed &seed, AbstractImage *image, int pw_type, std::string type, std::vector< float > &lambdas) :
        seed(seed), image(image), pw_type(pw_type), type(type), lambdas(lambdas) { }
    
public:
    
    class UnaryIterator {
    public:
        typedef UnaryIterator self_type;
        typedef std::pair<edgew, edgew> value_type;
        typedef value_type& reference;
        typedef value_type* pointer;
        typedef std::forward_iterator_tag iterator_category;
        typedef int difference_type;
        
        self_type operator++() { index++; return *this; }
        self_type operator++(int junk) { index++; return *this; }
        value_type operator*() { return graph.get_unary(index, lambda, lambda_bar); }
        bool operator==(const self_type& rhs) { return index == rhs.index; }
        bool operator!=(const self_type& rhs) { return index != rhs.index; }
    private:
        friend class AbstractGraph;
        AbstractGraph& graph;
        int index = 0;
        double lambda, lambda_bar;
        UnaryIterator(AbstractGraph& graph, double lambda, double lambda_bar)
                : graph(graph), lambda(lambda), lambda_bar(lambda_bar) { }
    };

    class PairwiseIterator {
    public:
        typedef PairwiseIterator self_type;
        typedef pw_edge value_type;
        typedef value_type& reference;
        typedef value_type* pointer;
        typedef std::forward_iterator_tag iterator_category;
        typedef int difference_type;

        self_type operator++() { index++; return *this; }
        self_type operator++(int junk) { index++; return *this; }
        value_type operator*() { return graph.image->get_pairwise(graph.pw_type, index); }
        bool operator==(const self_type& rhs) { return index == rhs.index; }
        bool operator!=(const self_type& rhs) { return index != rhs.index; }
    private:
        friend class AbstractGraph;
        AbstractGraph& graph;
        int index = 0;
        PairwiseIterator(AbstractGraph& graph): graph(graph) { }
    };

    AbstractImage *image;

    Seed &seed;
    std::string type;
    std::vector< float > &lambdas;
    int pw_type;

    virtual inline edgew get_lambda_s(int pixel) = 0;
    virtual inline edgew get_lambda_t(int pixel) = 0;
    virtual inline edgew get_nonlambda_s(int pixel) = 0;
    virtual inline edgew get_nonlambda_t(int pixel) = 0;

    inline edgew seed_inv_distance(int p, dist_fun distance) {
        edgew minim = infinity;
        for ( auto s: seed ) {
            minim = std::min(minim, distance(image->get_color(p), image->get_color(s)));
        }
        return exp(-minim / 2.0);
    }

    inline std::pair<edgew, edgew> get_unary(int sp, double lambda, double lambda_bar) {
//        auto dist = eucl_distance(image->get_color(sp), image->get_color(*seed.begin()));
//        std::cout << dist << " " << exp(-dist/256.0) << std::endl;
//        std::cout << get_lambda_s(sp) << " " << get_nonlambda_s(sp) << " " << get_lambda_t(sp) << " " << get_nonlambda_t(sp) << std::endl;
        return std::make_pair(get_lambda_s(sp) * lambda + get_nonlambda_s(sp),
                              get_lambda_t(sp) * lambda_bar + get_nonlambda_t(sp));
    }
    
    static inline edgew eucl_distance(color_t color1, color_t color2) {
        color_t w = color1 - color2;
        return sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
    }
    
    inline UnaryIterator get_unaries(double lambda, double lambda_bar){
        return UnaryIterator(*this, lambda, lambda_bar);
    }

    inline PairwiseIterator get_pairwise_it() {
        return PairwiseIterator(*this);
    }

    inline int getNumPairwise() {
        return image->getNumPairwise();
    }

};


//class SPGraph: public AbstractGraph {
//    pws &spixels;
//    PWEdges *pairwise;
//
//protected:
//    SPGraph(Seed &seed, SPixels &spixels, PWEdges *pairwise, std::string type, std::vector< float > &lambdas) :
//        AbstractGraph(seed, "UniformInternGraph", lambdas(lambdas)), spixels(spixels), pairwise(pairwise)  { }
//};
//
//
//class PixelGraph: public AbstractGraph {
//    PixelGraph(Seed &seed, SPixels &pixels, std::string type, std::vector< float > &lambdas) :
//        AbstractGraph(seed, "UniformInternGraph", lambdas(lambdas)), spixels(spixels), pairwise(pairwise)  { }
//};


typedef std::unique_ptr<AbstractGraph> (*create_graph_func)(Seed &, AbstractImage *image);

#endif