//
//  FuxinGraph.h
//  rigor
//
//  Created by Ajmal Kunnummal on 4/12/15.
//
//

#ifndef rigor_FuxinGraph_h
#define rigor_FuxinGraph_h

#include "AbstractGraph.h"
#include <memory>
#include <array>

namespace FuxinGraphs {

    enum PWType { UNIFORM, COLOR };

    static const std::vector<PWParams> PWPARAMS = { { 1, 3.5, 1e-3, 1000}, { 1, 3.5, 1e-3, 1000} };
    
    static const std::array< float, 2 > lambda_bounds = { { 20., 150. } };
    static const int NLAMBDAS = 20;
    static const float SMALLEST_LAMBDA = 0.001;
    
    template<typename T>
    class Logspace {
    private:
        T curValue, base;
        
    public:
        Logspace(T base) : curValue(SMALLEST_LAMBDA), base(base) {}
        
        T operator()() {
            return curValue *= base;
        }
    };
    
    template <typename T>
    void print_vector (std::vector<T> vec) {
        std::cout << "[ ";
        std::copy(vec.begin(), vec.end(),
                  std::ostream_iterator<T>(std::cout, " "));
        std::cout << "]" << std::endl;
    }
    
    inline auto getLamdas(){
        std::vector< std::vector< float > > lambdas(2);
        lambdas[UNIFORM].reserve(NLAMBDAS);
        lambdas[COLOR].reserve(NLAMBDAS);
        lambdas[UNIFORM].push_back(0);
        lambdas[COLOR].push_back(0);
        
        float baseu = pow(10, (log10(lambda_bounds[UNIFORM]) - log10(SMALLEST_LAMBDA)) / NLAMBDAS);
        float basec = pow(10, (log10(lambda_bounds[COLOR]  ) - log10(SMALLEST_LAMBDA)) / NLAMBDAS);
        
        std::generate_n(std::back_inserter(lambdas[UNIFORM]), NLAMBDAS, Logspace<float>(baseu));
        std::generate_n(std::back_inserter(lambdas[COLOR]), NLAMBDAS, Logspace<float>(basec));
        
//        print_vector(lambdas[COLOR ]);
        
        return lambdas;
    }
    
    static std::vector< std::vector< float > > fuxin_lambdas = getLamdas();

    class UniformInternGraph : public AbstractGraph {
    public:
        UniformInternGraph(Seed &seed, SPixels &spixels, std::vector<PWEdges> &pairwise) :
        AbstractGraph(seed, spixels, &pairwise[UNIFORM], "UniformInternGraph", fuxin_lambdas[UNIFORM]) {}
        
        inline edgew get_lambda_s(int sp) {
            return !seed.count(sp) * spixels[sp].size;
        }
        inline edgew get_lambda_t(int sp) {
            return !seed.count(sp) * spixels[sp].size / 2;
        }
        inline edgew get_nonlambda_s(int sp) {
            return seed.count(sp) * infinity;
        }
        inline edgew get_nonlambda_t(int sp) {
            return 0;
        }
    };

    class UniformExternGraph : public AbstractGraph {
    public:
        UniformExternGraph(Seed &seed, SPixels &spixels, std::vector<PWEdges> &pairwise) :
            AbstractGraph(seed, spixels, &pairwise[UNIFORM], "UniformExternGraph", fuxin_lambdas[UNIFORM]) {  }
        
        inline edgew get_lambda_s(int sp) {
            return !seed.count(sp) * spixels[sp].ext * spixels[sp].size;
        }
        inline edgew get_lambda_t(int sp) {
            return !seed.count(sp) * spixels[sp].size / 256;
        }
        inline edgew get_nonlambda_s(int sp) {
            return seed.count(sp) * infinity;
        }
        inline edgew get_nonlambda_t(int sp) {
            return spixels[sp].ext;
        }
    };

    class UniformExtern2Graph : public AbstractGraph {
    public:
        UniformExtern2Graph(Seed &seed, SPixels &spixels, std::vector<PWEdges> &pairwise) :
            AbstractGraph(seed, spixels, &pairwise[UNIFORM], "UniformExtern2Graph", fuxin_lambdas[UNIFORM]) {  }
        
        inline edgew get_lambda_s(int sp) {
            return !seed.count(sp) * spixels[sp].size;
        }
        inline edgew get_lambda_t(int sp) {
            return !seed.count(sp) * spixels[sp].size / 256;
        }
        inline edgew get_nonlambda_s(int sp) {
            return seed.count(sp) * infinity;
        }
        inline edgew get_nonlambda_t(int sp) {
            return 0;
        }
    };

    class ColorInternGraph : public AbstractGraph {
    private:
        float scaling = 220;
    public:
        ColorInternGraph(Seed &seed, SPixels &spixels, std::vector<PWEdges> &pairwise) :
            AbstractGraph(seed, spixels, &pairwise[COLOR], "ColorInternGraph", fuxin_lambdas[COLOR]) {  }
        
        inline edgew get_lambda_s(int sp) {
            return 0;
        }
        inline edgew get_lambda_t(int sp) {
            return !seed.count(sp);
        }
        inline edgew get_nonlambda_s(int sp) {
            return seed.count(sp) * infinity + !seed.count(sp) * this->seed_inv_distance(sp, AbstractGraph::eucl_distance) * scaling;
        }
        inline edgew get_nonlambda_t(int sp) {
            return 0;
        }
    };

    class ColorExternGraph : public AbstractGraph {
    public:
        ColorExternGraph(Seed &seed, SPixels &spixels, std::vector<PWEdges> &pairwise) :
            AbstractGraph(seed, spixels, &pairwise[COLOR], "ColorExternGraph", fuxin_lambdas[COLOR]) { }
        
        inline edgew get_lambda_s(int sp) {
            return !seed.count(sp) * spixels[sp].size;
        }
        inline edgew get_lambda_t(int sp) {
            return !seed.count(sp) * spixels[sp].size * !spixels[sp].ext;
        }
        inline edgew get_nonlambda_s(int sp) {
            return seed.count(sp) * infinity + !seed.count(sp) * this->seed_inv_distance(sp, AbstractGraph::eucl_distance);
        }
        inline edgew get_nonlambda_t(int sp) {
            return !seed.count(sp) * ( !spixels[sp].ext * this->seed_inv_distance(sp, AbstractGraph::eucl_distance) + spixels[sp].ext);
        }
    };

    class ColorExtern2Graph : public AbstractGraph {
    public:
        ColorExtern2Graph(Seed &seed, SPixels &spixels, std::vector<PWEdges> &pairwise) :
            AbstractGraph(seed, spixels, &pairwise[COLOR], "ColorExtern2Graph", fuxin_lambdas[COLOR] ) { }
        
        inline edgew get_lambda_s(int sp) {
            return !seed.count(sp) * spixels[sp].size;
        }
        inline edgew get_lambda_t(int sp) {
            return !seed.count(sp) * spixels[sp].size;
        }
        inline edgew get_nonlambda_s(int sp) {
            return seed.count(sp) * infinity + !seed.count(sp) * this->seed_inv_distance(sp, AbstractGraph::eucl_distance);
        }
        inline edgew get_nonlambda_t(int sp) {
            return !seed.count(sp) * this->seed_inv_distance(sp, AbstractGraph::eucl_distance);
        }
    };

    std::unique_ptr<AbstractGraph> uniform_intern_graph  ( Seed &seed, SPixels &spixels, std::vector<PWEdges> &pairwise );
    std::unique_ptr<AbstractGraph> uniform_extern_graph  ( Seed &seed, SPixels &spixels, std::vector<PWEdges> &pairwise );
    std::unique_ptr<AbstractGraph> uniform_extern2_graph ( Seed &seed, SPixels &spixels, std::vector<PWEdges> &pairwise );
    std::unique_ptr<AbstractGraph> color_intern_graph    ( Seed &seed, SPixels &spixels, std::vector<PWEdges> &pairwise );
    std::unique_ptr<AbstractGraph> color_extern_graph    ( Seed &seed, SPixels &spixels, std::vector<PWEdges> &pairwise );
    std::unique_ptr<AbstractGraph> color_extern2_graph   ( Seed &seed, SPixels &spixels, std::vector<PWEdges> &pairwise );
};


#endif
