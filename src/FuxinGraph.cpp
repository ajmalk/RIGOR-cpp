//
//  FuxinGraph.cpp
//  rigor
//
//  Created by Ajmal Kunnummal on 4/24/15.
//
//

#include "FuxinGraph.h"
#include <iostream>


namespace FuxinGraphs {
    std::unique_ptr<AbstractGraph> uniform_intern_graph(Seed &seed, SPixels &spixels, std::vector<PWEdges> &pairwise){
        return std::make_unique<UniformInternGraph>(seed, spixels, pairwise);
    }
    std::unique_ptr<AbstractGraph> uniform_extern_graph(Seed &seed, SPixels &spixels, std::vector<PWEdges> &pairwise){
        return std::make_unique<UniformExternGraph>(seed, spixels, pairwise);
    }
    std::unique_ptr<AbstractGraph> uniform_extern2_graph(Seed &seed, SPixels &spixels, std::vector<PWEdges> &pairwise){
        return std::make_unique<UniformExtern2Graph>(seed, spixels, pairwise);
    }
    std::unique_ptr<AbstractGraph> color_intern_graph(Seed &seed, SPixels &spixels, std::vector<PWEdges> &pairwise){
        return std::make_unique<ColorInternGraph>(seed, spixels, pairwise);
    }
    std::unique_ptr<AbstractGraph> color_extern_graph(Seed &seed, SPixels &spixels, std::vector<PWEdges> &pairwise){
        return std::make_unique<ColorExternGraph>(seed, spixels, pairwise);
    }
    std::unique_ptr<AbstractGraph> color_extern2_graph(Seed &seed, SPixels &spixels, std::vector<PWEdges> &pairwise){
        return std::make_unique<ColorExtern2Graph>(seed, spixels, pairwise);
    }
};
