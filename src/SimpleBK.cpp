//
//  SimpleBK.cpp
//  rigor
//
//  Created by Ajmal Kunnummal on 4/24/15.
//
//

#include "SimpleBK.h"
#include "graph.h"
#include "AbstractGraph.h"

using namespace std;

SimpleBK::SimpleBK(AbstractGraph &graph) :
bk_graph(RigorGraph(graph.spixels.size(), graph.pairwise.size())) {
        
    bk_graph.add_node(graph.spixels.size());
    
    AbstractGraph::UnaryIterator unaries = graph.get_unaries(1, 10);
    for (int i = 0; i < graph.spixels.size(); unaries++, i++) {
        std::pair<edgew, edgew> unary_pair = *unaries;
        bk_graph.add_tweights(i, unary_pair.first, unary_pair.second);
    }
    
    for (pw_edge edge: graph.pairwise) {
        bk_graph.add_edge(edge.a, edge.b, edge.w, edge.w);
    }
}

vector<SimpleBK> simplebk(vector< unique_ptr<AbstractGraph> >& graphs){
    vector<SimpleBK> bk_graphs;
    bk_graphs.reserve(graphs.size());
    for ( auto& graph: graphs){
        bk_graphs.emplace_back(*graph);
    }
    for ( auto& graph: bk_graphs){
        graph.run();
    }
    return bk_graphs;
}