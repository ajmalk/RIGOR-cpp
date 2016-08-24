// @authors:     Ahmad Humayun
// @contact:     ahumayun@cc.gatech.edu
// @affiliation: Georgia Institute of Technology
// @date:        Fall 2013 - Summer 2014

#include "bk_dynamicgraphs.h"

using namespace std;

vector<bool> nodynamic_param_maxflow(vector< unique_ptr<AbstractGraph> > &graphs) {
    const int num_sp = graphs[0]->image->getNumPixels();
    vector<bool> cuts;
    cuts.reserve(num_sp * graphs[0]->lambdas.size() * graphs.size());
    
//    cout << graphs.size() << endl;
    
    for ( auto &graph: graphs) {
//        cout << graph->type<< endl;
        for (int i = 0; i < graph->lambdas.size(); i++) {
    //        cout << lambda.first << " " << lambda.second << " ";
            GraphType bk_graph = GraphType(num_sp, graph->getNumPairwise());
            bk_graph.add_node(num_sp);
            
            /* add unary capacity edges (t-links) */
            AbstractGraph::UnaryIterator unaries = graph->get_unaries(graph->lambdas[i], graph->lambdas[graph->lambdas.size() - i - 1]);
            for (int i = 0; i < num_sp; unaries++, i++) {
                std::pair<edgew, edgew> unary_pair = *unaries;
                bk_graph.add_tweights(i, unary_pair.first, unary_pair.second);
            }
            
            /* add pairwise capacity edges */
            AbstractGraph::PairwiseIterator pairwiseIt = graph->get_pairwise_it();
            for (int i = 0; i < graph->getNumPairwise(); pairwiseIt++, i++) {
                pw_edge edge = *pairwiseIt;
                bk_graph.add_edge(edge.a, edge.b, edge.w, edge.w);
            }
            
            /* run max-flow */
            bk_graph.maxflow();
            
            bool all_src = true;
            bool all_eq = true;
            
            /* get the variables which changed to the src side fo the cut */
            
            if (cuts.size() >= num_sp){
                for (auto i = 0; i < num_sp; i++) {
                    bool in_src = !bk_graph.what_segment(i);
                    cuts.push_back(in_src);
                    all_src &= in_src;
                    all_eq &= cuts[cuts.size() - num_sp - 1] == in_src;
    //                cout << in_src;
                }
            } else {
                // handle the first one
                all_eq = false;
                for (auto i = 0; i < num_sp; i++) {
                    bool in_src = !bk_graph.what_segment(i);
    //                cout << in_src;
                    cuts.push_back(in_src);
                    all_src &= in_src;
                }
            }
    //        cout << " " << all_eq << endl;
            
//            if(all_src) {
//                cuts.erase(cuts.end() - num_sp, cuts.end());
//                break;
//            }
            
            if(all_eq) {
                cuts.erase(cuts.end() - num_sp, cuts.end());
            }
        }
    }
    return cuts;
}
