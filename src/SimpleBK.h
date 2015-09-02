//
//  SimpleBK.h
//  rigor
//
//  Created by Ajmal Kunnummal on 4/24/15.
//
//

#ifndef __rigor__SimpleBK__
#define __rigor__SimpleBK__

#include <vector>

#include "graph.h"

#include "AbstractGraph.h"

using namespace std;

typedef Graph<double, double, double> RigorGraph;

class SimpleBK: public GraphCut {
    RigorGraph bk_graph;
public:
    SimpleBK(AbstractGraph &graph);
    
    inline void run(){
        bk_graph.maxflow();
    }
    
    inline bool in_source_seg(int sp){
        return !bk_graph.what_segment(sp);
    }
};

vector<SimpleBK> simplebk(vector< unique_ptr<AbstractGraph> >& graphs);

#endif /* defined(__rigor__SimpleBK__) */
