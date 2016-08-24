//
//  Parameters.h
//  Rigor
//
//  Created by Ajmal Kunnummal on 3/2/15.
//  Copyright (c) 2015 ajmal. All rights reserved.
//

#ifndef __Rigor__Parameters__
#define __Rigor__Parameters__

#include "FuxinGraph.h"

class Parameters {
	friend class PropGenerator;

	bool debug = false;

    int num_sp = 1000;

    double compactness = 20;

    int seeds = 25;

    int seed_radius = 5;

    float scale_color_unary = 220;

    std::vector<create_graph_func> graph_types = {
//        FuxinGraphs::uniform_intern_graph,
//        FuxinGraphs::uniform_extern_graph,
//        FuxinGraphs::uniform_extern2_graph,
        FuxinGraphs::color_intern_graph,
//        FuxinGraphs::color_extern_graph,
//        FuxinGraphs::color_extern2_graph,
    };

    std::vector<PWParams> pw_params = FuxinGraphs::PWPARAMS;

    Parameters() {}

public:

    inline static Parameters getDefault() {
	    return Parameters();
	}

    inline static Parameters getDebug() {
	    Parameters param = Parameters();
	    param.debug = true;
	    param.seeds = 1;
	    param.seed_radius = 1;
	    param.num_sp = 1000;
	    return param;
	}
};

#endif /* defined(__Rigor__Parameters__) */
