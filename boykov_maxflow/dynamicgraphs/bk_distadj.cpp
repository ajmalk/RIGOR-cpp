// @authors:     Ahmad Humayun
// @contact:     ahumayun@cc.gatech.edu
// @affiliation: Georgia Institute of Technology
// @date:        Fall 2013 - Summer 2014

#include "bk_dynamicgraphs.h"


#include <fstream>

void distadj_param_maxflow(const size_t& num_vars,
                           const size_t& num_pairwise_edges,
                           const size_t& num_params,
                           const unarycaptype* const nonlambda_s,
                           const unarycaptype* const nonlambda_t,
                           const unarycaptype* const lambda_s,
                           const unarycaptype* const lambda_t,
                           const lambdaparamtype* const lambda_range,
                           const pairwisecaptype* const pairwise_edges,
                           const DistAdjParams* const distadj_params,
                           resulttype* const cuts,
                           secslisttype* const graphconst_times,
                           secslisttype* const maxflow_times,
                           counterlisttype* const num_growths,
                           counterlisttype* const num_augmentations,
                           counterlisttype* const num_adoptions,
                           flowlisttype* const flow_vals,
                           const int& seed_idx, const size_t& graph_type_idx,
                           const GraphType::SrcSeedList* const src_seed_lst,
                           const FGSeedMapType* const fg_seed_map,
                           const bool& debug_files)
{
  const int lambda_idx = 0;

  // This stores the set of additional variables just moved to the src side of
  // the cut.
  NodeLst new_src_cut_vars;

  // generate list of values storing the closest distance from each variable to
  // any of the source variables
  disttype* src_cut_dsts = new disttype[num_vars];
  init_min_src_dist(num_vars, distadj_params, graph_type_idx, src_seed_lst, 
                    src_cut_dsts);

  std::ofstream ofs, ofs2, ofs3;
  if (debug_files) {
    ofs.open ("test.txt", std::ofstream::out | std::ofstream::app);
    ofs << "SEED GRAPH: " << seed_idx << "\n";
    for (SrcSeedList_It it = src_seed_lst->begin(); it != src_seed_lst->end(); ++it) 
      ofs << "\t" << it->first+1 << "\n";
    for (std::size_t i = 0; i < num_vars; ++i) 
      ofs << src_cut_dsts[i] << "\n";
    ofs.close();

    ofs2.open ("cutdists.txt", std::ofstream::out | std::ofstream::app);
    ofs2 << "SEED GRAPH: " << seed_idx << "\n";
  
    ofs3.open ("unaries.txt", std::ofstream::out | std::ofstream::app);
    ofs3 << "SEED GRAPH: " << seed_idx << "\n";
  }

  const pairwisecaptype* const pairwise_edges_v = pairwise_edges + num_pairwise_edges;
  const pairwisecaptype* const pairwise_edges_cap = pairwise_edges_v + num_pairwise_edges;

  bool all_src = false;

  tbb::tick_count t0;

  for (unsigned int lambda_idx=0; lambda_idx < num_params && !all_src;
       ++lambda_idx) {
    t0 = tbb::tick_count::now();
    GraphType *g = new GraphType(num_vars, num_pairwise_edges);
    g->add_node(num_vars);

    if (debug_files) {
      ofs3 << "lambda_idx: " << lambda_idx << "\n";
    }

    /* add unary capacity edges (t-links) */
    unarycaptype s_cap, t_cap;
    for (GraphType::node_id var_i = 0; var_i < num_vars; ++var_i) {
      get_linear_param_caps(var_i, num_params, nonlambda_s, nonlambda_t,
                            lambda_s, lambda_t, lambda_range,
                            lambda_idx, distadj_params, src_cut_dsts,
                            graph_type_idx, s_cap, t_cap);

      g->add_tweights(var_i, s_cap, t_cap);

      if (debug_files) {
        ofs3 << s_cap - t_cap << "\n";
      }
    }

    /* add pairwise capacity edges */
    for (size_t i = 0; i < num_pairwise_edges; ++i) {
      g->add_edge((GraphType::node_id)pairwise_edges[i]-1,
                  (GraphType::node_id)pairwise_edges_v[i]-1,
                  pairwise_edges_cap[i], pairwise_edges_cap[i]);
    }
    double gc_time = (tbb::tick_count::now()-t0).seconds();

    /* run max-flow */
    t0 = tbb::tick_count::now();
    unarycaptype flow = g->maxflow();
    //std::cout << flow << std::endl;
    flow_vals->push_back(flow);
    double mf_time = (tbb::tick_count::now()-t0).seconds();

    all_src = true;


    new_src_cut_vars.clear();

    /* get the variables which changed to the src side fo the cut */
    for (GraphType::node_id var_i = 0; var_i < num_vars; ++var_i) {
      if (cuts[var_i] == 0) {
        if (g->what_segment(var_i) == GraphType::SOURCE) {
          cuts[var_i] = lambda_idx+1;
          // add to the list of new variables added to the cut
          new_src_cut_vars[var_i] = true;
        } else {
          all_src = false;
        }
      } else {
        GraphType::gassert(g->what_segment(var_i)==GraphType::SOURCE, "Nesting broken");
      }
    }

    adjust_min_src_dist(g, num_vars, distadj_params, graph_type_idx, 
                        &new_src_cut_vars, src_cut_dsts);

    if (debug_files && !new_src_cut_vars.empty()) {
      ofs2 << "lambda_idx: " << lambda_idx << "\n";
      for(NodeLst_it it = new_src_cut_vars.begin(); it != new_src_cut_vars.end(); ++it)  
        ofs2 << "\t" << it->first+1 << "\n";
      for (std::size_t i = 0; i < num_vars; ++i) 
        ofs2 << src_cut_dsts[i] << "\n";
    }

    /* add timers to the vectors */
    graphconst_times->push_back(gc_time);
    maxflow_times->push_back(mf_time);

    /* add the number of stages BK has to run to the vectors */
    num_growths->push_back(g->get_num_growths());
    num_augmentations->push_back(g->get_num_augmentations());
    num_adoptions->push_back(g->get_num_adoptions());

    delete g;
  }

  if (debug_files) {
    ofs2.close();
    ofs3.close();
  }
  /////////////////////////////////////////////////////////////////////////
  
  /////////////////////////////////////////////////////////////////////////

  //  std::cout << "\nTimings:\n--------" << std::endl;
  //  std::cout << "Graph construction time: "
  //            << graphconst_time.format(10, "%w") << "s\n";
  //  std::cout << "Max flow time: "
  //            << maxflow_time.format(10, "%w") << "s\n";

  delete []src_cut_dsts;
}


struct ParallelDistAdjCutComputation : public ParallelCutComputation {
  const size_t num_graph_types;
  const graphtypeidxtype* const graph_type_start_idx;
  const DistAdjParams* const distadj_params;
  const GraphType::FGSeedsType* const fg_seeds;
  const FGSeedMapType* const fg_seed_map;
  const bool debug_files;

  void operator()(const tbb::blocked_range<int>& range) const {
    GraphType::SrcSeedList seed_remove_list;

    /* compute cuts for each seed using the precomputation graph */
    for (int seed_idx = range.begin(); seed_idx != range.end(); ++seed_idx) {
      size_t graph_type_idx = 
        gIdx2GSubType(seed_idx, num_graph_types, graph_type_start_idx);

      unsigned long offsets = seed_idx*num_vars;
      distadj_param_maxflow(num_vars, num_pairwise_edges, num_params,
                            nonlambda_s+offsets, nonlambda_t+offsets,
                            lambda_s+offsets, lambda_t+offsets,
                            lambda_range, pairwise_edges, 
                            distadj_params, cuts+offsets,
                            &((*graphconst_time_all)[seed_idx]),
                            &((*maxflow_time_all)[seed_idx]),
                            &((*growths_all)[seed_idx]),
                            &((*augmentations_all)[seed_idx]),
                            &((*adoptions_all)[seed_idx]),
                            &((*flowvals_all)[seed_idx]),
                            seed_idx, graph_type_idx,
                            &((*fg_seeds)[seed_idx]), fg_seed_map,
                            debug_files);
    }
  }

  ParallelDistAdjCutComputation(const size_t& num_vars,
                                const size_t& num_pairwise_edges,
                                const size_t& num_params,
                                const size_t& num_graph_types,
                                const unarycaptype* const nonlambda_s,
                                const unarycaptype* const nonlambda_t,
                                const unarycaptype* const lambda_s,
                                const unarycaptype* const lambda_t,
                                const lambdaparamtype* const lambda_range,
                                const pairwisecaptype* const pairwise_edges,
                                const DistAdjParams* const distadj_params,
                                const graphtypeidxtype* const graph_type_start_idx,
                                resulttype* const cuts,
                                std::vector<secslisttype>* const graphconst_time_all,
                                std::vector<secslisttype>* const maxflow_time_all,
                                std::vector<counterlisttype>* const growths_all,
                                std::vector<counterlisttype>* const augmentations_all,
                                std::vector<counterlisttype>* const adoptions_all,
                                std::vector<flowlisttype>* const flowvals_all,
                                const GraphType::FGSeedsType* const fg_seeds,
                                const FGSeedMapType* const fg_seed_map,
                                const bool& debug_files)
      : ParallelCutComputation(num_vars, num_pairwise_edges, num_params,
                               nonlambda_s, nonlambda_t, lambda_s, lambda_t,
                               lambda_range, pairwise_edges, cuts,
                               graphconst_time_all, maxflow_time_all,
                               growths_all, augmentations_all, adoptions_all,
                               flowvals_all), 
        num_graph_types(num_graph_types), distadj_params(distadj_params),
        graph_type_start_idx(graph_type_start_idx), fg_seeds(fg_seeds), 
        fg_seed_map(fg_seed_map), debug_files(debug_files) {
  }
};


void distadj_param_maxflow_allseeds(const size_t& num_seeds,
                                    const size_t& num_vars,
                                    const size_t& num_pairwise_edges,
                                    const size_t& num_params,
                                    const size_t& num_graph_types,
                                    const unarycaptype* const nonlambda_s,
                                    const unarycaptype* const nonlambda_t,
                                    const unarycaptype* const lambda_s,
                                    const unarycaptype* const lambda_t,
                                    const lambdaparamtype* const lambda_range,
                                    const pairwisecaptype* const pairwise_edges,
                                    const DistAdjParams* const distadj_params,
                                    const graphtypeidxtype* const graph_type_start_idx,
                                    resulttype* const cuts,
                                    std::vector<metainfotype>* const metainfo,
                                    const bool& parallel, const bool& print_mode, 
                                    const bool& debug_files)
{
  /* print the settings for the max flows we are going to perform */
  if (print_mode) {
    std::cout << "dist.adj BK ||=" << parallel << " L=" << num_params << " ";
    int curr_num_seeds;
    for (size_t i=0; i < num_graph_types; ++i) {
      if (i < num_graph_types-1)
        curr_num_seeds = graph_type_start_idx[i+1] - graph_type_start_idx[i];
      else
        curr_num_seeds = num_seeds - graph_type_start_idx[i] + 1;

      std::string postfix = distadj_params ? distadj_params->desc_txts[i] : "";

      std::cout << " [#seeds=" << curr_num_seeds
                << postfix << "]";
    }
  }

  if (debug_files) {
    std::ofstream ofs;
    ofs.open ("test.txt", std::ofstream::out | std::ofstream::trunc);
    ofs << "dist.adj BK ||=" << parallel << ", #seeds=" << num_seeds << "\n";
    ofs.close();
    std::ofstream ofs2;
    ofs2.open ("cutdists.txt", std::ofstream::out | std::ofstream::trunc);
    ofs2 << "dist.adj BK ||=" << parallel << ", #seeds=" << num_seeds << "\n";
    ofs2.close();
    std::ofstream ofs3;
    ofs3.open ("unaries.txt", std::ofstream::out | std::ofstream::trunc);
    ofs3 << "dist.adj BK ||=" << parallel << ", #seeds=" << num_seeds << "\n";
    ofs3.close();
  }

  std::vector<secslisttype> graphconst_time_all(num_seeds, secslisttype());
  std::vector<secslisttype> maxflow_time_all(num_seeds, secslisttype());

  std::vector<counterlisttype> growths_all(num_seeds, counterlisttype());
  std::vector<counterlisttype> augmentations_all(num_seeds, counterlisttype());
  std::vector<counterlisttype> adoptions_all(num_seeds, counterlisttype());

  std::vector<flowlisttype> flowvals_all(num_seeds, flowlisttype());


  std::vector<GraphType::FGSeedsType> fg_seeds(1, GraphType::FGSeedsType());
  //std::vector<FGSeedMapType> fg_seed_map(1, FGSeedMapType());
  gather_seed_vars(num_seeds, num_vars, &(fg_seeds[0]), NULL,
                   nonlambda_s, nonlambda_t, false);

  /* create the class which can do parrallel seed cuts */
  ParallelDistAdjCutComputation parallelcut(num_vars, num_pairwise_edges,
                                            num_params, num_graph_types,
                                            nonlambda_s, nonlambda_t, 
                                            lambda_s, lambda_t,
                                            lambda_range, pairwise_edges,
                                            distadj_params, 
                                            graph_type_start_idx, cuts, 
                                            &graphconst_time_all,
                                            &maxflow_time_all, &growths_all,
                                            &augmentations_all,
                                            &adoptions_all, &flowvals_all,
                                            &(fg_seeds[0]), NULL,
                                            debug_files);

  if (parallel)
    /* compute max-flow/min-cut in parallel */
    tbb::parallel_for(tbb::blocked_range<int>(0, num_seeds), parallelcut);
  else
    /* compute max-flow/min-cut in serial */
    parallelcut(tbb::blocked_range<int>(0, num_seeds));

  /* a meta info row for each info type */
  metainfo->assign(7, metainfotype());

  /* collate results in the output vectors */
  for (int s=0; s < num_seeds; ++s) {
    for (int i=0; i < graphconst_time_all[s].size(); ++i) {
      (*metainfo)[0].push_back(graphconst_time_all[s][i]);
      (*metainfo)[1].push_back(maxflow_time_all[s][i]);
      (*metainfo)[2].push_back(growths_all[s][i]);
      (*metainfo)[3].push_back(augmentations_all[s][i]);
      (*metainfo)[4].push_back(adoptions_all[s][i]);
      (*metainfo)[5].push_back(flowvals_all[s][i]);
      (*metainfo)[6].push_back(s);
    }
  }
}
