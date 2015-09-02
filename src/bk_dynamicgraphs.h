// @authors:     Ahmad Humayun
// @contact:     ahumayun@cc.gatech.edu
// @affiliation: Georgia Institute of Technology
// @date:        Fall 2013 - Summer 2014

#ifndef _BK_DYNAMICGRAPHS_H_
#define _BK_DYNAMICGRAPHS_H_

#include <iostream>
#include <map>
#include <math.h>
#include <limits>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <tbb/tick_count.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#include "graph.h"
#include "AbstractGraph.h"

typedef double unarycaptype;
typedef double pairwisecaptype;
typedef double lambdaparamtype;
typedef double graphtypeidxtype;
typedef float disttype;
typedef unsigned short resulttype;

typedef float metainfosingletype;
typedef std::vector<float> metainfotype;
typedef std::vector<float> secslisttype;
typedef std::vector<unsigned int> counterlisttype;
typedef std::vector<unarycaptype> flowlisttype;

typedef Graph<pairwisecaptype, unarycaptype, unarycaptype> GraphType;

// FGSeedMapValType is an std::pair (graph index, index of this variable in the list of seeds for this graph)
typedef std::pair<int, int> FGSeedMapValType;
//   key: seed variable index
//   value: std::pair (graph index, index of this variable in the list of seeds for this graph)
// FGSeedMapType[i].first gives the graph index for which variable i is a seed
// FGSeedMapType[i].second gives the seed index in the list of seeds for the graph FGSeedMapType[i].first
typedef std::map<GraphType::node_id, FGSeedMapValType> FGSeedMapType;

typedef unsigned long long t_idx;
typedef GraphType::SrcSeedList::const_iterator SrcSeedList_It;
typedef std::map<GraphType::node_id, bool> NodeLst;
typedef NodeLst::const_iterator NodeLst_it;


const double INF_SEED_THRESH = 21475000000;

const unarycaptype DISTADJ_DFLT_CROSSOVER = 0.2;
const unarycaptype DISTADJ_DFLT_MAX_S = 1e7;
const unarycaptype DISTADJ_DFLT_MAX_T = 1e5;

/* Stores parameters for distance adjustment for all graph sub-types (solution 
 * for middle child problem) */
struct DistAdjParams {
  disttype *all_pairs_shrts_paths;
  std::vector<bool> distadj;
  std::vector<double> crossover_dist;
  std::vector<unarycaptype> max_s_delta;
  std::vector<unarycaptype> max_y_delta;
  std::vector<std::string> desc_txts;

  DistAdjParams() {
    all_pairs_shrts_paths = NULL;
  }
};

// This function returns the graph sub type index, given the seed graph index
inline size_t gIdx2GSubType(const size_t& seed_idx, const size_t& num_graph_types, 
                            const graphtypeidxtype* const graph_type_start_idx)
{
  /* identify which graph type this current seed belongs to */
  size_t graph_type_idx = 0;
  for (; graph_type_idx < num_graph_types-1; ++graph_type_idx) {
    graphtypeidxtype end_idx = graph_type_start_idx[graph_type_idx+1]-1;
    if (seed_idx < end_idx)
      break;
  }
  return graph_type_idx;
}


void multiseeddyn_param_maxflow();

void multiseeddyn_param_maxflow_allseeds();

void kohli_param_maxflow();

void kohli_param_maxflow_allseeds();

std::vector<bool> nodynamic_param_maxflow(std::vector< std::unique_ptr<AbstractGraph> > &graph);

void nodynamic_param_maxflow_allseeds();

void distadj_param_maxflow();

void distadj_param_maxflow_allseeds();

void update_cut(GraphType* const g, resulttype* const cuts,
                const unsigned int lambda_idx, unsigned int& in_src_cut,
                Block<GraphType::node_id>* const changed_list=NULL,
                const bool REV=false, NodeLst* const new_src_cut_vars=NULL);

void gather_seed_nums(const size_t& num_seeds, const size_t& num_vars,
                      std::vector<int>& fg_seed_nums,
                      const unarycaptype* const nonlambda_s,
                      const unarycaptype* const nonlambda_t);

void gather_seed_vars(const size_t& num_seeds, const size_t& num_vars,
                      GraphType::FGSeedsType* const fg_seeds,
                      FGSeedMapType* const fg_seed_map,
                      const unarycaptype* const nonlambda_s,
                      const unarycaptype* const nonlambda_t,
                      const bool& check_repeat_seed=true);

// function to get index in all_pairs_shrts_paths for d(i,j)
inline t_idx ij2idx(const t_idx& i, const t_idx& j, const std::size_t N)
{
  return i*N - ((i+1)*i)/2 + j - i - 1;
}

inline disttype ij2idx_v(const t_idx& i, const t_idx& j, 
                         const disttype* const all_pairs_shrts_paths, 
                         const std::size_t N)
{
  if (i == j)
    return 0;
  if (i < j)
    return all_pairs_shrts_paths[ij2idx(i, j, N)];
  else
    return all_pairs_shrts_paths[ij2idx(j, i, N)];
}

inline void iterate_set_min(const GraphType::node_id& src_id, 
                            const size_t& num_vars,
                            const disttype* const all_pairs_shrts_paths,
                            disttype* const src_cut_dsts);

void init_min_src_dist(const size_t& num_vars, 
                       const DistAdjParams* const distadj_params,
                       const size_t& graph_type_idx,
                       const GraphType::SrcSeedList* const src_seed_lst,
                       disttype* const src_cut_dsts);

/* WARNING: this function supposes that the cuts are nested, and provided in an
    increasing size order. That means that the cuts provided are in increasing
    lambda order */
void adjust_min_src_dist(const GraphType* const g,
                         const size_t& num_vars,
                         const DistAdjParams* const distadj_params,
                         const size_t& graph_type_idx,
                         const NodeLst* const new_src_cut_vars,
                         disttype* const src_cut_dsts);

struct ParallelCutComputation {
  const size_t num_vars, num_pairwise_edges, num_params;
  const unarycaptype* const nonlambda_s;
  const unarycaptype* const nonlambda_t;
  const unarycaptype* const lambda_s;
  const unarycaptype* const lambda_t;
  const lambdaparamtype* const lambda_range;
  const pairwisecaptype* const pairwise_edges;

  resulttype* const cuts;

  std::vector<secslisttype>* const graphconst_time_all;
  std::vector<secslisttype>* const maxflow_time_all;

  std::vector<counterlisttype>* const growths_all;
  std::vector<counterlisttype>* const augmentations_all;
  std::vector<counterlisttype>* const adoptions_all;

  std::vector<flowlisttype>* const flowvals_all;

  // The operator for doing the parallelization
  virtual void operator()(const tbb::blocked_range<int>& range) const = 0;

  ParallelCutComputation(const size_t& num_vars,
                         const size_t& num_pairwise_edges,
                         const size_t& num_params,
                         const unarycaptype* const nonlambda_s,
                         const unarycaptype* const nonlambda_t,
                         const unarycaptype* const lambda_s,
                         const unarycaptype* const lambda_t,
                         const lambdaparamtype* const lambda_range,
                         const pairwisecaptype* const pairwise_edges,
                         resulttype* cuts,
                         std::vector<secslisttype>* const graphconst_time_all,
                         std::vector<secslisttype>* const maxflow_time_all,
                         std::vector<counterlisttype>* const growths_all,
                         std::vector<counterlisttype>* const augmentations_all,
                         std::vector<counterlisttype>* const adoptions_all,
                         std::vector<flowlisttype>* const flowvals_all)
      : num_vars(num_vars), num_pairwise_edges(num_pairwise_edges),
        num_params(num_params), nonlambda_s(nonlambda_s),
        nonlambda_t(nonlambda_t), lambda_s(lambda_s), lambda_t(lambda_t),
        lambda_range(lambda_range), pairwise_edges(pairwise_edges),
        cuts(cuts), graphconst_time_all(graphconst_time_all),
        maxflow_time_all(maxflow_time_all), growths_all(growths_all),
        augmentations_all(augmentations_all), adoptions_all(adoptions_all),
        flowvals_all(flowvals_all) {
  }

  virtual ~ParallelCutComputation() {}
};

#endif // _BK_DYNAMICGRAPHS_H_
