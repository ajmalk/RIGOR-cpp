// @authors:     Ahmad Humayun
// @contact:     ahumayun@cc.gatech.edu
// @affiliation: Georgia Institute of Technology
// @date:        Fall 2013 - Summer 2014

#include "bk_dynamicgraphs.h"


void gather_seed_nums(const size_t& num_seeds, const size_t& num_vars,
                      std::vector<int>& fg_seed_nums,
                      const unarycaptype* const nonlambda_s,
                      const unarycaptype* const nonlambda_t)
{
  fg_seed_nums.assign(num_seeds, 0);

  size_t idx = 0;
  for (size_t s=0; s < num_seeds; ++s) {
    /* iterate over all nodes in the current seed graph */
    for (GraphType::node_id var_i=0; var_i < num_vars; ++var_i, ++idx) {
      /* see if the current variable is a seed pixel */
      if (nonlambda_s[idx] >= INF_SEED_THRESH) {
        ++fg_seed_nums[s];
      }
    }
  }
}

void gather_seed_vars(const size_t& num_seeds, const size_t& num_vars,
                      GraphType::FGSeedsType* const fg_seeds,
                      FGSeedMapType* const fg_seed_map,
                      const unarycaptype* const nonlambda_s,
                      const unarycaptype* const nonlambda_t,
                      const bool& check_repeat_seed)
{
  fg_seeds->assign(num_seeds, GraphType::SrcSeedList());

  size_t idx = 0;
  /* iterate over all seed graphs (each graph can have multiple seed variables) */
  for (size_t s=0; s < num_seeds; ++s) {
    GraphType::SrcSeedList curr_src_seeds;
    /* std::cout << "Seed:" << s; */
    /* iterate over all nodes in the current seed graph */
    for (GraphType::node_id var_i=0; var_i < num_vars; ++var_i, ++idx) {
      /* see if the current variable is a seed pixel */
      if (nonlambda_s[idx] >= INF_SEED_THRESH) {
        GraphType::NodeCap nic(nonlambda_s[idx], 0);    // NodeCap is a std::pair (src capacity, sink capacity)
        GraphType::SrcSeedNode ssn(var_i, nic);         // SrcSeedNode is a std::pair (variable idx, unary capacities)
        curr_src_seeds.push_back(ssn);                  // SrcSeedList is a vector of SrcSeedNode(s)

        // if caller asked for fg_seed_map
        if (fg_seed_map) {
          // if it was asked to check if this variable is not a seed in another graph
          if (check_repeat_seed)
            GraphType::gassert(fg_seed_map->find(var_i) == fg_seed_map->end(), 
                               (boost::format("Seed variable %d is repeated in multiple graphs") % var_i).str());
          // fg_seed_map
          //   key: seed variable index
          //   value: std::pair (graph index, index of this variable in the list of seeds for this graph)
          (*fg_seed_map)[var_i] = FGSeedMapValType(s, curr_src_seeds.size()-1);
        }

        //std::cout << "\tNode id:" << var_i << " : " << s << std::endl;
      }
    }

    // fg_seeds is a vector of vectors. Each vector at index s stores the list
    //  seeds for graph index s
    (*fg_seeds)[s] = curr_src_seeds;
  }
}


int change_cut_for_var(const GraphType* const g, resulttype* const cuts,
                        const GraphType::node_id& var_i,
                        const unsigned int& lambda_idx, const bool REV=false,
                        NodeLst* const new_src_cut_vars=NULL)
{
  char changed = 0;

  /* supposes that cuts can only grow */
  if (REV) {
    if (cuts[var_i] == 0) {
      if (g->what_segment(var_i) == GraphType::SINK) {
        cuts[var_i] = lambda_idx+1;
        changed = -1;
      }
    }
  } else {
    if (cuts[var_i] == 0) {
      if (g->what_segment(var_i) == GraphType::SOURCE) {
        cuts[var_i] = lambda_idx+1;
        changed = 1;
        if (new_src_cut_vars)
           (*new_src_cut_vars)[var_i] = true;
      }
    }
  }

  return changed;
}

void update_cut(GraphType* const g, resulttype* const cuts,
                const unsigned int lambda_idx, unsigned int& in_src_cut,
                Block<GraphType::node_id>* const changed_list, const bool REV,
                NodeLst* const new_src_cut_vars)
{
  if (new_src_cut_vars)
    new_src_cut_vars->clear();

  if (changed_list) {
    GraphType::node_id* ptr;
    for (ptr=changed_list->ScanFirst(); ptr; ptr=changed_list->ScanNext()) {
      GraphType::node_id var_i = *ptr;
      in_src_cut += change_cut_for_var(g, cuts, var_i, lambda_idx, REV, 
                                       new_src_cut_vars);
      g->remove_from_changed_list(var_i);
    }
    changed_list->Reset();
  } else {
    /* get the variables which changed to the src side fo the cut */
    for (GraphType::node_id var_i = 0; var_i < g->get_node_num(); ++var_i) {
      in_src_cut += change_cut_for_var(g, cuts, var_i, lambda_idx, REV, 
                                       new_src_cut_vars);
    }
  }
}

inline void iterate_set_min(const GraphType::node_id& src_id, 
                            const size_t& num_vars,
                            const disttype* const all_pairs_shrts_paths,
                            disttype* const src_cut_dsts)
{
  // iterate over all nodes and set if its minimum distance
  for (std::size_t i = 0; i < num_vars; ++i) {
    // get distance to current seed
    disttype curr_d = ij2idx_v(src_id, i, all_pairs_shrts_paths, num_vars);
    // if distance to current seed less than distance in array
    if (curr_d < src_cut_dsts[i])
      src_cut_dsts[i] = curr_d;
  }
}


void init_min_src_dist(const size_t& num_vars, 
                       const DistAdjParams* const distadj_params,
                       const size_t& graph_type_idx,
                       const GraphType::SrcSeedList* const src_seed_lst,
                       disttype* const src_cut_dsts)
{
  // incase distances wouldn't be used
  if (!distadj_params || !distadj_params->distadj[graph_type_idx])
    return;

  // initialize all distances to maximum
  std::fill_n(src_cut_dsts, num_vars, std::numeric_limits<disttype>::max());
  
  // iterate over all the fg seed nodes
  for (SrcSeedList_It it = src_seed_lst->begin(); it != src_seed_lst->end(); 
        ++it) {
    const GraphType::node_id seed_id = it->first;
    // iterate over all nodes and set if its minimum distance
    iterate_set_min(seed_id, num_vars, distadj_params->all_pairs_shrts_paths, 
                    src_cut_dsts);
  }
}


/* WARNING: this function supposes that the cuts are nested, and provided in an
    increasing size order. That means that the cuts provided are in increasing
    lambda order */
void adjust_min_src_dist(const GraphType* const g,
                         const size_t& num_vars,
                         const DistAdjParams* const distadj_params,
                         const size_t& graph_type_idx,
                         const NodeLst* const new_src_cut_vars,
                         disttype* const src_cut_dsts)
{
  // incase distances wouldn't be used
  if (!distadj_params || !distadj_params->distadj[graph_type_idx])
    return;
  
  // iterate over all the variables the got accepted in the src side of the cut
  for(NodeLst_it it = new_src_cut_vars->begin(); it != new_src_cut_vars->end(); 
        ++it) {
    const GraphType::node_id new_src_id =  it->first;
    // iterate over all nodes and adjust minimum distance
    iterate_set_min(new_src_id, num_vars, distadj_params->all_pairs_shrts_paths, 
                    src_cut_dsts); 
  }
}

