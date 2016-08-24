// @authors:     Ahmad Humayun
// @contact:     ahumayun@cc.gatech.edu
// @affiliation: Georgia Institute of Technology
// @date:        Fall 2014


#include "block.h"
#include "tests.h"
#include <boost/format.hpp>

struct CONST_VALS {
  int NUM_NODES_R;
  int NUM_NODES_C;
  int NUM_NODES;
  int MAX_CONN_RADIUS;
  float EDGE_PROB;
  float FG_NODE_PROB;
  int MIN_INIT_SRC, MAX_INIT_SRC;
  int GRAPH_TYPE;
  int MAX_PAIRWISE_CAP;
  int MAX_UNARY_CAP;
  int MAX_UNARY_DELTA;
  int MAX_SPECIAL_UNARY;
};

void gen_init_nest_graph(const CONST_VALS& const_vals, 
                         std::vector<int>& curr_lambda_cap,
                         std::vector<int>& curr_edge_cap,
                         std::vector<int>& fg_nodes,
                         std::vector< EdgeType >& edges,
                         std::vector<GraphType::termtype>& min_cut,
                         int& num_src_cut)
{
  int num_tries = 0;

  curr_lambda_cap.assign(const_vals.NUM_NODES, 0);

  while (num_tries++ < 10) {
    GraphType *g;
    double maxflow;

    // construct a random graph with NUM_NODES_RxNUM_NODES_C nodes
    if (const_vals.GRAPH_TYPE == 1) {
      // if grid graph
      create_rand_edges(edges, const_vals.EDGE_PROB, const_vals.NUM_NODES_R, 
                        const_vals.NUM_NODES_C);
    } else {
      create_rand_ho_edges(edges, const_vals.EDGE_PROB, const_vals.NUM_NODES_R, 
                           const_vals.NUM_NODES_C, const_vals.MAX_CONN_RADIUS);
    }
    curr_edge_cap.assign(edges.size(), 0);

    const int PCAP = const_vals.MAX_PAIRWISE_CAP;
    const int UCAP = const_vals.MAX_UNARY_CAP;

    /* create random edge and unary capacities */
    /* Note unary capacities are decided now, and then later it is decided whether that
       capacity is for src or sink. If node gets in fg_nodes, it is src capacity */
    for (int i = 0; i < edges.size(); ++i) curr_edge_cap[i] = rand() % PCAP;
    for (int i = 0; i < const_vals.NUM_NODES; ++i) curr_lambda_cap[i] = rand() % UCAP;

    fg_nodes.clear();

    /* create list of fg nodes */
    for (unsigned int i = 0; i < const_vals.NUM_NODES; ++i) {
      if (rand() <= RAND_MAX * const_vals.FG_NODE_PROB) {
        fg_nodes.push_back(i);
      }
    }

   /* create the graph using the edge and unary capacities generated above */
    g = create_graph(edges, curr_edge_cap, curr_lambda_cap,
                     fg_nodes, 0, const_vals.NUM_NODES, false);
    /* compute the graph cut and store the result */
    min_cut = compute_display_result(g, const_vals.NUM_NODES_R, 
                                     const_vals.NUM_NODES_C, maxflow, false);
    delete g;

    // check if there are atleast a certain number of variables on the src side of the cut
    num_src_cut = 0;
    for (int i = 0; i < const_vals.NUM_NODES; ++i) 
      num_src_cut += (min_cut[i] == GraphType::SOURCE);

    if (const_vals.MIN_INIT_SRC <= num_src_cut &&
        num_src_cut <= const_vals.MAX_INIT_SRC)
      break;
  }
}

bool gen_test_mutation_graph(const CONST_VALS& const_vals, 
                             const std::vector<int>& curr_lambda_cap,
                             const std::vector<int>& curr_edge_cap,
                             const std::vector<bool>& is_src_cap,
                             const std::vector< EdgeType >& edges,
                             const std::vector<GraphType::termtype>& prv_min_cut,
                             std::vector<int>& new_lambda_cap,
                             std::vector<int>& new_fg_nodes,
                             std::vector<GraphType::termtype>& min_cut,
                             int& num_src_cut)
{
  new_lambda_cap.assign(const_vals.NUM_NODES, 0);
  new_fg_nodes.clear();

  GraphType *g;
  double maxflow;

  /* increase src capacities and decrease sink capacities */
  for (int i = 0; i < const_vals.NUM_NODES; ++i) {
    /* first check if the node has src or sink capacity */
    bool node_src_cap = is_src_cap[i];
    /* check if the node was on the src side of the cut */
    bool was_src = (prv_min_cut[i] == GraphType::SOURCE);
    if (was_src) {
      // if it was in src before
      // increase src capacity and decrease sink capacity
      if (node_src_cap) {
        new_lambda_cap[i] = curr_lambda_cap[i] + (rand() % const_vals.MAX_UNARY_DELTA);
        new_fg_nodes.push_back(i);
      } else {
        int new_cap = curr_lambda_cap[i] - (rand() % const_vals.MAX_UNARY_DELTA);
        new_lambda_cap[i] = (new_cap < 0) ? 0 : new_cap;
      }
    } else {
      // if it was in sink side of the cut, do whatever
      new_lambda_cap[i] = rand() % const_vals.MAX_SPECIAL_UNARY;
      if (rand() <= RAND_MAX * const_vals.FG_NODE_PROB) {
        new_fg_nodes.push_back(i);
      }
    }
  }

 /* create the graph using the edge and unary capacities generated above */
  g = create_graph(edges, curr_edge_cap, new_lambda_cap,
                   new_fg_nodes, 0, const_vals.NUM_NODES, false);
  /* compute the graph cut and store the result */
  min_cut = compute_display_result(g, const_vals.NUM_NODES_R, 
                                   const_vals.NUM_NODES_C, maxflow, false);
  delete g;

  // check if nesting property is true
  bool nesting = true;
  num_src_cut = 0;
  for (int i = 0; i < const_vals.NUM_NODES; ++i) {
    num_src_cut += (min_cut[i] == GraphType::SOURCE);

    // if in previous solution it was on the src side of the cut, it needs
    // to be again on the src side of the cut for the nesting property to 
    // hold
    if ((prv_min_cut[i] == GraphType::SOURCE) && 
        (min_cut[i] != GraphType::SOURCE)) {
      nesting = false;
    }
  }
  
  return nesting;
}

void test_nesting(const int& graph_type) {
  CONST_VALS const_vals;
  const_vals.NUM_NODES_R = 8;
  const_vals.NUM_NODES_C = 8;
  const_vals.NUM_NODES = const_vals.NUM_NODES_R * const_vals.NUM_NODES_C;

  const_vals.GRAPH_TYPE = graph_type;
  if (graph_type == 1) {
    // if grid graph
    const_vals.EDGE_PROB = 0.6;
    const_vals.MAX_PAIRWISE_CAP = 20;
    const_vals.MAX_UNARY_CAP = 20;
    const_vals.MAX_UNARY_DELTA = 15;
  } else {
    // if higher order graph
    const_vals.EDGE_PROB = 0.3;
    const_vals.MAX_CONN_RADIUS = 5;
    const_vals.MAX_PAIRWISE_CAP = 20;
    const_vals.MAX_UNARY_CAP = 140;
    const_vals.MAX_UNARY_DELTA = 110;
  }
  const_vals.MAX_SPECIAL_UNARY = 100;
  const_vals.FG_NODE_PROB = 0.5;

  const_vals.MIN_INIT_SRC = const_vals.NUM_NODES*0.1;
  const_vals.MAX_INIT_SRC = const_vals.NUM_NODES*1.0;

  srand(time(NULL));

  const int MAX_TRIES = 10000000;
  const int CPYS_PER_GRAPH = 10;
  int curr_try = 0;

  int main_num_src_cut;
  int nest_num_src_cut;

  std::vector<int> main_lambda_cap;
  std::vector<int> nest_lambda_cap;
  std::vector<int> curr_edge_cap;
  std::vector<int> main_fg_nodes;
  std::vector<int> nest_fg_nodes;
  std::vector< EdgeType > edges;
  std::vector<GraphType::termtype> main_min_cut;
  std::vector<GraphType::termtype> nest_min_cut;

  std::cout << "Testing the nesting property" << std::endl;

  /* iterate over different edge graphs */
  while (curr_try++ < MAX_TRIES) {
    gen_init_nest_graph(const_vals, main_lambda_cap,
                        curr_edge_cap, main_fg_nodes, edges,
                        main_min_cut, main_num_src_cut);

    /* convert node indices in fg_nodes to a binary vector (faster search) */
    std::vector<bool> is_src_cap(const_vals.NUM_NODES, false);
    for (int i = 0; i < main_fg_nodes.size(); ++i) {
      is_src_cap[main_fg_nodes[i]] = true;
    }

    // generate multiple changes to the current graph, and compare results
    for (std::size_t c = 0; c < CPYS_PER_GRAPH; ++c) {
      // create a mutation of the current graph, and check nesting property
      bool is_nested = 
        gen_test_mutation_graph(const_vals, main_lambda_cap,
                                curr_edge_cap, is_src_cap, edges,
                                main_min_cut, nest_lambda_cap,
                                nest_fg_nodes, nest_min_cut,
                                nest_num_src_cut);

      std::cout << boost::format("#edges %d\t (%d)\t %+d") % edges.size() % main_num_src_cut %
                   (nest_num_src_cut - main_num_src_cut) << std::endl;
      if (!is_nested || nest_num_src_cut < main_num_src_cut) {
        std::cout << "Failed!" << std::endl;
        return;
      }
      if (nest_num_src_cut > const_vals.NUM_NODES || main_num_src_cut > const_vals.NUM_NODES) {
        std::cout << "Failed!" << std::endl;
        return;
      }
      if (nest_num_src_cut > const_vals.NUM_NODES || main_num_src_cut > const_vals.NUM_NODES) {
        std::cout << "Failed!" << std::endl;
        return;
      }
      // run max-flow/min-cut

      // check the nesting of solutions
    }
  }
}
