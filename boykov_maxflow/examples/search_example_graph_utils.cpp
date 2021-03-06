// @authors:     Ahmad Humayun
// @contact:     ahumayun@cc.gatech.edu
// @affiliation: Georgia Institute of Technology
// @date:        Fall 2013 - Summer 2014

#include "examples.h"

bool sort_pred(const EdgeType &left, const EdgeType &right) {
  if (left.first == right.first)
    return left.second < right.second;
  else
    return left.first < right.first;
}

void create_rand_ho_edges(std::vector< EdgeType >& edges,
                          const float EDGE_PROB, const int NUM_NODES_R,
                          const int NUM_NODES_C, const int MAX_RADIUS) {
  edges.clear();

  // iterate over all nodes and try to put long range edges to the neighbors
  // in a certain distance
  for (int node_id1 = 0; node_id1 < (NUM_NODES_R * NUM_NODES_C); ++node_id1) {
    int r = node_id1 % NUM_NODES_R;
    int c = node_id1 / NUM_NODES_R;
    int start_r = (r-MAX_RADIUS >= 0) ? r-MAX_RADIUS : 0;
    // for each node put edges in its neighborhood in northern hemisphere
    for (int ri = start_r; ri <= r; ++ri) {
      int start_c = (c-MAX_RADIUS >= 0) ? c-MAX_RADIUS : 0;
      int end_c = (c+MAX_RADIUS < NUM_NODES_C) ? c+MAX_RADIUS : NUM_NODES_C-1;
      for (int ci = start_c; ci <= end_c; ++ci) {
        // in case the neighbor is on the left row or is the same node itself
        // (the neighbors on the left would be covered when the node on left
        // search for neighbors to its right)
        if (ri == r && ci <= c)
          continue;
        int node_id2 = NUM_NODES_R * ci + ri;
        if (rand() <= RAND_MAX * EDGE_PROB) {
          edges.push_back(EdgeType(node_id1, node_id2));
          //std::cout << "(" << node_id1 << ", " << node_id2 << ")" << std::endl;
        }
      }
    }
  }
}

void create_rand_edges(std::vector< EdgeType >& edges,
                      const float EDGE_PROB, const int NUM_NODES_R,
                      const int NUM_NODES_C) {
  std::vector<bool> connected(NUM_NODES_R * NUM_NODES_C, false);

  const float NE_EXTRA_PROB = 0.3;

  edges.clear();

  bool created_se, created_prv_se = false;
  for (int node_id = 0; node_id < (NUM_NODES_R * NUM_NODES_C - NUM_NODES_R);
      ++node_id) {
    created_se = false;
    if (node_id % NUM_NODES_R == 0) {
      // only S, E, SE
      if (rand() <= RAND_MAX * EDGE_PROB) {
        edges.push_back(EdgeType(node_id, node_id + 1));
        connected[edges.back().first] =
            connected[edges.back().second] = true;
      }
      if (rand() <= RAND_MAX * EDGE_PROB) {
        edges.push_back(EdgeType(node_id, node_id + NUM_NODES_R));
        connected[edges.back().first] =
            connected[edges.back().second] = true;
      }
      if (rand() <= RAND_MAX * EDGE_PROB) {
        edges.push_back(EdgeType(node_id, node_id + NUM_NODES_R + 1));
        connected[edges.back().first] =
            connected[edges.back().second] = true;
        created_se = true;
      }
    } else if (node_id % NUM_NODES_R > 0
        && node_id % NUM_NODES_R < NUM_NODES_R - 1) {
      // S, NE, E, SE
      if (rand() <= RAND_MAX * EDGE_PROB) {
        edges.push_back(EdgeType(node_id, node_id + 1));
        connected[edges.back().first] =
            connected[edges.back().second] = true;
      }
      if (rand() <= RAND_MAX * (EDGE_PROB + NE_EXTRA_PROB) && !created_prv_se) {
        edges.push_back(EdgeType(node_id, node_id + NUM_NODES_R - 1));
        connected[edges.back().first] =
            connected[edges.back().second] = true;
      }
      if (rand() <= RAND_MAX * EDGE_PROB) {
        edges.push_back(EdgeType(node_id, node_id + NUM_NODES_R));
        connected[edges.back().first] =
            connected[edges.back().second] = true;
      }
      if (rand() <= RAND_MAX * EDGE_PROB) {
        edges.push_back(EdgeType(node_id, node_id + NUM_NODES_R + 1));
        connected[edges.back().first] =
            connected[edges.back().second] = true;
        created_se = true;
      }
    } else if (node_id % NUM_NODES_R == NUM_NODES_R - 1) {
      // NE, E
      if (rand() <= RAND_MAX * (EDGE_PROB + NE_EXTRA_PROB) && !created_prv_se) {
        edges.push_back(EdgeType(node_id, node_id + NUM_NODES_R - 1));
        connected[edges.back().first] =
            connected[edges.back().second] = true;
      }
      if (rand() <= RAND_MAX * EDGE_PROB) {
        edges.push_back(EdgeType(node_id, node_id + NUM_NODES_R));
        connected[edges.back().first] =
            connected[edges.back().second] = true;
      }
    }

    created_prv_se = created_se;
  }

  /* add south edges on the last column  */
  for (int node_id = (NUM_NODES_R * NUM_NODES_C - NUM_NODES_R);
      node_id < (NUM_NODES_R * NUM_NODES_C - 1); ++node_id) {
    // S
    if (rand() <= RAND_MAX * EDGE_PROB) {
      edges.push_back(EdgeType(node_id, node_id + 1));
      connected[edges.back().first] =
          connected[edges.back().second] = true;
    }
  }

  int num_dc = 0;
  /* iterate over all nodes to see which are disconnected */
  for (unsigned int node_id = 0; node_id < connected.size(); ++node_id) {
    if (!connected[node_id]) {
      int n_e_s_w_offset[] = {-1, NUM_NODES_R, +1, -NUM_NODES_R};
      bool n_e_s_w[] = {true, true, true, true};
      if (node_id % NUM_NODES_R == 0)
        n_e_s_w[0] = false;
      if (node_id / NUM_NODES_R == NUM_NODES_C - 1)
        n_e_s_w[1] = false;
      if (node_id % NUM_NODES_R == NUM_NODES_R - 1)
        n_e_s_w[2] = false;
      if (node_id / NUM_NODES_R == 0)
        n_e_s_w[3] = false;

      for (unsigned int direc = 0; direc < 4; ++direc) {
        if (n_e_s_w[direc]) {
          edges.push_back(EdgeType(node_id, node_id + n_e_s_w_offset[direc]));
          connected[edges.back().first] =
              connected[edges.back().second] = true;
          break;
        }
      }
    }
  }

  /* sort edges */
  std::sort(edges.begin(), edges.end(), sort_pred);

  assert(num_dc == 0);
}


GraphType* create_graph(const std::vector< EdgeType >& edges,
                        const std::vector<int>& curr_edge_cap,
                        const std::vector<int>& curr_lambda_cap,
                        const std::vector<int>& fg_nodes, const int INFTY,
                        const int NUM_NODES, const bool fg_cap_inf) {
  bool is_fg;
  GraphType *g = new GraphType(NUM_NODES, edges.size());
  g->add_node(NUM_NODES);
  for (unsigned int i = 0; i < NUM_NODES; ++i) {
    is_fg = false;
    for (unsigned int j = 0; j < fg_nodes.size(); ++j) {
      if (i == fg_nodes[j]) {
        is_fg = true;
        break;
      }
    }
    if (is_fg)
      g->add_tweights(i, /* capacities */(fg_cap_inf ? INFTY : curr_lambda_cap[i]), 0);
    else
      g->add_tweights(i, /* capacities */0, curr_lambda_cap[i]);
  }

  /* capacity edges */
  for (unsigned int i = 0; i < edges.size(); ++i) {
    g->add_edge(edges[i].first, edges[i].second, curr_edge_cap[i],
                curr_edge_cap[i]);
  }

  return g;
}


std::vector<GraphType::termtype> compute_display_result(GraphType* g,
                                                        const int NUM_NODES_R,
                                                        const int NUM_NODES_C,
                                                        double& flow, const bool disp) {
  std::vector<GraphType::termtype> result(g->get_node_num());

  /*
   for (unsigned int i = 0; i < g->get_node_num(); ++i)
   std::cout << g->get_trcap(i) << std::endl;
   */
  flow = g->maxflow();
  if (disp) std::cout << "Flow = " << flow << std::endl;
  if (disp) std::cout << "Minimum cut:" << std::endl;
  for (int j = 0; j < NUM_NODES_R; ++j) {
    for (int i = j; i < g->get_node_num(); i += NUM_NODES_R) {
      if (disp) {
        if (g->what_segment(i) == GraphType::SOURCE)
          std::cout << "\tS";
        else
          std::cout << "\tT";
      }
      result[i] = g->what_segment(i);
    }
    if (disp) std::cout << std::endl;
  }

  return result;
}
