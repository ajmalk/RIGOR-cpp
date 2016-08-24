// @authors:     Ahmad Humayun
// @contact:     ahumayun@cc.gatech.edu
// @affiliation: Georgia Institute of Technology
// @date:        Fall 2013 - Summer 2014


#include "dynamicgraphs/bk_dynamicgraphs.h"


void getDistAdjParams(const int& nrhs, const mxArray* prhs[], const size_t& num_vars,
                      const std::vector<std::string>& options_str_tok, 
                      DistAdjParams& distadj_params)
{
  const size_t DISTADJ_NPARAMS = 3;

  const size_t num_graph_types = options_str_tok.size();

  /* initialize parameters for all graph sub-types */
  distadj_params.distadj.assign(num_graph_types, false);
  distadj_params.crossover_dist.assign(num_graph_types, DISTADJ_DFLT_CROSSOVER);
  distadj_params.max_s_delta.assign(num_graph_types, DISTADJ_DFLT_MAX_S);
  distadj_params.max_y_delta.assign(num_graph_types, DISTADJ_DFLT_MAX_T);
  distadj_params.desc_txts.assign(num_graph_types, "");

  bool any_distadj = false;

  /* check if user wants distance adjustment with some graph (solution for
   * middle child problem) */
  for (size_t g_idx = 0; g_idx < options_str_tok.size(); ++g_idx) {
    if (options_str_tok[g_idx].find("distadj") != std::string::npos) {
      distadj_params.distadj[g_idx] = true;
      any_distadj = true;
    }
  }

  // if distadj option is not set, don't need to run the rest of the function
  if (!any_distadj) {
    return;
  }

  // check the all pairs shortest path input vector
  std::size_t num_pairs = (num_vars * (num_vars-1)) / 2;
  if (nrhs > 9 && (!mxIsSingle(prhs[9]) || mxGetNumberOfDimensions(prhs[9]) != 2 ||
                   mxGetM(prhs[9]) != num_pairs || mxGetN(prhs[9]) != 1)) {
    mexErrMsgTxt("9th argument (all_pairs_shrts_paths) should be an Nx(N-1)/2 "
                 "float vector");
  }

  if (nrhs > 9) {
    distadj_params.all_pairs_shrts_paths = (disttype*)mxGetData(prhs[9]);
  } else {
    mexErrMsgTxt("To run distadjtest or distadj option you need to provide"
                 " all_pairs_shrts_paths");
  }
  
  if (nrhs > 10) {
    if (!mxIsDouble(prhs[10]) || !(mxGetNumberOfElements(prhs[10]) == DISTADJ_NPARAMS || 
        mxGetNumberOfElements(prhs[10]) == DISTADJ_NPARAMS*num_graph_types)) {
      mexErrMsgTxt((boost::format("Extra parameters for distadj should be of "
                                 "double type and should be size %d or %d x "
                                 "#graph sub-types") % DISTADJ_NPARAMS 
                                 % DISTADJ_NPARAMS).str().c_str());
    }

    if (mxGetNumberOfElements(prhs[10]) == DISTADJ_NPARAMS) {
      distadj_params.crossover_dist.assign(num_graph_types, mxGetPr(prhs[10])[0]);
      distadj_params.max_s_delta.assign(num_graph_types,    mxGetPr(prhs[10])[1]);
      distadj_params.max_y_delta.assign(num_graph_types,    mxGetPr(prhs[10])[2]);
    } else {
      for (size_t g_idx = 0; g_idx < num_graph_types; ++g_idx) {
        size_t idx = g_idx*DISTADJ_NPARAMS;
        distadj_params.crossover_dist[g_idx] = mxGetPr(prhs[10])[idx];
        distadj_params.max_s_delta[g_idx]    = mxGetPr(prhs[10])[idx+1];
        distadj_params.max_y_delta[g_idx]    = mxGetPr(prhs[10])[idx+2];
      }
    }
  }

  for (size_t g_idx = 0; g_idx < num_graph_types; ++g_idx) {
    if (distadj_params.distadj[g_idx]) {
      distadj_params.desc_txts[g_idx] = (boost::format(", †: %.2f %.2e %.2e") 
                                % distadj_params.crossover_dist[g_idx] 
                                % distadj_params.max_s_delta[g_idx] 
                                % distadj_params.max_y_delta[g_idx]).str();
    }
  }
}


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  /* passed 11 variables (N is the number of variables, S is the number of
   * seeds). The last 3 arguments are optional:
   *
   *  prhs[0] (nonlambda_s): is a NxS matrix which is a gives the non-lambda
   *                         unaries from s to v_i
   *  prhs[1] (nonlambda_t): is a NxS matrix which is a gives the non-lambda
   *                         unaries from v_i to t
   *  prhs[2] (lambda_s): is a NxS matrix which is a gives the lambda unaries
   *                         from s to v_i. Parametric cuts would need to
   *                         need to take place for using different lambdas
   *  prhs[3] (lambda_t): is a NxS matrix which is a gives the lambda unaries
   *                         from v_i to t. Parametric cuts would need to
   *                         need to take place for using different lambdas
   *  prhs[4] (lambda_range): is a 1xL vector giving the lambda values used
   *                         to parameterize the graph over lambda_s and
   *                         lambda_t
   *  prhs[5] (pairwise_edges): is a Px3 matrix of pairwise edges. The first
   *                         two columns give the node ids to which the edge
   *                         is adjacent to. The third column gives the
   *                         capacity on the edge. Note that the node ids are
   *                         offset by +1
   *  prhs[6] (graph_type_start_idx): Are the starting indices (1+) for all
   *                         different graph sub-type 
   *  prhs[7] (options_str): is a string specifying options for running
   *                         maxflow
   *  prhs[8] (print_mode):  if set to false, the maxflow code does not print
   *                         any descriptive messages (just errors).
   *  prhs[9] (all_pairs_shrts_paths): is a (N*(N-1))/2 x 1 float vector,
   *                         all pair shortest paths in the graph. It can be
   *                         as a distance matrix for guiding parametric min
   *                         cut.
   *  prhs[10] (distadj_params): are the parameters to use for distance 
   *                         adjustment of unaries for parametric min-cut.
   *                         They are important for avoiding the middle child
   *                         problem.
   */

  unarycaptype *nonlambda_s, *nonlambda_t, *lambda_s, *lambda_t;
  lambdaparamtype *lambda_range;
  pairwisecaptype *pairwise_edges;
  graphtypeidxtype *graph_type_start_idx;
  std::string options_str;
  bool print_mode = true;
  DistAdjParams distadj_params;

  size_t num_vars, num_seeds, num_params, num_pairwise_edges, num_graph_types;

  /* verify number of inputs */
  if (nrhs < 8 || nrhs > 11) {
    mexErrMsgTxt("Needs 8 to 11 input variables: \n"
                 "bk_dynamicgraphs_mex(nonlambda_s, nonlambda_t, lambda_s, "
                 "lambda_t, lambda_range, pairwise_edges, "
                 "graph_type_start_idx, options_str, [print_mode], "
                 "[all_pairs_shrts_paths], [distadj parameters])\n");
  }
  /* validate the size and type of inputs */
  if (!mxIsDouble(prhs[0]) || mxGetNumberOfDimensions(prhs[0]) != 2) {
    mexErrMsgTxt("nonlambda_s should be a two dimensional double matrix");
  }
  num_vars = mxGetM(prhs[0]);
  num_seeds = mxGetN(prhs[0]);

  if (!mxIsDouble(prhs[1]) || mxGetNumberOfDimensions(prhs[1]) != 2 ||
      mxGetM(prhs[1]) != num_vars || mxGetN(prhs[1]) != num_seeds) {
    mexErrMsgTxt("nonlambda_t should be a two dimensional double matrix with "
                 "same size as nonlambda_s");
  }
  if (!mxIsDouble(prhs[2]) || mxGetNumberOfDimensions(prhs[2]) != 2 ||
      mxGetM(prhs[2]) != num_vars || mxGetN(prhs[2]) != num_seeds) {
    mexErrMsgTxt("lambda_s should be a two dimensional double matrix with "
                 "same size as nonlambda_s");
  }
  if (!mxIsDouble(prhs[3]) || mxGetNumberOfDimensions(prhs[3]) != 2 ||
      mxGetM(prhs[3]) != num_vars || mxGetN(prhs[3]) != num_seeds) {
    mexErrMsgTxt("lambda_t should be a two dimensional double matrix with "
                 "same size as nonlambda_s");
  }
  if (!mxIsDouble(prhs[4]) || mxGetNumberOfDimensions(prhs[4]) != 2 ||
      mxGetM(prhs[4]) != 1) {
    mexErrMsgTxt("lambda_range should be a 1 dimensional double row vector");
  }
  if (!mxIsDouble(prhs[5]) || mxGetNumberOfDimensions(prhs[5]) != 2 ||
      mxGetN(prhs[5]) != 3) {
    mexErrMsgTxt("pairwise_edges should be a 2 dimensional double matrix with "
                 "3 columns");
  }
  if (!mxIsDouble(prhs[6]) || mxGetNumberOfDimensions(prhs[6]) != 2 ||
      mxGetM(prhs[6]) != 1) {
    mexErrMsgTxt("graph_type_start_idx should be a 1 dimensional double row "
                 "vector");
  }
  if (!mxIsChar(prhs[7])) {
    mexErrMsgTxt("8th argument (options_str) should be a string");
  }
  if (nrhs > 8) {
    if (!mxIsLogical(prhs[8]) || (mxGetNumberOfElements(prhs[8]) != 1)) {
      mexErrMsgTxt("9th argument (print_mode) should be a single boolean value");
    } else {
      print_mode = mxGetLogicals(prhs[8])[0];
    }
  }
  if (nlhs > 2) {
    mexErrMsgTxt("Outputs only 2 variable\n");
  }

  nonlambda_s = (unarycaptype*)mxGetData(prhs[0]);
  nonlambda_t = (unarycaptype*)mxGetData(prhs[1]);
  lambda_s = (unarycaptype*)mxGetData(prhs[2]);
  lambda_t = (unarycaptype*)mxGetData(prhs[3]);

  lambda_range = (lambdaparamtype*)mxGetData(prhs[4]);
  num_params = mxGetN(prhs[4]);

  pairwise_edges = (pairwisecaptype*)mxGetData(prhs[5]);
  num_pairwise_edges = mxGetM(prhs[5]);

  graph_type_start_idx = (graphtypeidxtype*)mxGetData(prhs[6]);
  num_graph_types = mxGetN(prhs[6]);

  options_str = mxArrayToString(prhs[7]);

  /* get the options for each graph type by tokenizing the options_str */
  std::vector<std::string> options_str_tok;
  boost::split(options_str_tok, options_str, boost::is_any_of("|"));

  if (options_str_tok.size() > num_graph_types) {
    mexErrMsgTxt("Options strings should have same number of tokens or less "
                 "than the number graph sub-types");
  }
  /*
  for (size_t n = 0; n < num_graph_types; ++n)
    std::cout << graph_type_start_idx[n] << "\n";
  std::cout << std::endl;
  for (size_t n = 0; n < options_str_tok.size(); ++n)
    std::cout << "\"" << options_str_tok[n] << "\"\n";
  std::cout << std::endl;
  */

  /*
  std::cout << num_vars << ", " << num_seeds << std::endl;
  std::cout << num_pairwise_edges << std::endl;

  for (int i=0; i < num_params; ++i)
    std::cout << boost::format("%.4f, ") % lambda_range[i];
  std::cout << std::endl;
  */

  // get the "distance adjusted" / "middle child" parameters
  getDistAdjParams(nrhs, prhs, num_vars, options_str_tok, distadj_params);

  // initialize result to all zeros
  plhs[0] = mxCreateNumericMatrix(num_vars, num_seeds, mxUINT16_CLASS,
                                  mxREAL);
  resulttype* cuts = (resulttype*)mxGetData(plhs[0]);


  std::vector<metainfotype> metainfo;


  /* check if user wants parallelization or not (only first option dictates if
   * there will be parallelism) */
  bool parallel = true;
  if (options_str_tok[0].find("noparallel") != std::string::npos) {
    parallel = false;
  }

  /* check if the user wants to schedule the lambdas in reverse */
  std::vector<bool> rev_lambda(num_graph_types, false);
  for (size_t graph_idx = 0; graph_idx < options_str_tok.size(); ++graph_idx) {
    if (options_str_tok[graph_idx].find("rev") != std::string::npos) {
      rev_lambda[graph_idx] = true;
    }
  }

  /* run max flow according to options (the first option decides the method) */
  if (options_str_tok[0].find("nodynamic") != std::string::npos) {
    nodynamic_param_maxflow_allseeds(num_seeds, num_vars, num_pairwise_edges,
                                     num_params, num_graph_types, nonlambda_s,
                                     nonlambda_t, lambda_s, lambda_t,
                                     lambda_range, pairwise_edges,
                                     graph_type_start_idx, cuts, &metainfo,
                                     parallel, print_mode);
  } else if (options_str_tok[0].find("distadjtest") != std::string::npos) {
    distadj_param_maxflow_allseeds(num_seeds, num_vars, num_pairwise_edges,
                                   num_params, num_graph_types, nonlambda_s,
                                   nonlambda_t, lambda_s, lambda_t,
                                   lambda_range, pairwise_edges,
                                   &distadj_params, graph_type_start_idx, 
                                   cuts, &metainfo, parallel, print_mode, 
                                   true);
  } else if (options_str_tok[0].find("kohli") != std::string::npos) {
    /* check if user wants to reuse the trees or just want to reuse the graph */
    std::vector<bool> reusetrees(num_graph_types, false);
    for (size_t graph_idx = 0; graph_idx < options_str_tok.size(); ++graph_idx) {
      if (options_str_tok[graph_idx].find("reusetrees") != std::string::npos) {
        reusetrees[graph_idx] = true;
      }
    }

    kohli_param_maxflow_allseeds(num_seeds, num_vars, num_pairwise_edges,
                                 num_params, num_graph_types, nonlambda_s,
                                 nonlambda_t, lambda_s, lambda_t, lambda_range,
                                 pairwise_edges, &distadj_params, 
                                 graph_type_start_idx, cuts, &metainfo, 
                                 reusetrees, rev_lambda, parallel, print_mode);
  } else if (options_str_tok[0].find("multiseed") != std::string::npos) {
    /* check if user wants to build precomputation graph from the opposite end
     * of the lambda spectrum */
    std::vector<bool> opp_precomp_lambda(num_graph_types, false);
    for (size_t graph_idx = 0; graph_idx < options_str_tok.size(); ++graph_idx) {
      if (options_str_tok[graph_idx].find("opplambda") != std::string::npos) {
        opp_precomp_lambda[graph_idx] = true;
      }
    }

    multiseeddyn_param_maxflow_allseeds(num_seeds, num_vars, num_pairwise_edges,
                                        num_params, num_graph_types,
                                        nonlambda_s, nonlambda_t, lambda_s,
                                        lambda_t, lambda_range, pairwise_edges,
                                        &distadj_params, graph_type_start_idx, 
                                        cuts, &metainfo, rev_lambda, 
                                        opp_precomp_lambda, parallel, print_mode);
  } else {
    mexErrMsgTxt((boost::format("'%s' is not a valid option") %
                  options_str).str().c_str());
  }

  /* if the output asks for meta info */
  if (nlhs == 2) {
    if (metainfo.empty())
      /* if no meta info returned by max flow function */
      plhs[1] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
    else {
      /* move meta info to output matrix */
      plhs[1] = mxCreateNumericMatrix(metainfo.size(), metainfo[0].size(), mxSINGLE_CLASS, mxREAL);
      metainfosingletype* metainfo_out = (metainfosingletype*)mxGetData(plhs[1]);

      size_t indx = 0;
      for (int idx=0; idx < metainfo[0].size(); ++idx) {
        for (int info_idx=0; info_idx < metainfo.size(); ++info_idx, ++indx) {
          metainfo_out[indx] = metainfo[info_idx][idx];
        }
      }
    }
  }
}
