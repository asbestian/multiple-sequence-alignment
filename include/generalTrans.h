/**
 * @file generalTrans.h
 * @brief generalized transitivity constraint
 * @author Sebastian Schenker
 **/

#ifndef GENERAL_TRANSITIVITY_CONSTRAINT
#define GENERAL_TRANSITIVITY_CONSTRAINT

#include "graphMsa.h"
#include "graphVariable.h"
#include <cstddef>
#include <lemon/static_graph.h>
#include <lemon/preflow.h>
#include <map>
#include <string>
#include <utility>
#include <vector>

class generalTrans {

 public:
  
  generalTrans(const std::vector<std::size_t>& inputSizes, const std::map<edge,double>& edgesInSolution);
  ~generalTrans();

  bool violatedTransConsExist();
  void getEdgesInViolatedTransCons(const std::vector<std::size_t>& inputSizes, graphMsa::edges& positiveEdgesInViolatedTransCons, graphMsa::edges& negativeEdgesInViolatedTransCons);

 private:
  static const double compareToOne;
  double minCutValue_;
  unsigned int sourceNode_;
  unsigned int sourceNodeString_;
  unsigned int str2_;
  unsigned int str3_;
  std::vector< unsigned int > S_2_;
  std::vector< unsigned int > S_3_;
  std::vector< std::pair<int,int> > arcs_;
  lemon::StaticDigraph graph_;
  lemon::StaticDigraph::ArcMap<double>* capacities_;
  lemon::Preflow< lemon::StaticDigraph,lemon::StaticDigraph::ArcMap<double> >* flow_;
  
  void buildGraph(const std::vector<std::size_t>& inputSizes);
  
  void addSourceNodeCapacities(const std::map<edge,double>& edgesInSolution);
  
  void addSinkNodeCapacities(const std::vector<std::size_t>& inputSizes, const std::map<edge,double>& edgesInSolution);
  
  void addRemainingArcCapacities(const std::vector<std::size_t>& inputSizes, const std::map<edge,double>& edgesInSolution);
  
  void addCapacities(const std::vector<std::size_t>& inputSizes, const std::map<edge,double>& edgesInSolution);
  
  void getNodeCuts(const std::vector<std::size_t>& inputSizes);
  
  double computeMinCutValue(const std::vector<std::size_t>& inputSizes);
  
  bool computeMaxFlow(unsigned int str1, unsigned int str2, unsigned int str3, const std::vector<std::size_t>& inputSizes, const std::map<edge,double>& edgesInSolution, bool case_str2=false, bool case_str3=false); /* case_strI=true means that strI is taken as the primary string in the transitivity inequality */
  
  
};

#endif
