/**
 * @file mixedCycle.h
 * @brief mixed cycle constraint
 * @author Sebastian Schenker
 **/

#ifndef MIXED_CYCLE_CONSTRAINT
#define MIXED_CYCLE_CONSTRAINT

#include "graphVariable.h"
#include "graphMsa.h"
#include <cstddef>
#include <lemon/static_graph.h>
#include <lemon/dijkstra.h>
#include <map>
#include <string>
#include <vector>

class mixedCycle {

 public:

  mixedCycle(const std::vector<std::size_t>& input, const std::map<edge,double>& graphWeights);
  bool violatedMixedCycleExist();
  void getEdgesInViolatedMixedCycle(graphMsa::edges& edgesInMixedCycle);
  ~mixedCycle();

  
 private:
  static const double compareToOne;
  int valueOfMixedCycle_;
  int endNodeString_;  // string number where the mixed cycle (path) in graph instance starts
  int endNode_; // node number where the mixed cycle path in graph instance ends
  std::vector< std::pair<int,int> > arcs_;
  std::vector<std::size_t> inputSizes_; // vector of length of input strings
  lemon::StaticDigraph graph_; 
  lemon::StaticDigraph::ArcMap<double>* weightMap_;
  lemon::Dijkstra<lemon::StaticDigraph,lemon::StaticDigraph::ArcMap<double> >* dijkstra_;

  void buildGraphInstance(unsigned int numberOfStrings, unsigned int numberOfNodes, const std::map<edge,double>& arcWeights);
  void computeMaximalMixedCycle(unsigned int numberOfStrings);
  
  
};

#endif

