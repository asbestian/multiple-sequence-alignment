/**
 * @file edgeClique.h
 * @brief edge clique
 * @author Sebastian Schenker
 **/

#ifndef EDGE_CLIQUE
#define EDGE_CLIQUE

#include "graphVariable.h"
#include "graphMsa.h"
#include <lemon/static_graph.h>
#include <lemon/dijkstra.h>
#include <map>
#include <utility>
#include <vector>

class edgeClique {

 public:
  
  /* case_ji indicates whether we consider E_ij(lb<->le,1<->|s^j|) [case_ji = false] 
     or E_ji(lb<->le,1<->|s^i|) [case_ji = true] ;
     if case_ji == true, then call constructor edgeClique(strI,strJ,strI_length,strJ_length,...,true) 
     with strI < strJ */
  edgeClique(int strI, int strJ, int strI_length, int strJ_length, const std::map<edge,double>& edgesInSolution, bool case_ji=false);

  /* computes longest path ( based on E_ij(lb<->le,1<->|strJ|) ) in the pairgraph */
  double computeLongestPath(int lb, int le) const;
  void getEdgesInClique(int le, graphMsa::edges& edgesInClique) const;
  
  /* compute all pairs longest path; needed for edge weights for lifted mixed cycle constraints */
  void computeEdgeWeightsForLiftedMixedCycle(std::map<edge,double>& edgeWeightsForLiftedMixedCycle);

  ~edgeClique();

 private:

  static const double epsilon;
  int strI_;
  int strJ_;
  int strI_length_;
  int strJ_length_;
  bool case_ji_;
  std::vector< std::pair<int,int> > edges_;
  lemon::StaticDigraph pairGraph_;
  lemon::StaticDigraph::ArcMap<double>* edgeWeights_;
  lemon::Dijkstra<lemon::StaticDigraph,lemon::StaticDigraph::ArcMap<double> >* dijkstra_;

  int getHorizontalLength() const;
  int getVerticalLength() const;

};

#endif
