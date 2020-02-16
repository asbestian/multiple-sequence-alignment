/**
 * @file graphMsa.h
 * @brief class for representing graph of msa problem
 * @author Sebastian Schenker
 **/

#ifndef GRAPH_MSA
#define GRAPH_MSA

#include "graphVariable.h"
#include <cstddef>
#include <iosfwd>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

/// graphMsa

/* This class is used for creating and representing the graph for the
   multiple sequence alignment problem. In other words, edges between
   nodes and gap arcs are created. 
*/

class graphMsa {

 public:

  typedef std::vector<edge> edges;
  typedef std::vector<arc> arcs;
  typedef std::pair<std::size_t,std::size_t> indexPair;

  typedef std::map< indexPair, indexPair > edgeIndices;
  typedef std::map< indexPair, indexPair > arcIndices;
  
  explicit graphMsa(const std::vector<std::string>& inputStrings);

  const edges& getEdges() const;
  const arcs& getArcs() const;
  std::size_t numberOfStrings() const;
  std::size_t numberOfEdges() const;
  std::size_t numberOfArcs() const;
  indexPair indsForEdgesBetween(std::size_t stringI, std::size_t stringJ) const;
  indexPair indsForArcsBetween(std::size_t stringI, std::size_t stringJ) const;
  void printEdges(std::ostream& out) const;
  void printArcs(std::ostream& out) const;
  std::ostream& print(std::ostream& out) const;

 private:
  std::size_t numberOfStrings_;
  edges graphEdges_;
  arcs graphArcs_;
  edgeIndices edgesBetween_;
  arcIndices arcsBetween_;

};

std::ostream& operator<<(std::ostream& out, const graphMsa& gr);

#endif
