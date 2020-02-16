/**
 * @file graphMsa.cpp
 * @brief implements class for representing the graph of msa problem
 * @author Sebastian Schenker
 **/

#include "graphMsa.h"
#include "graphVariable.h"
#include <assert.h>
#include <cstddef>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace std;

graphMsa::graphMsa(const vector<string>& input) {
  
  //assert( !input.empty() );
  numberOfStrings_ = input.size();
  
  // create edge and arc variables and entries in indices maps
  for(size_t i=0; i<numberOfStrings_-1; ++i) // first string
    for(size_t j=i+1; j<numberOfStrings_; ++j) // second string
    {
      size_t strI_length = input[i].length();
      size_t strJ_length = input[j].length();
      
      for(size_t nodeI=0; nodeI<strI_length; ++nodeI) //letters in first string
      {
	// create arc variables arc(nodeI,nodeK,i,j)
	for(size_t nodeK=nodeI; nodeK<strI_length; ++nodeK)
	  graphArcs_.push_back( arc(nodeI,nodeK,i,j) );

	// create edge variables edge(nodeI,nodeJ,i,j)
	for(size_t nodeJ=0; nodeJ<strJ_length; ++nodeJ) //letters in second string
	  graphEdges_.push_back( edge(nodeI,nodeJ,i,j) );
      }
      
      // create arc variables arc(nodeJ,nodeK,j,i)
      for(size_t nodeJ=0; nodeJ<strJ_length; ++nodeJ)
	for(size_t nodeK=nodeJ; nodeK<strJ_length; ++nodeK)
	  graphArcs_.push_back( arc(nodeJ,nodeK,j,i) );
    }


  // create entries in edgesBetween
  size_t currentI = graphEdges_[0].strI(),
    currentJ = graphEdges_[0].strJ(),
    begin = 0,
    end = graphEdges_.size();
  for(size_t i=0; i<end; ++i)
  {
    if( graphEdges_[i].strI() != currentI || graphEdges_[i].strJ() != currentJ )
    {
      edgesBetween_[ indexPair(currentI,currentJ) ] = indexPair(begin,i);
      begin = i;
      currentI = graphEdges_[i].strI();
      currentJ = graphEdges_[i].strJ();
    }
  }
  edgesBetween_[ indexPair(currentI,currentJ) ] = indexPair(begin,end);

  // create entries in arcsBetween
  currentI = graphArcs_[0].strI();
  currentJ = graphArcs_[0].strJ();
  begin = 0;
  end = graphArcs_.size();
  for(size_t i=0; i<end; ++i)
  {
    if( graphArcs_[i].strI() != currentI || graphArcs_[i].strJ() != currentJ )
    {
      arcsBetween_[ indexPair(currentI,currentJ) ] = indexPair(begin,i);
      begin = i;
      currentI = graphArcs_[i].strI();
      currentJ = graphArcs_[i].strJ();
    }
  }
  arcsBetween_[ indexPair(currentI,currentJ) ] = indexPair(begin,end);
  
}

void graphMsa::printEdges(ostream& out) const
{
  out << "Edges:\n";
  edges::const_iterator eEnd(graphEdges_.end());
  for(edges::const_iterator it=graphEdges_.begin(); it!=eEnd; ++it)
    out << *it << "\n";
}

void graphMsa::printArcs(ostream& out) const
{
  out << "Arcs:\n";
  arcs::const_iterator aEnd(graphArcs_.end());
  for(arcs::const_iterator it=graphArcs_.begin(); it!=aEnd; ++it)
    out << *it << "\n";
}

ostream& graphMsa::print(ostream& out) const
{
  printEdges(out);
  printArcs(out);

  out << "Indices for edges between: \n";
  edgeIndices::const_iterator eIndsEnd(edgesBetween_.end());
  for(edgeIndices::const_iterator it=edgesBetween_.begin(); it!=eIndsEnd; ++it)
    out << (it->first).first << " - " << (it->first).second << " = (" << (it->second).first << "," << (it->second).second << ")\n";
  
  cout << "Indices for arcs between: \n";
  arcIndices::const_iterator aIndsEnd(arcsBetween_.end());
  for(arcIndices::const_iterator it=arcsBetween_.begin(); it!=aIndsEnd; ++it)
    out << (it->first).first << " - " << (it->first).second << " = (" << (it->second).first << "," << (it->second).second << ")\n";
  return out;
}

ostream& operator<<(ostream& out, const graphMsa& gr)
{
  return gr.print(out);
}

const graphMsa::edges& graphMsa::getEdges() const 
{
 return graphEdges_;
}

const graphMsa::arcs& graphMsa::getArcs() const 
{
  return graphArcs_;
}

size_t graphMsa::numberOfEdges() const
{ 
  return graphEdges_.size();
}

size_t graphMsa::numberOfArcs() const
{
  return graphArcs_.size();
}

size_t graphMsa::numberOfStrings() const 
{
  return numberOfStrings_;
}

graphMsa::indexPair graphMsa::indsForEdgesBetween(size_t strI, size_t strJ) const
{
  return edgesBetween_.find(make_pair(strI,strJ))->second;
}

graphMsa::indexPair graphMsa::indsForArcsBetween(size_t strI, size_t strJ) const
{
  return arcsBetween_.find(make_pair(strI,strJ))->second;
}


