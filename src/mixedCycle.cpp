/**
 * @file mixedCycle.cpp
 * @brief Implementation of mixed cycle constraint class
 * @author Sebastian Schenker
 **/

#include "mixedCycle.h"
#include "graphVariable.h"
#include "graphMsa.h"
#include <cassert>
#include <cstddef>
#include <map>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

using namespace std;
using namespace lemon;

const double mixedCycle::compareToOne = 1.0 + 1.0e-09;

mixedCycle::mixedCycle(const vector<size_t>& inputSizes, const map<edge,double>& arcWeights)
  : valueOfMixedCycle_(0), endNodeString_(-1), endNode_(-1), inputSizes_(inputSizes), weightMap_(0), dijkstra_(0)
{
  unsigned int numberOfNodes = accumulate(inputSizes.begin(),inputSizes.end(),0);
  buildGraphInstance(inputSizes.size(),numberOfNodes,arcWeights);
  computeMaximalMixedCycle(inputSizes.size());
}

mixedCycle::~mixedCycle()
{
  delete dijkstra_;
  delete weightMap_;
}

void mixedCycle::computeMaximalMixedCycle(unsigned int numberOfStrings)
{
  StaticDigraph::Node start, end;
  unsigned int added(0), length;
  for(unsigned int i=0; i<numberOfStrings; ++i)
  {
    length = inputSizes_[i];
    for(unsigned int j=0; j<length-1; ++j)
    {
      start = graph_.node(added + j+1);
      end = graph_.node(added + j);
      dijkstra_->run(start,end);
      valueOfMixedCycle_ = dijkstra_->dist(end);
      if( valueOfMixedCycle_ < 1 ) 
      {
	endNodeString_ = i; 
	endNode_ = added + j;
	return;
      }
    }
    added +=length;
  }
}

void mixedCycle::getEdgesInViolatedMixedCycle(graphMsa::edges& edgesInMixedCycle)
{
  assert( edgesInMixedCycle.empty() );
  assert( endNode_ > -1 );

  unsigned int numberOfStrings(inputSizes_.size()), 
    currentString(endNodeString_),
    nodeI(0), nodeJ(0), strI, strJ,
    added(0);
  for(unsigned int i=0; i<currentString; ++i)
    added += inputSizes_[i];

  SimplePath<StaticDigraph> path = dijkstra_->path(graph_.node(endNode_));
  for(SimplePath<StaticDigraph>::ArcIt it(path); it!=INVALID; ++it)
  {
    strI = currentString % numberOfStrings;

    strJ = (currentString+1) % numberOfStrings;

    nodeI = arcs_[ graph_.index(it) ].first - added;

    added += inputSizes_[currentString % numberOfStrings];

    if( strJ == 0 ) 
      added = 0; // we start to count the nodes from the top left again
    
    nodeJ = arcs_[ graph_.index(it) ].second - added;

    if( strJ == 0 ) // or strI == numberOfStrings-1
    {
      swap(nodeI,nodeJ);
      swap(strI,strJ);
    }

    edgesInMixedCycle.push_back( edge(nodeI,nodeJ,strI,strJ) ); 
    ++currentString;
  }
}

bool mixedCycle::violatedMixedCycleExist()
{
  return valueOfMixedCycle_ > compareToOne;
}

/* builds digraph as follows: (n stands for node)
   n_0  n_1 ...  n_|str1|-1        <-- nodes corresponds to first string

   n_|str1| n_|str1|+1  ...  n_|str1|+|str2|-1    <-- nodes corresponds to second string 

         .                         .
	 .            ...          .
	 .                         .

   n_|str1|+...+|strk-1| ...  n_|str1|+...+|strK|-1    <-- nodes corresponds to last string
   
   
   arcs: each node of level i is connected to every nodes of level i+1 (with i=k and i+1=1)
 */
void mixedCycle::buildGraphInstance(unsigned int numberOfStrings, unsigned int numberOfNodes, const map<edge,double>& arcWeights)
{
  unsigned int iSize(0), jSize(0), added(0);
  vector<double> weights;
  for(unsigned int l=0; l<numberOfStrings; ++l)
  {
    iSize = inputSizes_[l]; // string l
    jSize = inputSizes_[(l+1) % numberOfStrings]; // string l+1 with l=0 for l+1==k
    for(unsigned int i=0; i<iSize; ++i)
      for(unsigned int j=0; j<jSize; ++j)
      {
	/* weights[i] is the corresponding weight of arcs[i] */
	arcs_.push_back( make_pair(added+i,(added+iSize+j) % numberOfNodes ) ); // add arc
	assert( arcWeights.count( edge(i,j,l,(l+1) % numberOfStrings) ) > 0 );
	weights.push_back( arcWeights.find(edge(i,j,l,(l+1) % numberOfStrings))->second); // add weight
      }
    added += iSize;
  }

  graph_.build(numberOfNodes, arcs_.begin(), arcs_.end() );
  weightMap_ = new StaticDigraph::ArcMap<double>(graph_,numberOfNodes);
  
  for(StaticDigraph::ArcIt it(graph_); it!=INVALID; ++it)
  {
    assert( weights[ graph_.index(it) ] >= 0.0 );
    (*weightMap_)[it] = weights[ graph_.index(it) ];
  }
  
  dijkstra_ = new Dijkstra<StaticDigraph,StaticDigraph::ArcMap<double> >(graph_,*weightMap_);

}


