/**
 * @file edgeClique.cpp
 * @brief edge clique
 * @author Sebastian Schenker
 **/

#include "edgeClique.h"
#include "graphMsa.h"
#include "graphVariable.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <lemon/static_graph.h>
#include <lemon/concepts/maps.h>
#include <lemon/dijkstra.h>
#include <map>
#include <vector>
#include <utility>

using namespace std;
using namespace lemon;

const double edgeClique::epsilon = 1.0e-09;

int edgeClique::getHorizontalLength() const 
{
  if( case_ji_ )
    return strI_length_;
  else 
    return strJ_length_;
}

int edgeClique::getVerticalLength() const
{
  if( case_ji_ )
    return strJ_length_;
  else
    return strI_length_;
}
  

edgeClique::edgeClique(int strI, int strJ, int strI_length, int strJ_length, const map<edge,double>& edgesInSol, bool case_ji) 
  : strI_(strI), strJ_(strJ), strI_length_(strI_length), strJ_length_(strJ_length), case_ji_(case_ji), edgeWeights_(0), dijkstra_(0)
{

  /* we add |strI_length| many source nodes to the pairgraph; that is necessary because 
     dijkstra works with weigths on arcs whereas we are given weights on the nodes and 
     have to adjust the pairgraph; each source node is connected to rightmost node corresponding 
     to lb-vertex in strI; nodeweights correspond to weights on incoming edges */
  
  int numberOfNodes = strI_length_*strJ_length_; // without sources!

  /* pair<int,int> needed by lemon staticDigraph constructor; The pairs
  must be in a non-decreasing order with respect to their first
  values! 
  Let the first string be strI = s1 s2 s3 and the second string strJ = t1 t2 t3 t4 with strI < strJ.
  The source nodes are numbered as follows: source1=numberOfNodes, source2=numberOfNodes+1, ...,
  The the pairgraph looks like this:
  
  n_11<- n_10<- n_9 <- n_8 <--source3  s3
   ^      ^      ^      ^
   |      |      |      |
  n_7 <- n_6 <- n_5 <- n_4 <--source2  s2
   ^      ^      ^      ^ 
   |	  |      |      |
  n_3 <- n_2 <- n_1 <- n_0 <--source1  s1
   t1     t2     t3     t4 
   
   The corresponding edges are (0,1),(0,4),(1,2),(1,5),...,(source1,0),(source2,4),(source3,8) 

   If we have case_ji == true, we exchange the positions of strI and strJ: string strI will be vertically aligned and string strJ horizontally aligned. This is reflected in functions getHorizontalLength() and getVerticalLength();
  */
  int vertLength(getVerticalLength()),
    horLength(getHorizontalLength());
  
  for(int i=0; i<vertLength; ++i) // i corresponds to vertical string
    for(int j=0; j<horLength; ++j) // j corresponds to horizontal string 
    {
      if( j < horLength-1 ) // add horizontal arcs; none for leftmost nodes
	edges_.push_back( make_pair(i*horLength+j, i*horLength+(j+1)) );

      if( i < vertLength-1 ) // add vertical arcs; none for topmost nodes
	edges_.push_back( make_pair(i*horLength+j, (i+1)*horLength+j) );
    }

  /* add arcs connecting source nodes; source nodes = numberOfNodes,...,numberOfNodes+strI_length-1  */
  for(int i=0; i<vertLength; ++i)
    edges_.push_back( make_pair(numberOfNodes+i, i*horLength) );

  /* build pairgraph */
  pairGraph_.build(numberOfNodes+vertLength, edges_.begin(), edges_.end());
    
  /* declare and construct edge weights needed for Dijkstra algo */
  edgeWeights_ = new StaticDigraph::ArcMap<double>(pairGraph_,numberOfNodes);

  /* we have a maximization problem (longest path problem; 
     in order to use dijkstra we set the edge weights to numberOfNodes - x^*_e */
  int nodeIndex; // node number in pairGraph
  map<edge,double>::const_iterator end(edgesInSol.end()),
    e;
  for(e=edgesInSol.begin(); e!=end; ++e)
  {
    /* nodes in pairgraph are indexed from bottom right; first character in string t corresponds to leftmost node in pairgraph */
    if( case_ji_ )
      nodeIndex = e->first.nodeJ()*strI_length_ + (strI_length_-1 - e->first.nodeI());
    else
      nodeIndex = e->first.nodeI()*strJ_length_ + (strJ_length_-1 - e->first.nodeJ());

    for(StaticDigraph::InArcIt it(pairGraph_,pairGraph_.node(nodeIndex)); it!=INVALID; ++it)
      (*edgeWeights_)[it] -= e->second;
  }

  dijkstra_ = new Dijkstra<StaticDigraph,StaticDigraph::ArcMap<double> >(pairGraph_,*edgeWeights_);

}

edgeClique::~edgeClique()
{
  delete dijkstra_;
  delete edgeWeights_;
}

double edgeClique::computeLongestPath(int lb, int le) const
{
  int numberOfNodes(strI_length_*strJ_length_),
    s(numberOfNodes + lb),
    t(le*getHorizontalLength() + getHorizontalLength()-1);
  
  StaticDigraph::Node endnode(pairGraph_.node(t));
  dijkstra_->run(pairGraph_.node(s),endnode);

  SimplePath<StaticDigraph> path = dijkstra_->path(endnode);
  // value of longest path
  return path.length() * numberOfNodes - dijkstra_->dist(endnode);
}

void edgeClique::getEdgesInClique(int le, graphMsa::edges& edgesInClique) const
{
  //assert( edgesInClique.empty() );
  int t(le*getHorizontalLength() + getHorizontalLength()-1);

  StaticDigraph::Node endnode(pairGraph_.node(t));
  SimplePath<StaticDigraph> path = dijkstra_->path(endnode);
  
  int nodeIndex, horizIndex, vertIndex;
  for(SimplePath<StaticDigraph>::ArcIt it(path); it!=INVALID; ++it) 
  {
    nodeIndex = edges_[pairGraph_.index(it)].second;
    horizIndex = getHorizontalLength() - 1 - ( nodeIndex % getHorizontalLength() );
    vertIndex = nodeIndex / getHorizontalLength();

    if( case_ji_ )
      swap(horizIndex,vertIndex);

    edgesInClique.push_back( edge(vertIndex, horizIndex, strI_, strJ_) );

  }  
}

void edgeClique::computeEdgeWeightsForLiftedMixedCycle(map<edge,double>& edgeWeightsForLiftedMixedCycle)
{
  int numberOfNodes(strI_length_ * strJ_length_);
  dijkstra_->run(pairGraph_.node(numberOfNodes));   // all pairs dijskstra 

  int horizIndex, vertIndex, 
    strI(strI_), 
    strJ(strJ_);

  double weightValue;

  if( case_ji_ )
    swap(strI,strJ);

  for(int i=0; i<numberOfNodes; ++i)
  {
    SimplePath<StaticDigraph> p = dijkstra_->path(pairGraph_.node(i));

    /* weightValue = 1 - \max_{K_E \in \cal{E}_ij(l<->|strI|,1<->m)} \sum_{e \in K_E} x^*_e */ 
    weightValue = 1.0 - p.length()*numberOfNodes + dijkstra_->dist(pairGraph_.node(i));
    if( abs(weightValue) < epsilon )
      weightValue = 0.0;

    horizIndex = getHorizontalLength() - 1 - ( i % getHorizontalLength() );
    vertIndex = i / getHorizontalLength();

    assert( edgeWeightsForLiftedMixedCycle.count(edge(vertIndex,horizIndex,strI,strJ)) == 0 );
    edgeWeightsForLiftedMixedCycle.insert( make_pair( edge(vertIndex,horizIndex,strI,strJ), weightValue) );
  }
}



