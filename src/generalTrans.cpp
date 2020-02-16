/**
 * @file genralTrans.cpp
 * @brief Implementation of the general transitivity class
 * @author Sebastian Schenker
 **/

#include "generalTrans.h"
#include "graphMsa.h"
#include "graphVariable.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <lemon/static_graph.h>
#include <lemon/preflow.h>
#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace std;
using namespace lemon;

const double generalTrans::compareToOne = 1.0 + 1.0e-09;

bool generalTrans::violatedTransConsExist()
{
  return minCutValue_ > compareToOne;
}

generalTrans::generalTrans(const vector<size_t>& inputSizes, const map<edge,double>& edgesInSol)
  : minCutValue_(0), sourceNode_(0), sourceNodeString_(0), str2_(0), str3_(0), capacities_(0), flow_(0)
{
  unsigned int numberOfStrings(inputSizes.size());
 
  for(unsigned int i=0; i<numberOfStrings-2; ++i)
    for(unsigned int j=i+1; j<numberOfStrings-1; ++j)
      for(unsigned int k=j+1; k<numberOfStrings; ++k)
      {
	if( computeMaxFlow(i,j,k,inputSizes,edgesInSol) || computeMaxFlow(i,j,k,inputSizes,edgesInSol,true,false) || computeMaxFlow(i,j,k,inputSizes,edgesInSol,false,true) )
	  return;
      }
}

generalTrans::~generalTrans()
{
  delete flow_;
  delete capacities_;
}

void generalTrans::getEdgesInViolatedTransCons(const vector<size_t>& inputSizes, graphMsa::edges& positiveEdgesInViolatedTransCons, graphMsa::edges& negativeEdgesInViolatedTransCons )
{
  assert( positiveEdgesInViolatedTransCons.empty() );
  assert( negativeEdgesInViolatedTransCons.empty() );

  unsigned int str2_length(inputSizes[str2_]);
  
  vector<unsigned int>::const_iterator endS_2(S_2_.end()),
    itS_2;
  for(itS_2=S_2_.begin(); itS_2!=endS_2; ++itS_2)
  {
    if( sourceNodeString_ > str2_ )
      positiveEdgesInViolatedTransCons.push_back( edge(*itS_2-1,sourceNode_,str2_,sourceNodeString_) );
    else
      positiveEdgesInViolatedTransCons.push_back( edge(sourceNode_,*itS_2-1,sourceNodeString_,str2_) );
  }
  
  vector<unsigned int>::const_iterator endS_3(S_3_.end()),
    itS_3;
  for(itS_3=S_3_.begin(); itS_3!=endS_3; ++itS_3)
  {
    if( sourceNodeString_ > str3_ )
      positiveEdgesInViolatedTransCons.push_back( edge(*itS_3 - (str2_length+1),sourceNode_,str3_,sourceNodeString_) );
    else
      positiveEdgesInViolatedTransCons.push_back( edge(sourceNode_,*itS_3 - (str2_length+1),sourceNodeString_,str3_) );
  }

  for(itS_2=S_2_.begin(); itS_2!=endS_2; ++itS_2)
    for(itS_3=S_3_.begin(); itS_3!=endS_3; ++itS_3)
      negativeEdgesInViolatedTransCons.push_back( edge(*itS_2-1,*itS_3-(str2_length+1),str2_,str3_) );
}

void generalTrans::buildGraph(const vector<size_t>& inputSizes)  
{
  unsigned int str2_length(inputSizes[str2_]),
    str3_length(inputSizes[str3_]);

  int numberOfNodes(2+str2_length+str3_length); 
  /* the built graph looks as follows:
                   n_1          n_|str2|+1
               /   n_2          n_|str2|+2      \
 (source)   n_0-->  .              .          -->n_|str2|+|str3|+1  (sink)
               \    .              .            /
                    .              .
	         n_|str2\     n_|str2|+|str3|
	   
     n_0 is connected to nodes (n_1,...,n_|str2|) and nodes (n_|str2|+1,...,n_|str2|+|str3|) connected 
     to the sink. Disregarding n_0 and the sink, the subgraph corresponding to str2 and str3 is bipartite!
  */
  arcs_.clear();
  for(unsigned int i=1; i<=str2_length; ++i)
    arcs_.push_back( make_pair(0,i) ); // add arcs from source 
  
  for(unsigned int i=1; i<=str2_length; ++i)
    for(unsigned int j=str2_length+1; j<=str2_length+str3_length; ++j)
      arcs_.push_back( make_pair(i,j) );   // add arcs between str2 and str3

  for(unsigned int i=str2_length+1; i<=str2_length+str3_length; ++i)
    arcs_.push_back( make_pair(i,numberOfNodes-1) ); // add arcs between nodes corresponding to str2 and the sink
 
  /* construct graph */
  graph_.clear();
  graph_.build(numberOfNodes, arcs_.begin(), arcs_.end());
}


/* add capacities corresponding to source node; for enumeration of graph see buildGraph() */
void generalTrans::addSourceNodeCapacities(const map<edge,double>& edgesInSolution) 
{
  unsigned int nodeI, nodeJ, strI, strJ;
  map<edge,double>::const_iterator end(edgesInSolution.end()),
    it;
  for(it=edgesInSolution.begin(); it!=end; ++it)
  {
    nodeI = it->first.nodeI();
    nodeJ = it->first.nodeJ();
    strI = it->first.strI();
    strJ = it->first.strJ();

    if( sourceNodeString_ > str2_ ) // due to datastructure: strI() < strJ() for graphVariable
    {
      swap(nodeI,nodeJ);
      swap(strI,strJ);
    }

    if( strI == sourceNodeString_ && strJ == str2_ && nodeI == sourceNode_ ) // edge contributes to capacites
      (*capacities_)[graph_.arc(nodeJ)] = it->second;
  }
}

/* add capacities corresponding to sink node; for enumeration of graph see buildGraph() */
void generalTrans::addSinkNodeCapacities(const vector<size_t>& inputSizes, const map<edge,double>& edgesInSolution)
{
  unsigned int addTo(inputSizes[str2_]*(1+inputSizes[str3_])),
    nodeI, nodeJ, strI, strJ;
  map<edge,double>::const_iterator end(edgesInSolution.end()),
    it;
  for(it=edgesInSolution.begin(); it!=end; ++it)
  {
    nodeI = it->first.nodeI();
    nodeJ = it->first.nodeJ();
    strI = it->first.strI();
    strJ = it->first.strJ();
    
    if( sourceNodeString_ > str3_ ) // due to datastructure: strI() < strJ() for graphVariable
    {
      swap(nodeI,nodeJ);
      swap(strI,strJ);
    }
    
    if( strI == sourceNodeString_ && strJ == str3_ && nodeI == sourceNode_ ) // edge contributes to capacities
      (*capacities_)[graph_.arc(addTo+nodeJ)] = it->second;
  }
}

/* add capacities corresponding to str2 and str3; for enumeration of graph see buildGraph() */
void generalTrans::addRemainingArcCapacities(const vector<size_t>& inputSizes, const map<edge,double>& edgesInSolution)
{
  unsigned int str2_length(inputSizes[str2_]),
    str3_length(inputSizes[str3_]);
  map<edge,double>::const_iterator end(edgesInSolution.end()),
    it;
  for(it=edgesInSolution.begin(); it!=end; ++it)
    if( it->first.strI() == str2_ && it->first.strJ() == str3_ )
      (*capacities_)[graph_.arc(str2_length + it->first.nodeI()*str3_length + it->first.nodeJ())] = it->second;
}

void generalTrans::addCapacities(const vector<size_t>& inputSizes, const map<edge,double>& edgesInSolution)
{
  delete capacities_;
  capacities_ = new StaticDigraph::ArcMap<double>(graph_,0.0);

  addSourceNodeCapacities(edgesInSolution);
  addSinkNodeCapacities(inputSizes,edgesInSolution);
  addRemainingArcCapacities(inputSizes,edgesInSolution);
}


/* case_strI=true means that strI is taken as the source node string in the transitivity  inequality */
bool generalTrans::computeMaxFlow(unsigned int str1, unsigned int str2, unsigned int str3, const vector<size_t>& inputSizes, const map<edge,double>& edgesInSolution, bool case_str2, bool case_str3)
{
  unsigned int str1_length(inputSizes[str1]),
    str2_length(inputSizes[str2]),
    str3_length(inputSizes[str3]),
    endNode;
  
  if( case_str2 )
  {
    sourceNodeString_ = str2;
    str2_ = str1;
    str3_ = str3;
    endNode = str1_length+str3_length+1;
    buildGraph(inputSizes);

  }
  else if( case_str3 )
  {
    sourceNodeString_ = str3;
    str2_ = str1;
    str3_ = str2;
    endNode = str1_length+str2_length+1;
    buildGraph(inputSizes);
  }
  else
  {
    sourceNodeString_ = str1;
    str2_ = str2;
    str3_ = str3;
    endNode = str2_length+str3_length+1;
    buildGraph(inputSizes);
  }

  StaticDigraph::Node s = graph_.node(0),
    t = graph_.node(endNode);

  sourceNode_ = 0;
  while( sourceNode_ < inputSizes[sourceNodeString_] )
  {
    addCapacities(inputSizes,edgesInSolution);

    /* create max flow instance */
    delete flow_;
    flow_ = new Preflow< StaticDigraph,StaticDigraph::ArcMap<double> >(graph_,*capacities_,s,t);
    flow_->runMinCut();

    /* compute min cut value */
    minCutValue_ = computeMinCutValue(inputSizes);
    if( minCutValue_ > compareToOne )
      return true;
    ++sourceNode_;
  }
  return false;
}

void generalTrans::getNodeCuts(const vector<size_t>& inputSizes)
{
  unsigned int str2_length(inputSizes[str2_]),
    str3_length(inputSizes[str3_]);

  S_2_.clear();
  S_3_.clear();

  for(unsigned int i=1; i<=str2_length; ++i)
    if( flow_->minCut(graph_.node(i)) )
      S_2_.push_back(i);
 
  for(unsigned int i=str2_length+1; i<=str2_length+str3_length; ++i)
    if( !flow_->minCut(graph_.node(i)) )
      S_3_.push_back(i);
}

double generalTrans::computeMinCutValue(const vector<size_t>& inputSizes)
{
  getNodeCuts(inputSizes);
  
  double profit(0.0),
    cost(0.0);
  unsigned int str2_length(inputSizes[str2_]),
    str3_length(inputSizes[str3_]);
  
  vector<unsigned int>::const_iterator endS_2(S_2_.end()),itS_2;
  for(itS_2=S_2_.begin(); itS_2!=endS_2; ++itS_2)
    profit += (*capacities_)[graph_.arc(*itS_2-1)];
  
  vector<unsigned int>::const_iterator endS_3(S_3_.end()),itS_3;
  for(itS_3=S_3_.begin(); itS_3!=endS_3; ++itS_3)
    profit += (*capacities_)[graph_.arc(str2_length*(1+str3_length) + (*itS_3 - str2_length-1) )];

  for(itS_2=S_2_.begin(); itS_2!=endS_2; ++itS_2)
    for(itS_3=S_3_.begin(); itS_3!=endS_3; ++itS_3)
      cost += (*capacities_)[graph_.arc(str2_length + (*itS_2-1)*str3_length + (*itS_3-str2_length-1))];
  
  assert( profit - cost >= 0.0 );
  return profit - cost;
}


