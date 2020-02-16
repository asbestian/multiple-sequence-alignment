/**
 * @file maxClique.cpp
 * @brief maximal clique constraint
 * @author Sebastian Schenker
 **/

#include "maxClique.h"
#include "arcClique.h"
#include "edgeClique.h"
#include "graphVariable.h"
#include <cassert>
#include <cstddef>
#include <map>
#include <vector>

using namespace std;

const double maxClique::compareToOne = 1.0 + 1.0e-09;

void maxClique::copySolutions(unsigned int i, unsigned int j, const map<edge,double>& allSols, map<edge,double>& ijSols)
{
  assert( ijSols.empty() );
  map<edge,double>::const_iterator end(allSols.end()),
    it;
  for(it=allSols.begin(); it!=end; ++it)
  {
    if( it->first.strI() == i && it->first.strJ() == j )
      ijSols.insert(*it);
  }
}

void maxClique::copySolutions(unsigned int i, unsigned int j, const map<arc,double>& allSols, map<arc,double>& ijSols, map<arc,double>& jiSols)
{
  map<graphVariable,double>::const_iterator end(allSols.end()),
    it;
  for(it=allSols.begin(); it!=end; ++it)
  {
    if( it->first.strI() == i && it->first.strJ() == j )
      ijSols.insert(*it);
    else if( it->first.strI() == j && it->first.strI() == i )
      jiSols.insert(*it);
  }
}

const graphMsa::edges& maxClique::getEdgesInViolatedClique() const 
{
  return edgesInMaxClique_;
}

const graphMsa::arcs& maxClique::getArcsInViolatedClique() const
{
  return arcsInMaxClique_;
}

maxClique::maxClique(const vector<size_t>& inputSizes, const map<edge,double>& edgeSols, const std::map<arc,double>& arcSols, map<edge,double>& edgeWeightsForLiftedMixedCycle)
  : valueOfMaxClique_(0.0)
{
  edgeWeightsForLiftedMixedCycle.clear();

  unsigned int numberOfStrings(inputSizes.size());
  map<edge,double> ijEdgesInSol,
    ijArcsInSol,
    jiArcsInSol;
  bool ret(false);
  for(unsigned int i=0; i<numberOfStrings-1; ++i)
    for(unsigned int j=i+1; j<numberOfStrings; ++j)
    {
      ijEdgesInSol.clear();
      ijArcsInSol.clear();
      jiArcsInSol.clear();

      copySolutions(i,j,edgeSols,ijEdgesInSol); 
      copySolutions(i,j,arcSols,ijArcsInSol,jiArcsInSol);
      
      if( i==0 && j==numberOfStrings-1 )
	ret = violatedCliqueExist(i,j,inputSizes[i],inputSizes[j],ijEdgesInSol,ijArcsInSol,jiArcsInSol,false,true,edgeWeightsForLiftedMixedCycle);
      else if( j==i+1 )
	ret = violatedCliqueExist(i,j,inputSizes[i],inputSizes[j],ijEdgesInSol,ijArcsInSol,jiArcsInSol,true,false,edgeWeightsForLiftedMixedCycle);
      else
	ret = violatedCliqueExist(i,j,inputSizes[i],inputSizes[j],ijEdgesInSol,ijArcsInSol,jiArcsInSol,false,false,edgeWeightsForLiftedMixedCycle);

      if( ret )
	return;	
    }
}


bool maxClique::violatedCliqueExist(unsigned int strI, unsigned int strJ, unsigned int strI_length, unsigned int strJ_length, const map<edge,double>& ijEdgesInSol, const map<arc,double>& ijArcsInSol, const map<arc,double>& jiArcsInSol, bool addIJedgeWeights, bool addJIedgeWeights, map<edge,double>& edgeWeightsForLiftedMixedCycle)
{
  unsigned int current(0),
    limit(strI_length-1);
  
  /* check for violated max clique of form K_E = {}, K_A = A_ij(m+1 <-> m) */
  arcClique ijArcClique(strI_length,ijArcsInSol);
  while( current < limit )
  {
    valueOfMaxClique_ = ijArcClique.valueOfArcsSpanningEitherOr(current);
    if( valueOfMaxClique_ > compareToOne ) // violated max clique found
    {
      createArcsSpanningEitherOr(current,strI,strJ,strI_length);
      return true;
    }
    ++current;
  }
  
  /* check for violated max clique inequality of form K_E = {}, K_A = A_ji(m+1 <-> m) */
  arcClique jiArcClique(strJ_length,jiArcsInSol);
  current = 0;
  limit = strJ_length-1;
  while( current < limit )
  {
    valueOfMaxClique_ = jiArcClique.valueOfArcsSpanningEitherOr(current);
    if( valueOfMaxClique_ > compareToOne ) //violated max clique found
    {
      createArcsSpanningEitherOr(current,strJ,strI,strJ_length);
      return true;
    }
    ++current;
  }

  /* no violated max clique of form K_E={}, K_A = A_ij(m+1 <-> m) or K_A = A_ji(m+1<->m)
     was found; */

  /* now checking for violated max clique of form K_E = \cal{E}_ij(lb <-> le, 1 <->|s_j|),
     K_A = A_ij(lb <-> le) */

  edgeClique ijEdgeClique(strI,strJ,strI_length,strJ_length,ijEdgesInSol);
  for(unsigned int i=0; i<strI_length; ++i)
    for(unsigned int j=i; j<strI_length; ++j)
    {
      valueOfMaxClique_ = ijEdgeClique.computeLongestPath(i,j) + getArcValueSpanning_lb_le(i,j,ijArcsInSol);
      if( valueOfMaxClique_ > compareToOne )
      {
	/* create edgesInMaxClique_ */
	ijEdgeClique.getEdgesInClique(j,edgesInMaxClique_);
	/* create arcsInMaxClique_ */
	createArcsSpanning_lb_le(i,j,strI,strJ,strI_length);
	return true;
      }
    }

  /* No violated max clique found so far! 
     We consider now the case E_ji and A_ji!
     check for max clique of form K_E = \cal{E}_ji{lb <-> le,1 <-> |s_i|},
     K_A = A_ji(lb <-> le) with 1 <= lb <= le <= |s_j| */
  edgeClique jiEdgeClique(strI,strJ,strI_length,strJ_length,ijEdgesInSol,true);
  for(unsigned int i=0; i<strJ_length; ++i)
    for(unsigned int j=i; j<strJ_length; ++j)
    {
      valueOfMaxClique_ = jiEdgeClique.computeLongestPath(i,j) + getArcValueSpanning_lb_le(i,j,jiArcsInSol);
      if( valueOfMaxClique_ > compareToOne )
      {
	/* call with order strI, strJ because we work only with edges connecting strI and strJ with strI < strJ;
	 create edgesInMaxClique_ */
	jiEdgeClique.getEdgesInClique(j,edgesInMaxClique_); 
	/* create arcsInMaxClique_ */
	createArcsSpanning_lb_le(i,j,strJ,strI,strJ_length);
	return true;
      }
    }

  /* possibility to add edgeWeights only if no violated maximal clique was found */
  if( addIJedgeWeights )
  {
    ijEdgeClique.computeEdgeWeightsForLiftedMixedCycle(edgeWeightsForLiftedMixedCycle);
  }
  if( addJIedgeWeights )
  {
    jiEdgeClique.computeEdgeWeightsForLiftedMixedCycle(edgeWeightsForLiftedMixedCycle);
  }

  // no violated clique found
  return false;
}


double maxClique::getArcValueSpanning_lb_le(unsigned int lb, unsigned int le, const map<arc,double>& arcsInSol)
{
  double arcValue(0.0);
  map<arc,double>::const_iterator end(arcsInSol.end()),
    it;
  for(it=arcsInSol.begin(); it!=end; ++it)
    if( it->first.nodeI() <= lb && it->first.nodeJ() >= le)
      arcValue += it->second;

  return arcValue;
}

void maxClique::createArcsSpanning_lb_le(unsigned int lb, unsigned int le, unsigned int strI, unsigned int strJ, unsigned int strI_length)
{
  for(unsigned int i=0; i<=lb; ++i)
    for(unsigned int j=le; j<strI_length; ++j)
      arcsInMaxClique_.push_back( arc(i,j,strI,strJ) );
}

void maxClique::createArcsSpanningEitherOr(unsigned int m, unsigned int strI, unsigned int strJ, unsigned int strI_length)
{
  /* create arcs having m as endpoint */
  for(unsigned int i=0; i<=m; ++i)
    arcsInMaxClique_.push_back( arc(i,m,strI,strJ) );

  /* create arcs having m+1 as startpoint */
  for(unsigned int i=m+1; i<strI_length; ++i)
    arcsInMaxClique_.push_back( arc(m+1,i,strI,strJ) );

  /* create arcs spanning m <-> m+1 */
  createArcsSpanning_lb_le(m,m+1,strI,strJ,strI_length);
}

bool maxClique::violatedMaxCliqueExist() 
{
  return valueOfMaxClique_ > compareToOne;
}

