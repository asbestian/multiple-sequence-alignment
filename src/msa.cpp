/**
 * @file msa.cpp
 * @brief multiple sequence alignment problem
 * @author Sebastian Schenker
 **/

#include "msa.h"
#include "generalTrans.h"
#include "graphMsa.h"
#include "graphVariable.h"
#include "maxClique.h"
#include "mixedCycle.h"
#include <algorithm>
#include <cstddef>
#include <functional>
#include <ilcplex/ilocplex.h>
#include <ilconcert/iloexpression.h>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace std;

const double msa::epsilon = 1.0e-09;

void msa::addAssignmentConsToModel(const vector<size_t>& edgeInds, const vector<size_t>& arcInds)
{
  IloExpr ineq(env_);
  vector<size_t>::const_iterator end(edgeInds.end()),it;
  for(it = edgeInds.begin(); it!=end; ++it)
    ineq += edgeVars_[*it];
  
  end = arcInds.end();
  for(it = arcInds.begin(); it!= end; ++it)
    ineq += arcVars_[*it];

  basicModel_.add( ineq == 1.0 );
  ineq.end();
}


void msa::buildScoreMatrix(const char* scoreMatrixFile, map<pair<string,string>,int>& scoreMatrix)
{
  ifstream file(scoreMatrixFile);
  if( !file.good())
    throw invalid_argument("score matrix file could not good!");
  
  bool lineNotProcessed(true),
    notFinished;
  string letter,line,substring;
  vector<string> consideredLetters;
  size_t counter, pos, end,
    lineEnd;
  int score;
  while( !file.eof() )
  {
    getline(file,line);
    if( line.find("#",0) != string::npos || line.empty() )   // skip comments and empty lines
      continue;  

    pos = 0;
    lineEnd = line.size();
    while( line[pos] == ' ' ) // skip whitespace
      ++pos;
    end = line.find(' ',pos);
    substring = line.substr(pos,end-pos);

    if( lineNotProcessed )
      consideredLetters.push_back( substring );      
    else
      letter = substring;

    notFinished = true;
    counter = 0;
    while( notFinished )
    {
      pos = end + 1;
      while( line[pos] == ' ' ) // skip whitespace
	++pos;
      end = line.find(' ',pos);
      
      if( pos == lineEnd )
	break;
      else if( end == string::npos)
      {
	end = line.size();
	notFinished = false;
      }

      substring = line.substr(pos,end-pos);

      if( lineNotProcessed )
	consideredLetters.push_back( substring );
      else 
      {
	stringstream s(substring);
	s >> score;
	scoreMatrix[make_pair(letter,consideredLetters[counter])] = score;
	++counter;
      }
    }
    lineNotProcessed = false;
  }
}

msa::msa(const vector<string>& inputStrings, double lambda1, double lambda2, bool countGaps, char* scoreMatrixFile)
  : countGaps_(countGaps),
    objValue_(0.0), 
    graph_(inputStrings), 
    basicModel_(env_), 
    edgeVars_(env_,graph_.numberOfEdges(),0.0,1.0,ILOFLOAT), 
    arcVars_(env_,graph_.numberOfArcs(),0.0,1.0,ILOFLOAT), 
    edgeWeights_(env_,graph_.numberOfEdges()),
    arcWeights_(env_,graph_.numberOfArcs()), 
    inputSizes_(inputStrings.size())
{
  transform(inputStrings.begin(),inputStrings.end(),inputSizes_.begin(),mem_fun_ref(&string::size));

  unsigned int numberOfEdges(graph_.numberOfEdges()),
    numberOfArcs(graph_.numberOfArcs());

  graphMsa::edges graphEdges(graph_.getEdges());
  graphMsa::arcs graphArcs(graph_.getArcs());

  map< pair<string,string>, int> scoreMatrix;
  if( scoreMatrixFile )
    buildScoreMatrix(scoreMatrixFile,scoreMatrix);
  
  for(unsigned int i=0; i<numberOfEdges; ++i)
  {
    edgeWeights_[i] = edgeWeight(inputStrings,graphEdges[i],scoreMatrix);
    edgesToInds_.insert( make_pair(graphEdges[i],i) );
  }
  for(unsigned int i=0; i<numberOfArcs; ++i)
  {
    arcWeights_[i] = arcWeight(graphArcs[i]);
    arcsToInds_.insert( make_pair(graphArcs[i],i) );
  }

  //gaps_ = IloScalProd(arcWeights_,arcVars_);

  // add objective
  basicModel_.add(IloMaximize(env_,lambda1*IloScalProd(edgeWeights_,edgeVars_) - lambda2*IloScalProd(arcWeights_,arcVars_)));

  // add node assignment constraints
  addNodeAssignmentCons(inputSizes_);
  
}

void msa::printAlignment(const vector<string>& align)
{
  size_t width(align.size()),
    height(align[0].length());

  cout << "\n";

  for(size_t i=0; i<height; ++i)
  {
    for(size_t j=0; j<width; ++j)
      cout << align[j][i] << " ";
    cout << "\n";
  }
  cout << "\n";
}

void msa::buildAlignment(const vector<string>& inputStrings, vector<string>& alignment)
{
  size_t numberOfStrings(inputSizes_.size()),
    nodeJ,
    strJ;
    
  vector< vector<bool> > charProcessed; 
  for(size_t i=1; i<numberOfStrings; ++i)
    charProcessed.push_back( vector<bool>(inputSizes_[i],false) ); /* later use: true if index was aligned, false if it was not yet aligned */

  map<edge,double>::const_iterator itEnd(edgeSolutions_.end()),
    it;
  for(it=edgeSolutions_.begin(); it!=itEnd; ++it) 
  {
    if( it->first.strI() == 0 )
    {
      nodeJ = it->first.nodeJ();
      strJ = it->first.strJ();
      /* align nodeJ to nodeI in alignment */
      alignment[it->first.nodeI()][strJ] = inputStrings[strJ][nodeJ];
      /* set nodeJ index from stringJ to true */
      charProcessed[strJ-1][nodeJ] = true;
    }
  }

  size_t numberOfCharsBeforeNewColumn,
    end;
  vector<string>::iterator insertIter;
  for(size_t i=1; i<numberOfStrings; ++i) 
  {    
    end = inputSizes_[i];
    for(size_t j=0; j<end; ++j)
    {
      if( charProcessed[i-1][j] )
	continue;
      else
      {
	string newColumn(numberOfStrings,'_');
	newColumn[i] = inputStrings[i][j];

	for(it=edgeSolutions_.begin(); it!=itEnd; ++it)
	  if( it->first.strI() == i && it->first.nodeI() == j )
	  {
	    strJ = it->first.strJ();
	    nodeJ = it->first.nodeJ();
	    newColumn[strJ] = inputStrings[strJ][nodeJ];
	    charProcessed[strJ-1][nodeJ] = true;
	  }

	numberOfCharsBeforeNewColumn = 0;
	insertIter = alignment.begin();
	while( numberOfCharsBeforeNewColumn < j)
	{
	  if( (*insertIter)[i] != '_' )
	    ++numberOfCharsBeforeNewColumn;
	  ++insertIter;
	}
	alignment.insert(insertIter,newColumn);
      }
    }
  }
}

void msa::visualizeAlignment(const vector<string>& inputStrings)
{
  size_t numberOfStrings(inputStrings.size()),
    str0Length(inputStrings[0].length());
  
  vector<string> alignment(str0Length,string(numberOfStrings,'_'));
  for(size_t i=0; i<str0Length; ++i)
    alignment[i][0] = inputStrings[0][i];

  buildAlignment(inputStrings,alignment);
  printAlignment(alignment);
}

double msa::getObjValue() const
{
  return objValue_;
}

void msa::printSolution()  
{
  cout << "EDGES IN SOLUTION:\n";
  map<graphVariable,double>::iterator end(edgeSolutions_.end()),
    it;
  for(it=edgeSolutions_.begin(); it!=end; ++it)
    cout << it->first << "\n";

  cout << "ARCS IN SOLUTION:\n";
  end = arcSolutions_.end();
  for(it=arcSolutions_.begin(); it!=end; ++it)
    cout << it->first << "\n";
}

bool msa::computeAlignment()
{
  cout << "computeAlignment\n";
  IloModel model(env_);
  model.add(basicModel_);
  
  // add epsilon constraint w.r.t number of indels/gaps
  //  model.add( gaps_ <= epsilonConstraint );
  
  IloCplex cplex(model);
  cplex.setOut(env_.getNullStream());
  cplex.setWarning(env_.getNullStream());
  cplex.setParam(IloCplex::NumericalEmphasis, 1);
  cplex.setParam(IloCplex::EpRHS, 1e-9);
  cplex.setParam(IloCplex::TiLim,6000);
  cplex.exportModel("test.lp");

  map<edge,double> edgeWeightsForLiftedMixedCycle;
  bool separationNotFinished(true),
    furtherCheckNeeded(true);
  cout << "separation in progress...\n";
  while( separationNotFinished )
  {
    cplex.solve();
    if(cplex.getStatus() != IloAlgorithm::Optimal && cplex.getStatus() != IloAlgorithm::Feasible)
      return false;

    /* determine positive solutions in current linear program */
    determineObjectiveAndSolutions(cplex);

    maxClique clique(inputSizes_,edgeSolutions_,arcSolutions_,edgeWeightsForLiftedMixedCycle);
    if( clique.violatedMaxCliqueExist() ) 
    {
      //cout << "violated max clique exist" << endl;
      separationNotFinished = true;
      addMaxCliqueConsToLP(clique.getEdgesInViolatedClique(),clique.getArcsInViolatedClique(),model);
      continue;
    }

    mixedCycle cycle(inputSizes_,edgeWeightsForLiftedMixedCycle);
    if( cycle.violatedMixedCycleExist() )
    {
      //cout << "violated mixed cycle exist" << endl;
      separationNotFinished = true;
      graphMsa::edges edgesInVioMixCycle;
      cycle.getEdgesInViolatedMixedCycle(edgesInVioMixCycle);
      addMixedCycleConsToLP(edgesInVioMixCycle,model,inputSizes_.size()-1.0);
      continue;
    }

    generalTrans trans(inputSizes_,edgeSolutions_);
    if( trans.violatedTransConsExist() )
    {
      //cout << "violated general trans exist" << endl;
      separationNotFinished = true;
      graphMsa::edges posEdgesInVioTransCons,
	negEdgesInVioTransCons;
      trans.getEdgesInViolatedTransCons(inputSizes_,posEdgesInVioTransCons,negEdgesInVioTransCons);
      addGeneralTransConsToLP(posEdgesInVioTransCons,negEdgesInVioTransCons,model);
      continue; 
    }
     
    if( furtherCheckNeeded && solutionIsFractional() )
    {
      furtherCheckNeeded = false;
      cout << "solution is fractional...   ...consider now integer program\n";
      model.add(IloConversion(env_,edgeVars_,ILOBOOL));
      model.add(IloConversion(env_,arcVars_,ILOBOOL));
      separationNotFinished = true;
      continue;
    }

    separationNotFinished = false;
    cout << "separation finished\n";
  }

  cplex.end();
  model.end();

  return true;
}

msa::~msa() 
{
  arcWeights_.end();
  edgeWeights_.end();
  arcVars_.end();
  edgeVars_.end();
  basicModel_.end();
  env_.end();
}

bool msa::solutionIsFractional()
{
  map<graphVariable,double>::const_iterator end(edgeSolutions_.end()),
    it;
  for(it=edgeSolutions_.begin(); it!=end; ++it)
    if( it->second > epsilon && it->second < 1.0-epsilon )
      return true;

  end = arcSolutions_.end();
  for(it=arcSolutions_.begin(); it!=end; ++it)
    if( it->second > epsilon && it->second < 1.0-epsilon )
      return true;

  return false;
							     
}

void msa::addGeneralTransConsToLP(const graphMsa::edges& posEdges, const graphMsa::edges& negEdges, IloModel& model)
{
  IloExpr cons(env_);
  graphMsa::edges::const_iterator end(posEdges.end()),
    it;
  for(it=posEdges.begin(); it!=end; ++it)
    cons += edgeVars_[edgesToInds_.find(*it)->second];

  end = negEdges.end();
  for(it=negEdges.begin(); it!=end; ++it)
    cons -= edgeVars_[edgesToInds_.find(*it)->second];
  
  model.add(cons <= 1.0);
  cons.end();
  
}

void msa::addMixedCycleConsToLP(const graphMsa::edges& edgesInViolatedCons, IloModel& model,  double rhs)
{
  IloExpr cons(env_);
  graphMsa::edges::const_iterator end(edgesInViolatedCons.end()),
    it;
  for(it=edgesInViolatedCons.begin(); it!=end; ++it)
    cons += edgeVars_[edgesToInds_.find(*it)->second];
  
  model.add(cons <= rhs);
  cons.end();
}


void msa::addMaxCliqueConsToLP(const graphMsa::edges& edgesInViolatedCons, const graphMsa::arcs& arcsInViolatedCons, IloModel& model)
{
  IloExpr cons(env_);

  graphMsa::edges::const_iterator eEnd(edgesInViolatedCons.end()),
    eIt;
  for(eIt=edgesInViolatedCons.begin(); eIt!=eEnd; ++eIt)
    cons += edgeVars_[edgesToInds_.find(*eIt)->second];
   
  graphMsa::arcs::const_iterator aEnd(arcsInViolatedCons.end()),
    aIt;
  for(aIt=arcsInViolatedCons.begin(); aIt!=aEnd; ++aIt)
    cons += arcVars_[arcsToInds_.find(*aIt)->second];

  model.add(cons <= 1.0);
  cons.end();
}

void msa::determineObjectiveAndSolutions(const IloCplex& cp)
{
  objValue_ = cp.getObjValue();
  
  edgeSolutions_.clear();
  arcSolutions_.clear();

  size_t numberOfEdges(graph_.numberOfEdges());
  size_t numberOfArcs(graph_.numberOfArcs());
  
  graphMsa::edges graphEdges(graph_.getEdges());
  graphMsa::arcs graphArcs(graph_.getArcs());

  double value;
  for(size_t i=0; i<numberOfEdges; ++i)
    if( (value=cp.getValue(edgeVars_[i])) > epsilon )
    {
      edgeSolutions_.insert( make_pair(graphEdges[i],value) );
    }
  
  for(size_t i=0; i<numberOfArcs; ++i)
    if( (value=cp.getValue(arcVars_[i])) > epsilon )
    {
      arcSolutions_.insert( make_pair(graphArcs[i],value) );
    }
}

/* add constraints of form \sum_{m=1}^|s_j| x_{v^i_l,v^j_m} + 
   + \sum_{a \in A_ij(l<->l) y_a = 1, i,j=1,...,k;i!=j;l=1,...,|s^i| */
void msa::addNodeAssignmentCons(const vector<size_t>& inputSizes)
{
  vector<size_t> edgeIndicesInCons;
  vector<size_t> arcIndicesInCons;
  size_t numberOfStrings(inputSizes.size());
  for(size_t i=0; i<numberOfStrings-1; ++i)
    for(size_t j=i+1; j<numberOfStrings; ++j)
    {
      pair<size_t,size_t> edgeInds = graph_.indsForEdgesBetween(i,j);
      
      /* consider assignment w.r.t. string I (case_ij == true) */
      for(size_t k=0; k<inputSizes[i]; ++k)
      {
	pair<size_t,size_t> arcInds = graph_.indsForArcsBetween(i,j);
	getIndsForNodeAssignmentCons(k,edgeInds.first,edgeInds.second,arcInds.first,arcInds.second,graph_.getEdges(),graph_.getArcs(),edgeIndicesInCons,arcIndicesInCons,true);
	
	addAssignmentConsToModel(edgeIndicesInCons,arcIndicesInCons);
	
	edgeIndicesInCons.clear();
	arcIndicesInCons.clear();
      }

      /* consider assignment w.r.t. string J (case_ij == false) */
      for(size_t k=0; k<inputSizes[j]; ++k)
      {
	pair<size_t,size_t> arcInds = graph_.indsForArcsBetween(j,i);
	getIndsForNodeAssignmentCons(k,edgeInds.first,edgeInds.second,arcInds.first,arcInds.second,graph_.getEdges(),graph_.getArcs(),edgeIndicesInCons,arcIndicesInCons,false);
	
	addAssignmentConsToModel(edgeIndicesInCons,arcIndicesInCons);
	
	edgeIndicesInCons.clear();
	arcIndicesInCons.clear();
      }
    }
}

/* needed edges and arcs for node assignment constraint for node I between strI and strJ */
void msa::getIndsForNodeAssignmentCons(size_t consideredNode, size_t edgeBegin, size_t edgeEnd, size_t arcBegin, size_t arcEnd, const graphMsa::edges& graphEdges, const graphMsa::arcs& graphArcs, vector<size_t>& edgeIndsInCons, vector<size_t>& arcIndsInCons, bool case_ij)
{
  assert( edgeIndsInCons.empty() );
  assert( arcIndsInCons.empty() );
  
  size_t i, node;
  for(i=edgeBegin; i<edgeEnd; ++i)
  {
    if( case_ij )
      node = graphEdges[i].nodeI();
    else
      node = graphEdges[i].nodeJ();

    if( node == consideredNode )
      edgeIndsInCons.push_back(i);
  }

  for(i=arcBegin; i<arcEnd; ++i)
  {
    if( graphArcs[i].nodeI() <= consideredNode && graphArcs[i].nodeJ() >= consideredNode )
      arcIndsInCons.push_back(i);
  }

}

int msa::arcWeight(const arc& a, int valueIndel)
{
  return countGaps_ ? valueIndel : valueIndel*(a.nodeJ()-a.nodeI()+1);
}

int msa::edgeWeight(const vector<string>& input, edge e, map<pair<string,string>,int>& scoreMatrix) 
{
  if( scoreMatrix.empty() )
    return input[e.strI()][e.nodeI()] == input[e.strJ()][e.nodeJ()] ? 1 : -1;
  else
  {    
    string s1 = input[e.strI()].substr(e.nodeI(),1),
      s2 = input[e.strJ()].substr(e.nodeJ(),1);
    return scoreMatrix.at( make_pair(s1,s2) );
  }
}

