/**
 * @file msa.h
 * @brief multiple sequence alignment problem
 * @author Sebastian Schenker
 **/

#ifndef MSA_CLASS
#define MSA_CLASS

#include "graphMsa.h"
#include <cstddef>
#include <ilcplex/ilocplex.h>
#include <ilconcert/iloexpression.h>
#include <string>
#include <vector>

class msa {

 public:
  msa(const std::vector<std::string>& inputStrings, double lambda1, double lambda2, bool countGaps, char* scoreMatrix=nullptr); // countGaps=true -> count gaps, else count indels

  virtual ~msa();

  bool computeAlignment();

  double getObjValue() const;

  void visualizeAlignment(const std::vector<std::string>& inputStrings);
  
  void printSolution();

 private:
  static const double epsilon;
  
  bool countGaps_;
  double objValue_;
  graphMsa graph_;
  IloEnv env_;
  IloModel basicModel_;
  IloNumVarArray edgeVars_;
  IloNumVarArray arcVars_;
  IloNumArray edgeWeights_;
  IloNumArray arcWeights_;
  //IloExpr gaps_;

  std::map<edge,double> edgeSolutions_;
  std::map<arc,double> arcSolutions_;  
  std::map<edge,unsigned int> edgesToInds_;
  std::map<arc,unsigned int> arcsToInds_;
  std::vector<std::size_t> inputSizes_;

  /* each node in graph has to be aligned to edge or arc */
  void addNodeAssignmentCons(const std::vector<std::size_t>& inputSizes);

  void addAssignmentConsToModel(const std::vector<std::size_t>& edgeIndicies, const std::vector<std::size_t>& arcIndices);

  void getIndsForNodeAssignmentCons(std::size_t consideredNode, std::size_t edgeBegin, std::size_t edgeEnd, std::size_t arcBegin, std::size_t arcEnd, const graphMsa::edges& graphEdges, const graphMsa::arcs& graphArcs, std::vector<std::size_t>& edgeIndsInCons, std::vector<std::size_t>& arcIndsInCons, bool case_ij);

  void buildScoreMatrix(const char* scoreMatrixFile, std::map< std::pair<std::string,std::string>,int>& scoreMatrix);

  void determineObjectiveAndSolutions(const IloCplex& cp);
  
  bool solutionIsFractional();
  
  void buildAlignment(const std::vector<std::string>& inputStrings, std::vector<std::string>& alignment);
  
  void printAlignment(const std::vector<std::string>& alignment);
  
  void addMaxCliqueConsToLP(const graphMsa::edges& edgesInViolatedMaxClique, const graphMsa::arcs& arcsInViolatedMaxClique, IloModel& model);

  void addMixedCycleConsToLP(const graphMsa::edges& edgesInViolatedMixedCycle, IloModel& model, double rhs);
 
  void addGeneralTransConsToLP(const graphMsa::edges& posEdgesInViolatedTransCons, const graphMsa::edges& negEdgesInViolatedTransCons, IloModel& model);

  int edgeWeight(const std::vector<std::string>& inputStrings, edge e, std::map< std::pair<std::string,std::string>,int>& scoreMatrix);
  int arcWeight(const arc& a, int valueIndel=1); 
  


};

#endif
