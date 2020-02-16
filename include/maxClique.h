/**
 * @file maxClique.h
 * @brief maximal clique constraint
 * @author Sebastian Schenker
 **/

#ifndef MAX_CLIQUE_CONSTRAINT
#define MAX_CLIQUE_CONSTRAINT

#include "graphMsa.h"
#include "graphVariable.h"
#include <cstddef>
#include <map>
#include <vector>

class maxClique {

 public: 
  
  maxClique(const std::vector<std::size_t>& inputSizes, const std::map<edge,double>& edgeSolutions, const std::map<arc,double>& arcSolutions, std::map<edge,double>& edgeWeightsForLiftedMixedCycle);

  bool violatedMaxCliqueExist();
  const graphMsa::edges& getEdgesInViolatedClique() const;
  const graphMsa::arcs& getArcsInViolatedClique() const;

 private: 
  static const double compareToOne;

  double valueOfMaxClique_;
  graphMsa::arcs arcsInMaxClique_;
  graphMsa::edges edgesInMaxClique_;

  bool violatedCliqueExist(unsigned int strI, unsigned int strJ, unsigned int strI_length, unsigned int strJ_length, const std::map<edge,double>& ijEdgesInSolution, const std::map<arc,double>& ijArcsInSolution, const std::map<arc,double>& jiArcsInSolution, bool addIJedgeWeights, bool addJIedgeWeights, std::map<edge,double>& edgeWeightsForLiftedMixedCycle);

  void copySolutions(unsigned int i, unsigned int j, const std::map<edge,double>& allSolutions, std::map<edge,double>& ijSols);
  void copySolutions(unsigned int i, unsigned int j, const std::map<arc,double>& allSolutions, std::map<arc,double>& ijSols, std::map<arc,double>& jiSols);

  double getArcValueSpanning_lb_le(unsigned int lb, unsigned int le, const std::map<arc,double>& arcsInSol);
  /* pushes arcs that span either m or m+1 into arcsInMaxClique */
  void createArcsSpanningEitherOr(unsigned int m, unsigned int strI, unsigned int strJ, unsigned int strI_length); 
  /* pushes arcs that span lb <-> le into arcsInMaxClique */
  void createArcsSpanning_lb_le(unsigned int lb,unsigned int le, unsigned int strI, unsigned int strJ, unsigned int strI_length);
};

#endif

