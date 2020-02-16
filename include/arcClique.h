/**
 * @file arcClique.h
 * @brief arc clique class
 * @author Sebastian Schenker
 **/

#ifndef ARC_CLIQUE
#define ARC_CLIQUE

#include "graphVariable.h"
#include <map>
#include <vector>

class arcClique {

 public:
  
  arcClique(int str_length, const std::map<arc,double>& arcsInSolution);

  double valueOfArcsSpanningEitherOr(int m) const;
  
 private:
  
  std::vector<double> omega_; // omega[q] = value of arcsInSolution with q as startnode
  std::vector<double> pi_; // pi[q] = value of arcsInSolution with q as endnode
  std::vector< std::vector<double> > sigma_;

};

#endif
