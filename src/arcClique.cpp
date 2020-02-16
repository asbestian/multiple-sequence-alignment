/** 
 * @file arcClique.cpp
 * @brief Implementation of arc clique class
 * @author Sebastian Schenker
 **/

#include "arcClique.h"
#include "graphVariable.h"
#include <map>
#include <vector>

using namespace std;

arcClique::arcClique(int str_length, const map<arc,double>& arcsInSol)
  : omega_(str_length,0.0), pi_(str_length,0.0), sigma_(str_length,vector<double>(str_length,0.0)) 
{
  
  int start,end;
  /* help vector */
  vector< vector<double> > piHelp(str_length,pi_);
  
  /* compute values for piHelp, pi and omega: 
     pi[q] = \sum_{l=1}^q y^*_(v^i_l,v^i_q)
     omega[q] = \sum_{m=q}^|strI| y^*_(v^i_q,v^i_m): 
  */
  map<arc,double>::const_iterator itEnd(arcsInSol.end()),
    it;
  for(it=arcsInSol.begin(); it!=itEnd; ++it)
  {
    start = it->first.nodeI();
    end = it->first.nodeJ();
    /* helper vector: for arc with endnode p and startnode i in arcsInSolution piHelp contains 
       at position [p,i] solution value */
    piHelp[end][start] = it->second;
    /* correctness can be seen if above equations for pi and omega are expandeded */  
    pi_[end] += it->second;
    omega_[start] += it->second; 
  }
  
  /* compute values for sigma:
     sigma[lb,le] := \sum_{a \in A_ij(lb<->le)} y^*_a
   */
  sigma_[0][0] = omega_[0];
  /* sigma[0,i] = sigma[0,i-1] - piHelp[i-1,0] */
  for(int i=1; i<str_length; ++i)
    sigma_[0][i] = sigma_[0][i-1] - piHelp[i-1][0];

  /* sigma[p,q] = sigma[p-1,q] + piHelp[q,p] + piHelp[q+1,p] +...+ piHelp[|strI],p] */
  for(int i=1; i<str_length; ++i)
    for(int j=i; j<str_length; ++j)
    {
      sigma_[i][j] = sigma_[i-1][j];
      for(int k=j; k<str_length; ++k)
	sigma_[i][j] += piHelp[k][i];
    }
}

double arcClique::valueOfArcsSpanningEitherOr(int m) const 
{
  return sigma_[m][m+1] + pi_[m] + omega_[m+1];
}






