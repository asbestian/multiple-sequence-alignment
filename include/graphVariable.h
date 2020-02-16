/**
 * @file graphVariable.h
 * @brief represents variable (edge or arc) in msa graph model
 * @author Sebastian Schenker
 **/

#ifndef GRAPH_VARIABLE
#define GRAPH_VARIABLE

/// graphVariable

/* Represents a variable in the graph model of the msa problem. A
 variable is either an edge connecting two nodes or an (gap) arc. An
 edge is defined by nodeI, nodeJ, stringI and stringJ meaning that the
 edge connects nodeI of stringI with nodeJ from stringJ. An arc is
 defined by nodeI, nodeJ, stringI and stringJ meaning that the
 characters between nodeI and nodeJ of stringI are not aligned to
 any character of stringJ.  */

#include <cstddef>
#include <iosfwd>

class graphVariable {

 public: 
  graphVariable(std::size_t nodeI, std::size_t nodeJ, std::size_t stringI, std::size_t stringJ);

  std::size_t nodeI() const;
  std::size_t nodeJ() const;
  std::size_t strI() const;
  std::size_t strJ() const;
  std::ostream& print(std::ostream& out) const;
  
 private:
  std::size_t nodeI_,
    nodeJ_,
    strI_,
    strJ_;
  
};

typedef graphVariable edge;
typedef graphVariable arc;

std::ostream& operator<<(std::ostream& out, const graphVariable& var);

bool operator<(const graphVariable& lhs, const graphVariable& rhs);


#endif
