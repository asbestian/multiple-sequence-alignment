/**
 * @file graphVariable.cpp
 * @brief represents variable (edge or arc) in msa graph model
 * @author Sebastian Schenker
 **/

#include "graphVariable.h"
#include <cstddef>
#include <iostream>

using namespace std;

bool operator<(const graphVariable& lhs, const graphVariable& rhs)
{
  if( lhs.strI() < rhs.strI() )
    return true;
  else if( lhs.strI() > rhs.strI() )
    return false;
  else // lhs.strI() == rhs.strI()
  {
    if( lhs.strJ() < rhs.strJ() )
      return true;
    else if( lhs.strJ() > rhs.strJ() )
      return false;
    else // lhs.strI()==rhs.strI() && lhs.strJ()==rhs.strJ() 
    {
      if( lhs.nodeI() < rhs.nodeI() )
	return true;
      else if( lhs.nodeI() > rhs.nodeI() )
	return false;
      else // lhs.strI()==rhs.strI() && lhs.strJ()==rhs.strJ() && lhs.nodeI()==rhs.nodeI()
      {
	if( lhs.nodeJ() < rhs.nodeJ() )
	  return true;
	else
	  return false;
      }
    }
  }
}

ostream& graphVariable::print(ostream& out) const {
  return out << "[" << nodeI_ << ", " << nodeJ_ << ", " << strI_ << ", " << strJ_ << "]";
}

ostream& operator<<(ostream& out, const graphVariable& var) {
  return var.print(out);
}

graphVariable::graphVariable(size_t nI, size_t nJ, size_t sI, size_t sJ) 
  : nodeI_(nI), nodeJ_(nJ), strI_(sI), strJ_(sJ) 
{}

size_t graphVariable::nodeI() const {
  return nodeI_;
}

size_t graphVariable::nodeJ() const {
  return nodeJ_;
}

size_t graphVariable::strI() const {
  return strI_;
}

size_t graphVariable::strJ() const {
  return strJ_;
}
