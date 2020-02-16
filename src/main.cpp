/**
 * @file msa.cpp
 * @brief multiple sequence alignment problem
 * @author Sebastian Schenker
 **/

#include "msa.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

int main(int argc, char* argv[]) {

  if( argc < 5 )
  {
    throw invalid_argument("you need to give at least four command line arguments: ./msa \"filename\" \"weight1\" \"weight2\" \"g or i\" \"(scoreMatrixFile)\" ");
  }
  
  // read input strings from given file 
  ifstream inputFile(argv[1]);
  if( !inputFile.good() )
    throw invalid_argument("input file could not be read!");

  vector<string> inputStrings;
  copy(istream_iterator<string>(inputFile), istream_iterator<string>(), back_inserter(inputStrings));
  if( inputStrings.size() < 3 )
    throw domain_error("input needs to be at least 3 strings");

  istringstream inputArg2(argv[2]),
    inputArg3(argv[3]);
    
  double lambda1,lambda2;
  inputArg2 >> lambda1;
  inputArg3 >> lambda2;

  msa* instance;
  if( argv[4][0] == 'g' )
    if( argc == 5 )
      instance = new msa(inputStrings,lambda1,lambda2,true,argv[5]); // consider gaps and scoreMatrixFile
    else
      instance = new msa(inputStrings,lambda1,lambda2,true); // consider gaps and matches/mismatches
  else if( argv[4][0] == 'i' )
    if( argc == 5 )
      instance = new msa(inputStrings,lambda1,lambda2,false,argv[5]); // consider indels and scoreMatrixFile
    else
      instance = new msa(inputStrings,lambda1,lambda2,false); // consider indels and matches/mismatches
  else
    throw domain_error("you need to specify whether gaps or indels are counted: g or i as input parameter");

  instance->computeAlignment();
  cout << "objective value = " << instance->getObjValue() << endl;
  instance->printSolution();
  instance->visualizeAlignment(inputStrings);

  delete instance;
  instance = 0;

  return 0;
}
