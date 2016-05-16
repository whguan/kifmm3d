/* Kernel Independent Fast Multipole Method
   Copyright (C) 2004 Lexing Ying, New York University

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.  */

/*!
 This file provides an example test for the sequential fmm3d.
 To generate the executable, type "make tt" and then run tt options_file
*/

#include "fmm3d.hpp"

using namespace std;

/*! parse the options file provided as input to tt
 * OptionsCreate takes a string and a character map and maps options file
 * names to their values to be parsed later
 */
int optionsCreate(const char* optionsFile, map<string,string>& options)
{
  options.clear();
  ifstream fileIn(optionsFile);
  if(fileIn.good()==false) {
	 cerr<<"wrong option file"<<endl;	 exit(1);
  }
  string name;  fileIn>>name;
  while(fileIn.good()) {
	 char content[100];	 fileIn.getline(content, 99);
	 options[name] = string(content);
	 fileIn>>name;
  }
  fileIn.close();
  return 0;
}

/*! Main function operates in several steps:
 * First, the options file as provided to optionsCreate is parsed.
 * From this, we get the number of sources and target points (numSrc, numTrg),
 * the kernel type (kt - see kernel3d.hpp and kernel3d.cpp), etc.  Random locations
 * are generated for source and target locations.  Random sources densities are allocated based
 * on the degrees of freedom of the specified kernel, and target value space is allocated as well.
 * Then, space for the fmm is allocated and variables are set.  Fmm is then this run and checked/clocked.
 */
int main(int argc, char** argv)
{
  srand48( (long)time(NULL) );  //srand48( 0 );
  iA(argc==2);
  map<string,string> optionsMap;
  optionsCreate(argv[1], optionsMap);
  
  ///1. allocate random data
  map<string,string>::iterator mapIdx;
  mapIdx = optionsMap.find("-numSrc"); assert(mapIdx!=optionsMap.end());
  /*! numSrc = Number of source points */
  int numSrc;  { istringstream ss((*mapIdx).second);  ss >> numSrc; }
  mapIdx = optionsMap.find("-numTrg"); assert(mapIdx!=optionsMap.end());
  /*! numTrg = Number of target points */
  int numTrg;  { istringstream ss((*mapIdx).second);  ss >> numTrg; }
  mapIdx = optionsMap.find("-kt"); assert(mapIdx!=optionsMap.end());
  /* ht = Kernel Type.  See kernel3d.hpp */
  int kernelType;  { istringstream ss((*mapIdx).second);  ss>>kernelType; }
    
  vector<double> tmp(2);	 tmp[0] = 1;	 tmp[1] = 0.25; //coefs in the kernel, work for all examples
  /*! Declare kernel as knl of type kt and coefficients temp */
  Kernel3d knl(kernelType, tmp);

  /*! Allocate random data for source positions (srcPos) and target positions (trgPos) */
  DblNumMat srcPos(3, numSrc);
  for(int i=0; i<numSrc; i++) {
	 srcPos(0,i) = (2.0*drand48()-1.0);
	 srcPos(1,i) = (2.0*drand48()-1.0);
	 srcPos(2,i) = (2.0*drand48()-1.0);
  }
  DblNumMat trgPos(3, numTrg);
  for(int i=0; i<numTrg; i++) {
	 trgPos(0,i) = (2.0*drand48()-1.0);
	 trgPos(1,i) = (2.0*drand48()-1.0);
	 trgPos(2,i) = (2.0*drand48()-1.0);
  }

  /*! srcDOF = source degree of freedom.  See kernel3d.hpp */
  int srcDOF = knl.srcDOF();
  /*! trgDOF = target degree of freedom.  See kernel3d.hpp */
  int trgDOF = knl.trgDOF();
  /*! srcDen = source density values of size sDOF*numSrc */
  DblNumVec srcDen(srcDOF * numSrc);
  for(int i=0; i<numSrc; i++) {
	 for(int d=0; d<srcDOF; d++)
		srcDen(d + i*srcDOF) = drand48(); //(2.0*drand48()-1.0);
  }
  /* trgVal = target values of size trgDOF *numTrg */
  DblNumVec trgVal(trgDOF * numTrg);
  
  ///2. allocate fmm 
  clock_t clockZero, clockOne;
  
  FMM3d* fmm = new FMM3d("fmm3d_");
  fmm->srcPos()=&srcPos;  fmm->srcNor()=&srcPos;
  fmm->trgPos()=&trgPos;
  fmm->center() = Point3(0,0,0); // CENTER OF THE TOPLEVEL BOX
  fmm->rootLevel() = 0;         // 2^(-rootlvl) is the RADIUS OF THE TOPLEVEL BOX
  fmm->knl() = knl;
  
  clockZero = clock();
  iC( fmm->setup(optionsMap) );
  clockOne = clock();  cout<<"fmm setup used "<<double(clockOne-clockZero)/CLOCKS_PER_SEC<<"secs "<<endl;
  
  ///3. run fmm
  for(int i=0; i<3; i++) {
	 clockZero = clock();
	 iC( fmm->evaluate(srcDen, trgVal) );
	 clockOne = clock();  cout<<"fmm eval used "<<double(clockOne-clockZero)/CLOCKS_PER_SEC<<"secs "<<endl;
  }
  
  ///4. check
  clockZero = clock();
  double relativeError;
  iC( fmm->check(srcDen, trgVal, 20, relativeError) );
  cout << "relative error: " << relativeError << endl;
  clockOne = clock();  cout<<"fmm check used "<<double(clockOne-clockZero)/CLOCKS_PER_SEC<<"sec "<<endl;
  
  delete fmm;
  
  return 0;
}
