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
#include "dense3d.hpp"
#include "common/vecmatop.hpp"

/*! Constructor calls KnlMat3d constructor solely */
Dense3d::Dense3d(const string& p): KnlMat3d(p)
{
}

Dense3d::~Dense3d()
{
}

// ----------------------------------------------------------------------
int Dense3d::setup(map<string,string>&)
{
  iA(_srcPos!=NULL && _srcNor!=NULL && _trgPos!=NULL);
  iA((*_srcPos).m()==dim() && (*_trgPos).m()==dim());  //nothing to do
  return 0;
}

// ---------------------------------------------------------------------- 
int Dense3d::evaluate(const DblNumVec& srcDen, DblNumVec& trgVal) 
{
  //-----------------------------------
  iA(srcDen.m()==srcDOF()*(*_srcPos).n());  iA(trgVal.m()==trgDOF()*(*_trgPos).n());
  
  int dim  = this->dim();
  int srcDOF = this->srcDOF();
  int trgDOF = this->trgDOF();
  /* Number of sources */
  int numSrc = (*_srcPos).n();
  /* Number of targets */
  int numTrg = (*_trgPos).n();
  
  DblNumMat inter(trgDOF, numSrc*srcDOF);
  for(int i=0; i<numTrg; i++) {
	 DblNumMat onePosMat(dim, 1, false, (*_trgPos).clmdata(i));
	 DblNumVec oneValVec(trgDOF, false, trgVal.data()+trgDOF*i);
	 iC( _knl.kernel((*_srcPos), (*_srcNor), onePosMat, inter) );
	 iC( dgemv(1.0, inter, srcDen, 0.0, oneValVec) );
  }
  
  return 0;
}
