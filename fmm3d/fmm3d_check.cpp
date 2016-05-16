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
#include "fmm3d.hpp"
#include "common/vecmatop.hpp"

using std::cerr;
using std::endl;

/* ********************************************************************** */
int FMM3d::check(const DblNumVec& srcden, DblNumVec& trgVal, int numChk, double& rerr)
{
  iA(_trgPos->n()>0); //have things to check
  int dim = this->dim();
  int srcDOF = this->srcDOF();
  int trgDOF = this->trgDOF();
  int srcNum = _srcPos->n();
  int trgNum = _trgPos->n();   //int lclnum = this->lclnum();
  
  //1. select point to check
  vector<int> chkVec(numChk);
  for(int k=0; k<numChk; k++) {
	 chkVec[k] = int( floor(drand48()*trgNum) );	 iA(chkVec[k]>=0 && chkVec[k]<trgNum);
  }
  DblNumMat ChkPos(dim, numChk);
  DblNumVec chkVal(numChk*trgDOF);
  DblNumVec chkDen(numChk*trgDOF);
  for(int k=0; k<numChk; k++) {
	 for(int i=0; i<dim; i++)		ChkPos(i, k) = (*_trgPos)(i, chkVec[k]);
	 for(int i=0; i<trgDOF;i++)		chkVal(k*trgDOF+i) = trgVal(chkVec[k]*trgDOF+i);
  }
  
  //2. compute and gather
  DblNumMat inter(trgDOF, srcNum*srcDOF);
  for(int i=0; i<numChk; i++) {
	 DblNumMat oneChkPos(dim, 1, false, ChkPos.clmdata(i));
	 DblNumVec onechkDen(trgDOF, false, chkDen.data()+i*trgDOF);
	 iC( _knl.kernel(*_srcPos, *_srcNor, oneChkPos, inter) );
	 iC( dgemv(1.0, inter, srcden, 0.0, onechkDen) );
  }
  
  //3. distribute to individual
  for(int k=0; k<numChk; k++)
	 for(int i=0; i<trgDOF; i++)
		chkDen(k*trgDOF+i) -= chkVal(k*trgDOF+i);
  
  double vn = 0;  double en = 0;
  for(int k=0; k<numChk; k++)
	 for(int i=0; i<trgDOF; i++) {
		vn += chkVal(k*trgDOF+i) * chkVal(k*trgDOF+i);
		en += chkDen(k*trgDOF+i) * chkDen(k*trgDOF+i);
	 }
  vn = sqrt(vn);
  en = sqrt(en);
  
  //cerr<<"relative error: "<<en/vn<<endl;  cerr<<"error norm    : "<<en<<endl;  cerr<<"solution norm : "<<vn<<endl;
  rerr = en/vn; //relative error
  
  return 0;
}



