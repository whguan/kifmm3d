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

using std::cerr;
using std::endl;

using std::istringstream;

// ---------------------------------------------------------------------- 
int FMM3d::setup(map<string,string>& opts)
{
  //-----------------------------------------------------
  map<string,string>::iterator mi;
  mi = opts.find("-" + prefix() + "np"); iA(mi!=opts.end());
  { istringstream ss((*mi).second); ss>>_np; }
    
  //-----------------------------------------------------
  iA(_srcPos!=NULL && _srcNor!=NULL && _trgPos!=NULL);
  iA((*_srcPos).m()==dim() && (*_trgPos).m()==dim());  //nothing to do
  //1. _let
  _let = new Let3d(prefix()+"let3d_");
  _let->srcPos()=_srcPos;  _let->trgPos()=_trgPos;  _let->center()=_center;  _let->rootLevel()=_rootLevel;
  iC( _let->setup(opts) );
  //2. decide _eq_mm and _mul_mm, and get matmgnt based on that
  switch(_knl.kernelType()) {
	 //laplace kernels
  case KNL_LAP_S_U: _knl_mm = Kernel3d(KNL_LAP_S_U, _knl.coefs()); break;
  case KNL_LAP_D_U: _knl_mm = Kernel3d(KNL_LAP_S_U, _knl.coefs()); break;
	 //stokes kernels
  case KNL_STK_S_U: _knl_mm = Kernel3d(KNL_STK_F_U, _knl.coefs()); break;
  case KNL_STK_S_P: _knl_mm = Kernel3d(KNL_LAP_S_U, vector<double>()); break;
  case KNL_STK_D_U: _knl_mm = Kernel3d(KNL_STK_F_U, _knl.coefs()); break;
  case KNL_STK_D_P: _knl_mm = Kernel3d(KNL_LAP_S_U, vector<double>()); break;
	 //navier kernels
  case KNL_NAV_S_U: _knl_mm = Kernel3d(KNL_NAV_S_U, _knl.coefs()); break;
  case KNL_NAV_D_U: _knl_mm = Kernel3d(KNL_NAV_S_U, _knl.coefs()); break;
	 //others
  case KNL_SQRTLAP: _knl_mm = Kernel3d(KNL_SQRTLAP, _knl.coefs()); break;
  case KNL_EXP    : _knl_mm = Kernel3d(KNL_EXP    , _knl.coefs()); break;
  default: iA(0);
  }
  _mul_mm = 1; //for the time being
  
  _matmgnt  = MatMgnt3d::getmmptr(_knl_mm, _np);
  //3. self setup
  iC( srcData() );
  iC( trgData() );
  //-----------------------------------------------------
  return (0);
}
// ---------------------------------------------------------------------- 
int FMM3d::srcData()
{
  //1. create vecs
  int srcNodeCnt = _let->srcNodeCnt();
  int srcExaCnt = _let->srcExaCnt();
  _srcExaPos.resize(dim(), srcExaCnt);
  _srcExaNor.resize(dim(), srcExaCnt);
  _srcExaDen.resize(srcExaCnt * srcDOF());
  _srcUpwEquDen.resize(srcNodeCnt * datSze(UE));
  _srcUpwChkVal.resize(srcNodeCnt * datSze(UC));
  
  //2. gather the Position using the Pos scatter
  vector<int> ordVec;  iC( _let->upwOrderCollect(ordVec) );
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if(_let->tag(gNodeIdx) & LET_SRCNODE) { //contributor
		if(_let->terminal(gNodeIdx)==true) { //terminal cbtr
		  DblNumMat srcExaPos(this->srcExaPos(gNodeIdx));
		  DblNumMat srcExaNor(this->srcExaNor(gNodeIdx));
		  vector<int>& curVecIdxs = _let->node(gNodeIdx).srcOwnVecIdxs();
		  for(int k=0; k < curVecIdxs.size(); k++) {
			 int poff = curVecIdxs[k];
			 for(int d=0; d < dim(); d++) {
				srcExaPos(d,k) = (*_srcPos)(d,poff);//Pos
				srcExaNor(d,k) = (*_srcNor)(d,poff);//Nor
			 }
		  }
		}
	 }
  }
  return (0);
}
// ---------------------------------------------------------------------- 
int FMM3d::trgData()
{
  //1. create vecs
  int trgNodeCnt = _let->trgNodeCnt();
  int trgExaCnt = _let->trgExaCnt();
  _trgExaPos.resize(dim(), trgExaCnt);
  _trgExaVal.resize(trgExaCnt * trgDOF());
  _trgDwnEquDen.resize(trgNodeCnt * datSze(DE));
  _trgDwnChkVal.resize(trgNodeCnt * datSze(DC));
  
  //2. gather data from _trgPos
  vector<int> ordVec; iC( _let->upwOrderCollect(ordVec) );
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if(_let->tag(gNodeIdx) & LET_TRGNODE) {
		if(_let->terminal(gNodeIdx)==true) {
		  DblNumMat trgExaPos(this->trgExaPos(gNodeIdx));
		  vector<int>& curVecIdxs = _let->node(gNodeIdx).trgOwnVecIdxs();
		  for(int k=0; k<curVecIdxs.size(); k++) {
			 int poff = curVecIdxs[k];
			 for(int d=0; d<dim(); d++)
				trgExaPos(d,k) = (*_trgPos)(d, poff);
		  }
		}
	 }
  }
  
  //3. allocate ENExt
  _nodeVec.resize( _let->nodeVec().size() );
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if(_let->tag(gNodeIdx) & LET_TRGNODE) {
		//V
		Let3d::Node& gg = _let->node(gNodeIdx);
		_nodeVec[gNodeIdx].VinNum() = gg.Vnodes().size();
		for(vector<int>::iterator vi=gg.Vnodes().begin(); vi!=gg.Vnodes().end(); vi++) {
		  _nodeVec[*vi].VotNum() ++;
		}
	 }
  }
  
  return (0);
}
