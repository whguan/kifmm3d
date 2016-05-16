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

// ---------------------------------------------------------------------- 
int FMM3d::evaluate(const DblNumVec& srcDen, DblNumVec& trgVal)
{
  _matmgnt->report();
  //-----------------------------------
  iA(srcDen.m()==srcDOF()*(*_srcPos).n());  iA(trgVal.m()==trgDOF()*(*_trgPos).n());
  
  //cerr<<"fmm src and trg numbers "<<pglbnum(_srcPos)<<" "<<pglbnum(_trgPos)<<endl;
  int srcDOF = this->srcDOF();
  int trgDOF = this->trgDOF();
  
  //1. zero out Vecs
  setvalue(trgVal, 0.0);
  
  setvalue(_srcExaDen, 0.0);
  setvalue(_srcUpwEquDen, 0.0);
  setvalue(_srcUpwChkVal, 0.0);
  
  setvalue(_trgExaVal, 0.0);
  setvalue(_trgDwnEquDen, 0.0);
  setvalue(_trgDwnChkVal, 0.0);
  clock_t ck0, ck1;
  //CLOCKING;
  ck0 = clock();
  vector<int> ordVec; iC( _let->upwOrderCollect(ordVec) ); //BOTTOM UP

  //2. for cbtr, load ExaDen
  ck0 = clock();
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if(_let->tag(gNodeIdx) & LET_SRCNODE) {
		if(_let->terminal(gNodeIdx)==true) {
		  DblNumVec srcExaDen(this->srcExaDen(gNodeIdx));
		  vector<int>& curVecIdxs = _let->node(gNodeIdx).srcOwnVecIdxs();
		  for(int k=0; k<curVecIdxs.size(); k++) {
			 int poff = curVecIdxs[k];
			 for(int d=0; d<srcDOF; d++) {
				srcExaDen(k*srcDOF+d) = srcDen(poff*srcDOF+d);
			 }
		  }
		}
	 }
  }
  //ck1 = clock();  cerr<<"load  "<<double(ck1-ck0)/CLOCKS_PER_SEC<<endl;  
  
  //3. up computation
  ck0 = clock();
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if(_let->tag(gNodeIdx) & LET_SRCNODE) {		//GNTra gnt = _let->gNodeIdx2gnt(gNodeIdx);
		//if(_let->depth(gNodeIdx)>=2) {
		if(_let->depth(gNodeIdx)>=0) {
		  DblNumVec srcUpwChkValgNodeIdx(srcUpwChkVal(gNodeIdx));
		  DblNumVec srcUpwEquDengNodeIdx(srcUpwEquDen(gNodeIdx));
		  if(_let->terminal(gNodeIdx)==true) {
			 //S2M - Source -> Multipole Exapnsion
			 iC( SrcEqu2UpwChk_dgemv(srcExaPos(gNodeIdx), srcExaNor(gNodeIdx), _let->center(gNodeIdx), _let->radius(gNodeIdx), srcExaDen(gNodeIdx), srcUpwChkValgNodeIdx) );
		  } else {
			 //M2M - Multipole -> Multipole
			 for(int a=0; a<2; a++) {
				for(int b=0; b<2; b++) {
				  for(int c=0; c<2; c++) {
					 Index3 idx(a,b,c);
					 int chi = _let->child(gNodeIdx, idx);
					 if(_let->tag(chi) & LET_SRCNODE) {
						iC( _matmgnt->UpwEqu2UpwChk_dgemv(_let->depth(chi)+_rootLevel, idx, srcUpwEquDen(chi), srcUpwChkValgNodeIdx) );
					 }
				  }
				}
			 }
		  }
		  //M2M - Multipole -> Multipole
		  iC( _matmgnt->UpwChk2UpwEqu_dgemv(_let->depth(gNodeIdx)+_rootLevel, srcUpwChkValgNodeIdx, srcUpwEquDengNodeIdx) );
		}
	 }
  }
  ordVec.clear();  iC( _let->dwnOrderCollect(ordVec) );
  //U - list contribution calculation
  ck0 = clock();
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if(_let->tag(gNodeIdx) & LET_TRGNODE) { //evaluator
		if( _let->terminal(gNodeIdx)==true ) { //terminal	
		  Let3d::Node& curNode = _let->node(gNodeIdx);
		  DblNumVec trgExaValgNodeIdx(trgExaVal(gNodeIdx));
		  DblNumMat trgExaPosgNodeIdx(trgExaPos(gNodeIdx));
		  for(vector<int>::iterator vi=curNode.Unodes().begin(); vi!=curNode.Unodes().end(); vi++) {
			 //S2T - source -> target
			 iC( SrcEqu2TrgChk_dgemv(srcExaPos(*vi), srcExaNor(*vi), trgExaPosgNodeIdx, srcExaDen(*vi), trgExaValgNodeIdx) );
		  }
		}
	 }
  }
  //V - list contribution calculation
  ck0 = clock();
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if( _let->tag(gNodeIdx) & LET_TRGNODE) { //eValuator		//GNTra gnt = _let->gNodeIdx2gnt(gNodeIdx);
		Point3 gNodeIdxCtr(_let->center(gNodeIdx));
		double D = 2.0 * _let->radius(gNodeIdx);
		DblNumVec trgDwnChkVal(this->trgDwnChkVal(gNodeIdx));
		Let3d::Node& curNode = _let->node(gNodeIdx);
		for(vector<int>::iterator vi=curNode.Vnodes().begin(); vi!=curNode.Vnodes().end(); vi++) {
		  Point3 viCtr(_let->center(*vi));
		  Index3 idx;
		  for(int d=0; d<dim(); d++){
			 idx(d) = int(round( (viCtr[d]-gNodeIdxCtr[d])/D ));
		  }
		  Node& srcPtr = node(*vi);
		  Node& trgPtr = node(gNodeIdx);
		  if(srcPtr.VotCnt()==0) {
			 srcPtr.effDen().resize( _matmgnt->effDatSze(UE) );			 setvalue(srcPtr.effDen(), 0.0);//1. resize effDen
			 iC( _matmgnt->plnDen2EffDen(_let->depth(gNodeIdx)+_rootLevel, srcUpwEquDen(*vi),  srcPtr.effDen()) );			 //2. transform from upeDen to effDen
		  }
		  if(trgPtr.VinCnt()==0) {
			 trgPtr.effVal().resize( _matmgnt->effDatSze(DC) );			 setvalue(trgPtr.effVal(), 0.0);			 //1. resize effVal
		  }
		  //M2L - multipole -> local
		  iC( _matmgnt->UpwEqu2DwnChk_dgemv(_let->depth(gNodeIdx)+_rootLevel, idx, srcPtr.effDen(), trgPtr.effVal()) );
		  
		  srcPtr.VotCnt()++;
		  trgPtr.VinCnt()++;
		  if(srcPtr.VotCnt()==srcPtr.VotNum()) {
			 srcPtr.effDen().resize(0);			 //1. resize effDen to 0
			 srcPtr.VotCnt()=0;
		  }
		  if(trgPtr.VinCnt()==trgPtr.VinNum()) {
			 iC( _matmgnt->effVal2PlnVal(_let->depth(gNodeIdx)+_rootLevel, trgPtr.effVal(), trgDwnChkVal) );			 //1. transform from effVal to dncVal
			 trgPtr.effVal().resize(0); //2. resize effVal to 0
			 trgPtr.VinCnt()=0;
		  }
		}
	 }
  }
  //W - list contrubtion calculation
  ck0 = clock();
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if( _let->tag(gNodeIdx) & LET_TRGNODE ) {
		if( _let->terminal(gNodeIdx)==true ) {
		  DblNumVec trgExaVal_gNodeIdx(this->trgExaVal(gNodeIdx));
		  Let3d::Node& curNode = _let->node(gNodeIdx);
		  for(vector<int>::iterator vi=curNode.Wnodes().begin(); vi!=curNode.Wnodes().end(); vi++) {
			 if(_let->terminal(*vi) && _let->node(*vi).srcExaNum()*srcDOF<_matmgnt->plnDatSze(UE)) { //use Exa instead
				//S2T - source -> target
				iC( SrcEqu2TrgChk_dgemv(srcExaPos(*vi), srcExaNor(*vi), trgExaPos(gNodeIdx), srcExaDen(*vi), trgExaVal_gNodeIdx) );
			 } else {
				//M2T - multipole -> target
				int vni = *vi;				
				iC( UpwEqu2TrgChk_dgemv(_let->center(vni), _let->radius(vni), trgExaPos(gNodeIdx), srcUpwEquDen(*vi), trgExaVal_gNodeIdx) );
			 }
		  }
		}
	 }
  }
  //X - list contrubtion calculation
  ck0 = clock();
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if( _let->tag(gNodeIdx) & LET_TRGNODE) {	
		Let3d::Node& curNode = _let->node(gNodeIdx);
		DblNumVec trgExaVal_gNodeIdx(trgExaVal(gNodeIdx));
		DblNumVec trgDwnChkVal_gNodeIdx(trgDwnChkVal(gNodeIdx));
		for(vector<int>::iterator vi=curNode.Xnodes().begin(); vi!=curNode.Xnodes().end(); vi++) {
		  if(_let->terminal(gNodeIdx) && _let->node(gNodeIdx).trgExaNum()*trgDOF<_matmgnt->plnDatSze(DC)) { //use Exa instead
			 iC( SrcEqu2TrgChk_dgemv(srcExaPos(*vi), srcExaNor(*vi), trgExaPos(gNodeIdx), srcExaDen(*vi), trgExaVal_gNodeIdx) );
		  } else {
			 //S2L - source -> local
			 iC( SrcEqu2DwnChk_dgemv(srcExaPos(*vi), srcExaNor(*vi), _let->center(gNodeIdx), _let->radius(gNodeIdx), srcExaDen(*vi), trgDwnChkVal_gNodeIdx) );
		  }
		}
	 }
  }
  //7. combine
  ck0 = clock();
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if( _let->tag(gNodeIdx) & LET_TRGNODE ) { //eValuator		
		if(_let->depth(gNodeIdx)>=3) {
		  int pargNodeIdx = _let->parent(gNodeIdx);	
		  Index3 chdIdx( _let->path2Node(gNodeIdx)-2 * _let->path2Node(pargNodeIdx) );
		  //L2L - local -> local
		  DblNumVec trgDwnChkVal_gNodeIdx(trgDwnChkVal(gNodeIdx));
		  iC( _matmgnt->DwnEqu2DwnChk_dgemv(_let->depth(pargNodeIdx)+_rootLevel, chdIdx, trgDwnEquDen(pargNodeIdx), trgDwnChkVal_gNodeIdx) );
		}
		if(_let->depth(gNodeIdx)>=2) {
		  //L2L - local -> local
		  DblNumVec trgDwnEquDen_gNodeIdx(trgDwnEquDen(gNodeIdx));
		  iC( _matmgnt->DwnChk2DwnEqu_dgemv(_let->depth(gNodeIdx)+_rootLevel, trgDwnChkVal(gNodeIdx), trgDwnEquDen_gNodeIdx) );
		}
		if(_let->terminal(gNodeIdx)) {
		  //L2T - local -> target
		  DblNumVec trgExaVal_gNodeIdx(trgExaVal(gNodeIdx));
		  iC( DwnEqu2TrgChk_dgemv(_let->center(gNodeIdx), _let->radius(gNodeIdx), trgExaPos(gNodeIdx), trgDwnEquDen(gNodeIdx), trgExaVal_gNodeIdx) );
		}
	 }
  }
  
  //8. save trgExaVal
  ck0 = clock();
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if( _let->tag(gNodeIdx) & LET_TRGNODE ) {
		if( _let->terminal(gNodeIdx)==true ) {
		  DblNumVec trgExaVal(this->trgExaVal(gNodeIdx));
		  vector<int>& curVecIdxs = _let->node(gNodeIdx).trgOwnVecIdxs();
		  for(int k=0; k<curVecIdxs.size(); k++) {
			 int poff = curVecIdxs[k];
			 for(int d=0; d<trgDOF; d++) {
				trgVal(poff*trgDOF+d) = trgExaVal(k*trgDOF+d);
			 }
		  }
		}
	 }
  }
  return (0);
}



