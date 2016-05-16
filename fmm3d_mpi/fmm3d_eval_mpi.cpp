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
#include "fmm3d_mpi.hpp"
#include "common/vecmatop.hpp"

using std::cerr;
using std::endl;

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "FMM3d_MPI::evaluate"
int FMM3d_MPI::evaluate(Vec srcDen, Vec trgVal)
{
  //begin  //ebiLogInfo( "multiply.............");
  //-----------------------------------
  //cerr<<"fmm src and trg numbers "<<pglbnum(_srcPos)<<" "<<pglbnum(_trgPos)<<endl;
  int tmp;
  pC( VecGetSize(srcDen,&tmp) );  pA(tmp==srcDOF()*procGlbNum(_srcPos));
  pC( VecGetSize(trgVal,&tmp) );  pA(tmp==trgDOF()*procGlbNum(_trgPos));
  
  int srcDOF = this->srcDOF();
  int trgDOF = this->trgDOF();
  
  //1. zero out vecs.  This includes all global, contributor, user, evaluator vectors.
  PetscScalar zero=0.0;
  pC( VecSet(trgVal, zero) );
  pC( VecSet(_glbSrcExaDen, zero) );
  pC( VecSet(_glbSrcUpwEquDen, zero) );
  pC( VecSet(_ctbSrcExaDen, zero) );
  pC( VecSet(_ctbSrcUpwEquDen, zero) );
  pC( VecSet(_ctbSrcUpwChkVal, zero) );
  pC( VecSet(_usrSrcExaDen, zero) );
  pC( VecSet(_usrSrcUpwEquDen, zero) );
  pC( VecSet(_evaTrgExaVal, zero) );  
  pC( VecSet(_evaTrgDwnEquDen, zero) );
  pC( VecSet(_evaTrgDwnChkVal, zero) );
  
  clock_t ck0, ck1;
  vector<int> ordVec;
  pC( _let->upwOrderCollect(ordVec) ); //BOTTOM UP collection of nodes

  //2. for contributors, load exact densities
  int procLclStart, procLclEnd; _let->procLclRan(_srcPos, procLclStart, procLclEnd);
  double* darr; pC( VecGetArray(srcDen, &darr) );
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if(_let->node(gNodeIdx).tag() & LET_CBTRNODE) {
		if(_let->terminal(gNodeIdx)==true) {
		  DblNumVec ctbSrcExaDen(this->ctbSrcExaDen(gNodeIdx));
		  vector<int>& curVecIdxs = _let->node(gNodeIdx).ctbSrcOwnVecIdxs();
		  for(int k=0; k<curVecIdxs.size(); k++) {
			 int poff = curVecIdxs[k] - procLclStart;
			 for(int d=0; d<srcDOF; d++) {
				ctbSrcExaDen(k*srcDOF+d) = darr[poff*srcDOF+d];
			 }
		  }
		}
	 }
  }
  pC( VecRestoreArray(srcDen, &darr) );
    
  //SCATTER
  pC( VecScatterBegin( _ctbSrcExaDen, _glbSrcExaDen,    ADD_VALUES, SCATTER_FORWARD, _ctb2GlbSrcExaDen) );
  
  //3. up computation
    for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if( _let->node(gNodeIdx).tag() & LET_CBTRNODE) {
		if(_let->depth(gNodeIdx)>=0) {
		  DblNumVec ctbSrcUpwChkValgNodeIdx(ctbSrcUpwChkVal(gNodeIdx));
		  DblNumVec ctbSrcUpwEquDengNodeIdx(ctbSrcUpwEquDen(gNodeIdx));
		  if(_let->terminal(gNodeIdx)==true) {
			 //S2M
			 pC( SrcEqu2UpwChk_dgemv(ctbSrcExaPos(gNodeIdx), ctbSrcExaNor(gNodeIdx), _let->center(gNodeIdx), _let->radius(gNodeIdx), ctbSrcExaDen(gNodeIdx), ctbSrcUpwChkValgNodeIdx) );
		  } else {
			 //M2M
			 for(int a=0; a<2; a++) for(int b=0; b<2; b++) for(int c=0; c<2; c++) {
				Index3 idx(a,b,c);
				int chi = _let->child(gNodeIdx, idx);
				if(_let->node(chi).tag() & LET_CBTRNODE) {
				  pC( _matmgnt->UpwEqu2UpwChk_dgemv(_let->depth(chi)+_rootLevel, idx, ctbSrcUpwEquDen(chi), ctbSrcUpwChkValgNodeIdx) );
				}
			 }
		  }
		  //M2M
		  pC( _matmgnt->UpwChk2UpwEqu_dgemv(_let->depth(gNodeIdx)+_rootLevel, ctbSrcUpwChkValgNodeIdx, ctbSrcUpwEquDengNodeIdx) );
		}
	 }
  }
  
  //4. vectbscatters
  //SCATTER
  pC( VecScatterBegin( _ctbSrcUpwEquDen, _glbSrcUpwEquDen,    ADD_VALUES, SCATTER_FORWARD, _ctb2GlbSrcUpwEquDen) );
  pC( VecScatterEnd(   _ctbSrcExaDen, _glbSrcExaDen,    ADD_VALUES, SCATTER_FORWARD, _ctb2GlbSrcExaDen) );
  pC( VecScatterEnd(   _ctbSrcUpwEquDen, _glbSrcUpwEquDen,    ADD_VALUES, SCATTER_FORWARD, _ctb2GlbSrcUpwEquDen) );
  
  pC( VecScatterBegin( _glbSrcExaDen, _usrSrcExaDen, INSERT_VALUES, SCATTER_REVERSE, _usr2GlbSrcExaDen) );
  pC( VecScatterBegin( _glbSrcUpwEquDen, _usrSrcUpwEquDen, INSERT_VALUES, SCATTER_REVERSE, _usr2GlbSrcUpwEquDen) );
  pC( VecScatterEnd(   _glbSrcExaDen, _usrSrcExaDen, INSERT_VALUES, SCATTER_REVERSE, _usr2GlbSrcExaDen) );
  
  ordVec.clear();  pC( _let->dwnOrderCollect(ordVec) );
  //U
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if( _let->node(gNodeIdx).tag() & LET_EVTRNODE) {
		if( _let->terminal(gNodeIdx)==true ) { //terminal
		  DblNumVec evaTrgExaValgNodeIdx(evaTrgExaVal(gNodeIdx));
		  DblNumMat evaTrgExaPosgNodeIdx(evaTrgExaPos(gNodeIdx));
		  for(vector<int>::iterator vi=_let->node(gNodeIdx).Unodes().begin(); vi!=_let->node(gNodeIdx).Unodes().end(); vi++) {
			 //S2T
			 pC( SrcEqu2TrgChk_dgemv(usrSrcExaPos(*vi), usrSrcExaNor(*vi), evaTrgExaPosgNodeIdx, usrSrcExaDen(*vi), evaTrgExaValgNodeIdx) );
		  }
		}
	 }
  }
  pC( VecScatterEnd(   _glbSrcUpwEquDen, _usrSrcUpwEquDen, INSERT_VALUES, SCATTER_REVERSE, _usr2GlbSrcUpwEquDen) );
  
  //V
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if( _let->node(gNodeIdx).tag() & LET_EVTRNODE ) { //evaluator		
		Point3 gNodeIdxctr(_let->center(gNodeIdx));
		double D = 2.0 * _let->radius(gNodeIdx);
		
		DblNumVec evaTrgDwnChkVal(this->evaTrgDwnChkVal(gNodeIdx));
		for(vector<int>::iterator vi=_let->node(gNodeIdx).Vnodes().begin(); vi!=_let->node(gNodeIdx).Vnodes().end(); vi++) {
		  Point3 victr(_let->center(*vi));
		  Index3 idx;		  for(int d=0; d<dim(); d++)			 idx(d) = int(round( (victr[d]-gNodeIdxctr[d])/D ));
		  
		  Node& srcnode = node(*vi);
		  Node& trgnode = node(gNodeIdx);
		  if(srcnode.vLstOthCnt()==0) {
			 srcnode.effDen().resize( _matmgnt->effDatSze(UE) );			 setvalue(srcnode.effDen(), 0.0);//1. resize effDen
			 pC( _matmgnt->plnDen2EffDen(_let->depth(gNodeIdx)+_rootLevel, usrSrcUpwEquDen(*vi),  srcnode.effDen()) );			 //2. transform from UpwEquDen to effDen
		  }
		  if(trgnode.vLstInCnt()==0) {
			 trgnode.effVal().resize( _matmgnt->effDatSze(DC) );			 setvalue(trgnode.effVal(), 0.0);			 //1. resize effVal
		  }
		  //M2L		  
		  pC( _matmgnt->UpwEqu2DwnChk_dgemv(_let->depth(gNodeIdx)+_rootLevel, idx, srcnode.effDen(), trgnode.effVal()) );
		  
		  srcnode.vLstOthCnt()++;
		  trgnode.vLstInCnt()++;
		  if(srcnode.vLstOthCnt()==srcnode.vLstOthNum()) {
			 srcnode.effDen().resize(0);			 //1. resize effDen to 0
			 srcnode.vLstOthCnt()=0;
		  }
		  if(trgnode.vLstInCnt()==trgnode.vLstInNum()) {
			 pC( _matmgnt->effVal2PlnVal(_let->depth(gNodeIdx)+_rootLevel, trgnode.effVal(), evaTrgDwnChkVal) );			 //1. transform from effval to DwnChkVal
			 trgnode.effVal().resize(0); //2. resize effVal to 0
			 trgnode.vLstInCnt()=0;
		  }
		}
	 }
  }

  //W
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if( _let->node(gNodeIdx).tag() & LET_EVTRNODE) {
		if( _let->terminal(gNodeIdx)==true ) {
		  DblNumVec evaTrgExaVal_gNodeIdx(this->evaTrgExaVal(gNodeIdx));
		  for(vector<int>::iterator vi=_let->node(gNodeIdx).Wnodes().begin(); vi!=_let->node(gNodeIdx).Wnodes().end(); vi++) {
			 if(_let->terminal(*vi) && _let->node(*vi).usrSrcExaNum()*srcDOF<_matmgnt->plnDatSze(UE)) { //use Exa instead
				//S2T
				pC( SrcEqu2TrgChk_dgemv(usrSrcExaPos(*vi), usrSrcExaNor(*vi), evaTrgExaPos(gNodeIdx), usrSrcExaDen(*vi), evaTrgExaVal_gNodeIdx) );
			 } else {
				//M2T
				int vni = *vi;		
				pC( UpwEqu2TrgChk_dgemv(_let->center(vni), _let->radius(vni), evaTrgExaPos(gNodeIdx), usrSrcUpwEquDen(*vi), evaTrgExaVal_gNodeIdx) );
			 }
		  }
		}
	 }
  }

  //X
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if( _let->node(gNodeIdx).tag() & LET_EVTRNODE) {
		DblNumVec evaTrgExaVal_gNodeIdx(evaTrgExaVal(gNodeIdx));
		DblNumVec evaTrgDwnChkVal_gNodeIdx(evaTrgDwnChkVal(gNodeIdx));
		for(vector<int>::iterator vi=_let->node(gNodeIdx).Xnodes().begin(); vi!=_let->node(gNodeIdx).Xnodes().end(); vi++) {
		  if(_let->terminal(gNodeIdx) && _let->node(gNodeIdx).evaTrgExaNum()*trgDOF<_matmgnt->plnDatSze(DC)) { //use Exa instead
			 pC( SrcEqu2TrgChk_dgemv(usrSrcExaPos(*vi), usrSrcExaNor(*vi), evaTrgExaPos(gNodeIdx), usrSrcExaDen(*vi), evaTrgExaVal_gNodeIdx) );
		  } else {
			 //S2L
			 pC( SrcEqu2DwnChk_dgemv(usrSrcExaPos(*vi), usrSrcExaNor(*vi), _let->center(gNodeIdx), _let->radius(gNodeIdx), usrSrcExaDen(*vi), evaTrgDwnChkVal_gNodeIdx) );
		  }
		}
	 }
  }

  //7. combine
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if( _let->node(gNodeIdx).tag() & LET_EVTRNODE ) { //evaluator	
		if(_let->depth(gNodeIdx)>=3) {
		  int pargNodeIdx = _let->parent(gNodeIdx);	
		  Index3 chdidx( _let->path2Node(gNodeIdx)-2 * _let->path2Node(pargNodeIdx) );
		  //L2L
		  DblNumVec evaTrgDwnChkVal_gNodeIdx(evaTrgDwnChkVal(gNodeIdx));
		  pC( _matmgnt->DwnEqu2DwnChk_dgemv(_let->depth(pargNodeIdx)+_rootLevel, chdidx, evaTrgDwnEquDen(pargNodeIdx), evaTrgDwnChkVal_gNodeIdx) );
		}
		if(_let->depth(gNodeIdx)>=2) {
		  //L2L
		  DblNumVec evaTrgDwnEquDen_gNodeIdx(evaTrgDwnEquDen(gNodeIdx));
		  pC( _matmgnt->DwnChk2DwnEqu_dgemv(_let->depth(gNodeIdx)+_rootLevel, evaTrgDwnChkVal(gNodeIdx), evaTrgDwnEquDen_gNodeIdx) );
		}
		if(_let->terminal(gNodeIdx)) {
		  //L2T
		  DblNumVec evaTrgExaVal_gNodeIdx(evaTrgExaVal(gNodeIdx));
		  pC( DwnEqu2TrgChk_dgemv(_let->center(gNodeIdx), _let->radius(gNodeIdx), evaTrgExaPos(gNodeIdx), evaTrgDwnEquDen(gNodeIdx), evaTrgExaVal_gNodeIdx) );
		}
	 }
  }
  
  //8. save tdtExaVal
  _let->procLclRan(_trgPos, procLclStart, procLclEnd);
  double* varr; pC( VecGetArray(trgVal, &varr) );
  for(int i=0; i<ordVec.size(); i++) {
	 int gNodeIdx = ordVec[i];
	 if( _let->node(gNodeIdx).tag() & LET_EVTRNODE ) {
		if( _let->terminal(gNodeIdx)==true ) {
		  DblNumVec evaTrgExaVal(this->evaTrgExaVal(gNodeIdx));
		  vector<int>& curVecIdxs = _let->node(gNodeIdx).evaTrgOwnVecIdxs();
		  for(int k=0; k<curVecIdxs.size(); k++) {
			 int poff = curVecIdxs[k] - procLclStart;
			 for(int d=0; d<trgDOF; d++) {
				varr[poff*trgDOF+d] = evaTrgExaVal(k*trgDOF+d);
			 }
		  }
		}
	 }
  }
  pC( VecRestoreArray(trgVal, &varr) );
  
  //----------------
  pC( MPI_Barrier(mpiComm()) );  //check vLstInCnt, vLstOthCnt
  return(0);
}



