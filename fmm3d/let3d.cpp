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
#include "let3d.hpp"

using std::min;
using std::max;
using std::set;
using std::queue;
using std::ofstream;
using std::cerr;
using std::endl;
using std::istringstream;

//-----------------------------------------------------------------------
Let3d::Let3d(const string& p):
  ComObject(p), _srcPos(NULL), _trgPos(NULL), _center(0.0), _rootLevel(0),
  _ptsMax(150), _maxLevel(10)
{
}

Let3d::~Let3d()
{
}

// ---------------------------------------------------------------------- 
int Let3d::setup(map<string,string>& optionsMap)
{
  //------------
  map<string,string>::iterator mapindex;
  mapindex = optionsMap.find("-" + prefix() + "ptsMax"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>_ptsMax; }
  mapindex = optionsMap.find("-" + prefix() + "maxLevel"); iA(mapindex!=optionsMap.end());
  { istringstream ss((*mapindex).second); ss>>_maxLevel; }
  //------------
  iC( srcData() );
  iC( trgData() );
  return 0;
}

// ---------------------------------------------------------------------- 
int Let3d::srcData()
{
  //-----------------------------------------
  ///gdata
  DblNumMat& pos = *(_srcPos);  iA( pos.m()==dim() );
  
  vector<Node>& nodeVec = this->nodeVec(); nodeVec.clear();
  vector< vector<int> > vecIdxs;  
  //local src number, the number of src point in each box
  vector<int> lclSrcNumVec; //glb 
  
  nodeVec.push_back( Node(-1,-1, Index3(0,0,0), 0) );
  vecIdxs.push_back( vector<int>() );
  vector<int>& curVecIdxs = vecIdxs[0];
  Point3 bbmin(center()-Point3(radius()));
  Point3 bbmax(center()+Point3(radius()));
  for(int k=0; k<pos.n(); k++) {
	 Point3 tmp(pos.clmdata(k));	 
	 iA(tmp>=bbmin && tmp<=bbmax);//LEXING: IMPORANT
	 curVecIdxs.push_back(k);
  }
  lclSrcNumVec.push_back( curVecIdxs.size() );
  
  int level = 0;
  int arrBeg = 0;
  int arrEnd = 1;
  int arrCnt = 0;
  while(arrBeg < arrEnd) {
	 //1.
	 arrCnt = arrEnd;
	 for(int k=arrBeg; k < arrEnd; k++) {
		//---
		if( lclSrcNumVec[k]>ptsMax() && level<maxLevel()-1 ) {
		  nodeVec[k].child() = arrCnt;
		  arrCnt = arrCnt + pow2(dim());
		  //children's ess		  
		  for(int a=0; a<2; a++) {
			 for(int b=0; b<2; b++) {
				for(int c=0; c<2; c++) {
				  nodeVec.push_back( Node(k,-1, 2*nodeVec[k].path2Node()+Index3(a,b,c), nodeVec[k].depth()+1) ); //par, chd
				  vecIdxs.push_back( vector<int>() );
				  lclSrcNumVec.push_back( 0 );
				}
			 }
		  }
		  //children's vector of indices
		  Point3 centerCurNode( center(k) ); //get center of current node
		  for(vector<int>::iterator vecIdxsIt=vecIdxs[k].begin(); vecIdxsIt!=vecIdxs[k].end(); vecIdxsIt++) {
			 Point3 tmp(pos.clmdata(*vecIdxsIt));
			 Index3 idx;
			 for(int j=0; j<dim(); j++){
				idx(j) = (tmp(j) >= centerCurNode(j));
			 }
			 int chdGNodeIdx = child(k, idx);
			 vecIdxs[chdGNodeIdx].push_back(*vecIdxsIt);
		  }
		  vecIdxs[k].clear(); //VERY IMPORTANT
		  //children's lsm		  
		  for(int a=0; a<2; a++) {
			 for(int b=0; b<2; b++) {
				for(int c=0; c<2; c++) {
				  int chdGNodeIdx = child( k, Index3(a,b,c) );
				  lclSrcNumVec[chdGNodeIdx] = vecIdxs[chdGNodeIdx].size();
				}
			 }
		  }
		}
	 }
	 level++;
	 arrBeg = arrEnd;
	 arrEnd = arrCnt;
  }
  _level = level; //SET LEVEL

  //ordering of the boxes, in top-down or bottom-up fashion
  vector<int> orderBoxesVec	;  iC( dwnOrderCollect(orderBoxesVec	) );
  //set other parts of essvec
  int cnt = 0;
  int sum = 0;
  for(int i=0; i < orderBoxesVec	.size(); i++) {
	 int gNodeIdx = orderBoxesVec	[i];
	 if(lclSrcNumVec[gNodeIdx]>0) {
		nodeVec[gNodeIdx].tag() = nodeVec[gNodeIdx].tag() | LET_SRCNODE;
		nodeVec[gNodeIdx].srcNodeIdx() = cnt;
		cnt++;
		if(nodeVec[gNodeIdx].child()==-1) {
		  nodeVec[gNodeIdx].srcExaBeg() = sum;
		  nodeVec[gNodeIdx].srcExaNum() = lclSrcNumVec[gNodeIdx];
		  sum += lclSrcNumVec[gNodeIdx];
		  nodeVec[gNodeIdx].srcOwnVecIdxs() = vecIdxs[gNodeIdx];
		}
	 }
  }
  _srcNodeCnt = cnt;  _srcExaCnt = sum; //SET S cnts
  return 0;
}

// ---------------------------------------------------------------------- 
int Let3d::trgData()
{
  //-----------------------------------------
  //edata
  DblNumMat& pos = *(_trgPos);  iA( pos.m()==dim() );
  
  vector<Node>& nodeVec = this->nodeVec();
  vector< vector<int> > vecIdxs; vecIdxs.resize(nodeVec.size() );
  vector<int> lclSrcNumVec;           lclSrcNumVec.resize(nodeVec.size(), 0);
  
  vector<int>& curVecIdxs = vecIdxs[0];
  Point3 bbmin(center()-Point3(radius()));
  Point3 bbmax(center()+Point3(radius()));  //cerr<<" bbx "<<bbmin<<" "<<bbmax<<endl;
  for(int k=0; k < pos.n(); k++) {
	 Point3 tmp(pos.clmdata(k));
	 iA(tmp>=bbmin && tmp<=bbmax);	 //LEXING: IMPORTANT
	 curVecIdxs.push_back(k);
  }
  lclSrcNumVec[0] = curVecIdxs.size();
  
  vector<int> orderBoxesVec	;  iC( dwnOrderCollect(orderBoxesVec	) );
  for(int i=0; i < orderBoxesVec	.size(); i++) {
	 int gNodeIdx = orderBoxesVec[i];
	 Node& curNode = nodeVec[gNodeIdx];
	 vector<int>& curVecIdxs = vecIdxs[gNodeIdx];	 
	 if(curNode.child()!=-1) { //not terminal
		//children's vecIdxs
		Point3 curCenter( center(gNodeIdx) );
		for(vector<int>::iterator curVecIdxsIt=curVecIdxs.begin();curVecIdxsIt !=curVecIdxs.end(); curVecIdxsIt++) {
		  Point3 tmp(pos.clmdata(*curVecIdxsIt));
		  Index3 idx;
		  for(int j=0; j<dim(); j++)
			 idx(j) = (tmp(j)>=curCenter(j));
		  int chdGNodeIdx = child(gNodeIdx, idx);
		  vector<int>& chdVecIdxs = vecIdxs[chdGNodeIdx];
		  chdVecIdxs.push_back(*curVecIdxsIt);
		}
		curVecIdxs.clear(); //VERY IMPORTANT
		//children's lsm		//		for(int ord=0; ord<8; ord++) {
		for(int a=0; a<2; a++) {
		  for(int b=0; b<2; b++) {
			 for(int c=0; c<2; c++) {
				int chdGNodeIdx = child(gNodeIdx, Index3(a,b,c));
				lclSrcNumVec[chdGNodeIdx] = vecIdxs[chdGNodeIdx].size();
			 }
		  }
		}
	 }
  }
  //set EVTR
  int cnt = 0;
  int sum = 0;
  for(int i=0; i<orderBoxesVec	.size(); i++) {
	 int gNodeIdx = orderBoxesVec	[i];
	 if(lclSrcNumVec[gNodeIdx]>0) { //evtr node
		nodeVec[gNodeIdx].tag() = nodeVec[gNodeIdx].tag() | LET_TRGNODE;
		nodeVec[gNodeIdx].trgNodeIdx() = cnt;
		cnt ++;
		if(nodeVec[gNodeIdx].child()==-1) { //terminal
		  nodeVec[gNodeIdx].trgExaBeg() = sum;
		  nodeVec[gNodeIdx].trgExaNum() = lclSrcNumVec[gNodeIdx];
		  sum += lclSrcNumVec[gNodeIdx];
		  nodeVec[gNodeIdx].trgOwnVecIdxs() = vecIdxs[gNodeIdx];
		}
	 }
  }
  _trgNodeCnt = cnt;  _trgExaCnt = sum;
  
  //set USER
  for(int i=0; i < orderBoxesVec	.size(); i++) {
	 int gNodeIdx = orderBoxesVec	[i];
	 if(nodeVec[gNodeIdx].tag() & LET_TRGNODE) { //a evtr		
		iC( calgnext(gNodeIdx) );
	 }
  }
  //cerr<<usndecnt<<" "<<usextcnt<<endl;
  return 0;
}

// ---------------------------------------------------------------------- 
int Let3d::calgnext(int gNodeIdx)
{
  vector<Node>&  nodeVec = this->nodeVec();
  
  set<int> Uset, Vset, Wset, Xset;
  int curGNodeIdx = gNodeIdx;
  if(root(curGNodeIdx)==false) {
	 int parGNodeIdx = parent(curGNodeIdx);
	 
	 Index3 minIdx(0);
	 Index3 maxIdx(pow2(depth(curGNodeIdx)));

	 for(int i=-2; i<4; i++) {
		for(int j=-2; j<4; j++)	{
		  for(int k=-2; k<4; k++) {
			 Index3 tryPath( 2*path2Node(parGNodeIdx) + Index3(i,j,k) );
			 if(tryPath >= minIdx && tryPath <  maxIdx && tryPath != path2Node(curGNodeIdx)) {	
				int resGNodeIdx = findgnt(depth(curGNodeIdx), tryPath);
				bool adj = adjacent(resGNodeIdx, curGNodeIdx);
				if( depth(resGNodeIdx) < depth(curGNodeIdx) ) {
				  if(adj){
					 if(terminal(curGNodeIdx)){
						Uset.insert(resGNodeIdx);
					 }
					 else { ; }
				  }
				  else {
					 Xset.insert(resGNodeIdx);
				  }
				}
				if( depth(resGNodeIdx)==depth(curGNodeIdx) ) {
				  if(!adj) {
					 Index3 bb(path2Node(resGNodeIdx)-path2Node(curGNodeIdx));
					 assert( bb.linfty()<=3 );
					 Vset.insert(resGNodeIdx);
				  }
				  else {
					 if(terminal(curGNodeIdx)) {
						queue<int> rest;
						rest.push(resGNodeIdx);
						while(rest.empty()==false) {
						  int fntGNodeIdx = rest.front(); rest.pop();					 //int fntgNodeIdx = fntgnt.gNodeIdx();
						  if(adjacent(fntGNodeIdx, curGNodeIdx)==false) {
							 Wset.insert( fntGNodeIdx );
						  }
						  else {
							 if(terminal(fntGNodeIdx)) {
								Uset.insert(fntGNodeIdx);
							 }
							 else { 
								for(int a=0; a<2; a++) {
								  for(int b=0; b<2; b++) {
									 for(int c=0; c<2; c++) {
										rest.push( child(fntGNodeIdx, Index3(a,b,c)) );
									 }
								  }
								}
							 }
						  }
						}
					 }
				  }
				}
			 }
		  }
		}
	 }
  }
  if(terminal(curGNodeIdx))
	 Uset.insert(curGNodeIdx);
  
  for(set<int>::iterator si=Uset.begin(); si!=Uset.end(); si++)
	 if(nodeVec[*si].tag() & LET_SRCNODE)		nodeVec[gNodeIdx].Unodes().push_back(*si);
  for(set<int>::iterator si=Vset.begin(); si!=Vset.end(); si++)
	 if(nodeVec[*si].tag() & LET_SRCNODE)		nodeVec[gNodeIdx].Vnodes().push_back(*si);
  for(set<int>::iterator si=Wset.begin(); si!=Wset.end(); si++)
	 if(nodeVec[*si].tag() & LET_SRCNODE)		nodeVec[gNodeIdx].Wnodes().push_back(*si);
  for(set<int>::iterator si=Xset.begin(); si!=Xset.end(); si++)
	 if(nodeVec[*si].tag() & LET_SRCNODE)		nodeVec[gNodeIdx].Xnodes().push_back(*si);
  
  return (0);
}
// ---------------------------------------------------------------------- 
int Let3d::dwnOrderCollect(vector<int>& orderBoxesVec	)
{
  orderBoxesVec.clear();
  for(int i=0; i < nodeVec().size(); i++)
	 orderBoxesVec.push_back(i);
  iA(orderBoxesVec.size()==nodeVec().size());
  return 0;
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "Let3d::upwardOrderCollect"
int Let3d::upwOrderCollect(vector<int>& orderBoxesVec	)
{
  orderBoxesVec.clear();
  for(int i=nodeVec().size()-1; i>=0; i--)
	 orderBoxesVec.push_back(i);
  iA(orderBoxesVec.size()==nodeVec().size());
  return 0;
}
// ---------------------------------------------------------------------- 
int Let3d::child(int gNodeIdx, const Index3& idx)
{
  assert(idx>=Index3(0) && idx<Index3(2));
  return node(gNodeIdx).child() + (idx(0)*4+idx(1)*2+idx(2));
}
Point3 Let3d::center(int gNodeIdx) //center of a node
{
  Point3 ll( center() - Point3(radius()) );
  int tmp = pow2(depth(gNodeIdx));
  Index3 pathLcl(path2Node(gNodeIdx));
  Point3 res;
  for(int d=0; d<dim(); d++) {
	 res(d) = ll(d) + (2*radius()) * (pathLcl(d)+0.5) / double(tmp);
  }
  return res;
}
double Let3d::radius(int gNodeIdx) //radius of a node
{
  return radius()/double(pow2(depth(gNodeIdx)));
}
// ---------------------------------------------------------------------- 
int Let3d::findgnt(int wntdepth, const Index3& wntpath)
{
  int cur = 0;  //cerr<<"GOOD "<<path(cur)<<"     ";
  Index3 leftpath(wntpath);
  while(depth(cur)<wntdepth && terminal(cur)==false) {
	 int dif = wntdepth-depth(cur);
	 int tmp = pow2(dif-1);
	 Index3 choice( leftpath/tmp );
	 leftpath -= choice*tmp;
	 cur = child(cur, choice);	 //cur = child(cur, IdxIter(0,2,true,choice) );	 //cerr<<path(cur)<<"["<<choice<<" "<<tmp<<"]"<<"     ";
  }  //cerr<<endl;
  return cur;
}
// ---------------------------------------------------------------------- 
bool Let3d::adjacent(int me, int yo)
{
  int md = max(depth(me),depth(yo));
  Index3 one(1);
  Index3 mecenter(  (2*path2Node(me)+one) * pow2(md - depth(me))  );
  Index3 yocenter(  (2*path2Node(yo)+one) * pow2(md - depth(yo))  );
  int meradius = pow2(md - depth(me));
  int yoradius = pow2(md - depth(yo));
  Index3 dif( abs(mecenter-yocenter) );
  int radius  = meradius + yoradius;
  return
	 ( dif <= Index3(radius) ) && //not too far
	 ( dif.linfty() == radius ); //at least one edge touch
}

// ---------------------------------------------------------------------- 
int Let3d::print()
{
  cerr<<_nodeVec.size()<<" "<<_level<<endl;
  cerr<<_srcNodeCnt <<" "<< _srcExaCnt << endl;
  cerr<<_trgNodeCnt <<" "<< _trgExaCnt << endl;

  for(int i=0; i<_nodeVec.size(); i++) {
	 cerr<<node(i).parent()<<" "<<node(i).child()<<endl;
  }
  
  for(int i=0; i<_nodeVec.size(); i++) {
	 //cerr<<node(i).Unodes().size()<<" "<<node(i).Vnodes().size()<<" "<<node(i).Wnodes().size()<<" "<<node(i).Xnodes().size()<<" "<<endl;
	 Node& n = node(i);
	 cerr<<n.srcNodeIdx()<<" "<<n.srcExaBeg()<<" "<<n.srcExaNum()<<" "
		  <<n.trgNodeIdx()<<" "<<n.trgExaBeg()<<" "<<n.trgExaNum()<<endl;
  }
  return 0;
}
