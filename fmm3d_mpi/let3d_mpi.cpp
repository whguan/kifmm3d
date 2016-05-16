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
#include "let3d_mpi.hpp"

using std::min;
using std::max;
using std::set;
using std::queue;
using std::ofstream;
using std::cerr;
using std::endl;

//-----------------------------------------------------------------------
Let3d_MPI::Let3d_MPI(const string& p):
  ComObject_MPI(p), _srcPos(NULL), _trgPos(NULL), _ctr(0.0), _rootLevel(0),
  _ptsMax(150), _maxLevel(10)
{
}

Let3d_MPI::~Let3d_MPI()
{
}

// ----------------------------------------------------------------------
/* Setup builds the basic structures needed for processing source data and building target data */
#undef __FUNCT__
#define __FUNCT__ "Let3d_MPI::setup"
int Let3d_MPI::setup()
{
  //begin
  //--------
  PetscTruth flg = PETSC_FALSE;
  /* Get the maximum number of points per leaf level box and max number of levels in the tree */
  pC( PetscOptionsGetInt( prefix().c_str(), "-ptsMax",   &_ptsMax,   &flg) ); pA(flg==true);
  pC( PetscOptionsGetInt( prefix().c_str(), "-maxLevel", &_maxLevel, &flg) ); pA(flg==true);
  //--------
  pC( srcData() );  
  pC( trgData() );  
  return(0);
}

// ----------------------------------------------------------------------
/* This function computes variables for source data computation
 * as well as constructing and partitioning the local essential tree
 * among processors
 */
#undef __FUNCT__
#define __FUNCT__ "Let3d_MPI::srcData"
int Let3d_MPI::srcData()
{
  /* Set the sample DENSITY to 1/10 th of the original positions */
  double DENSITY = 1.0/10.0;
  
  /* sample data with density - local sample positions vector */
  vector<Point3> lclSamPosVec;
  {
	 Vec pos = _srcPos;
    double* posArr; pC( VecGetArray(pos, &posArr) );
	 /* number of local positions for this processor */
	 int numLclPos = procLclNum(pos);
	 /* number of local sample positions for this processor */
	 int numLclSamPos = int(ceil(procLclNum(pos) * DENSITY));
	 /* sample a fraction of the original sources and
	  * store these points in lclSamPosVec - local sample positions vector */
	 for(int k=0; k<numLclSamPos; k++) {
		int tmpIdx = int(floor(drand48() * numLclPos));
		Point3 tmppt(posArr + tmpIdx*dim());
		lclSamPosVec.push_back( tmppt );
	 }
	 pA(lclSamPosVec.size()==numLclSamPos); /* Verify that the local sample positions vector is the size we want */
	 pC( VecRestoreArray(pos, &posArr) ); /* Restore the full global vector of positions */
  } /* end of this scope - psArr, numLclcPos, numLclSamPos, pos no longer available */

  
  
  /* all procs send sample positions to proc 0 */
  vector<Point3> allSamPosVec; /* will store all sample positions for all processors */
  {
	 int numLcl = lclSamPosVec.size(); /* number = numLclSamPos */
	 int sizLcl = numLcl * dim(); /* local size = number of local sample positions times the dimension (3) */

	 /* vector storing size of local vectors.  i.e., if mpisize()==2
	  * sizLclVec is a vector of length 2, where sizLclVec[i] stores the
	  * number of sample positions at the processor of rank i */
	 vector<int> sizLclVec(mpiSize()); 
	 /* MPI_Gather information:  http://www-unix.mcs.anl.gov/mpi/www/www3/MPI_Gather.html
	 * Here, we're gathering all sizLcl values at the buffer located in memory at &(sizLclVec[0]) onr int at a time */
	 pC( MPI_Gather(&sizLcl, 1, MPI_INT, &(sizLclVec[0]), 1, MPI_INT, 0, mpiComm()) );
	 
	 int sizGlb = 0; /* sizGlb will collect all sizLcl values */
	 for(int i=0; i<sizLclVec.size(); i++){
		sizGlb += sizLclVec[i];
	 }
	 /* number of global sample positions for all processors */
	 int numGlb = sizGlb / dim(); 
	 /* displacement global vector - tells how far into allSamPosVec need to go to get to processor i's samples */
	 vector<int> disGlbVec(mpiSize()); 
	 int cnt = 0;
	 for(int i=0; i<sizLclVec.size(); i++) {
		disGlbVec[i] = cnt;
		cnt += sizLclVec[i];	 }
	 pA(cnt==sizGlb);
	 
	 if(mpiRank()==0)		allSamPosVec.resize(numGlb); /* allSamPosVec stores numGlb points of type Point3 */
	 /* Information on MPI_Gatherv http://www-unix.mcs.anl.gov/mpi/www/www3/MPI_Gatherv.html
	 * Like MPI_Gather - all lclSamPosVec points are sent to the allSamPosVec vector.  sizLclVec[i] and disGlbVec[i] will
	 * indicate how to navigate allSamPosVec for processor i */
	 pC( MPI_Gatherv(&(lclSamPosVec[0]), sizLcl, MPI_DOUBLE, &(allSamPosVec[0]), &(sizLclVec[0]), &(disGlbVec[0]), MPI_DOUBLE, 0, mpiComm()) );
	 
	 lclSamPosVec.clear();
  } /* end of this scope - sizLclcVec and disGlbVec no longer available */
  
  /* processor 0  builds the tree 
	* The tree is build based on sample positions since we will
	* have a max number of points per box.  Sampling positions for
	* LET building makes the build faster */
  vector<Node>& nodeVec = this->nodeVec(); nodeVec.clear();
  if(mpiRank()==0) {
	 /* vector of vector of indices */
	 vector< vector<int> > vecIdxs;
	 /* local source number - number of source points in each box */
	 vector<int> lclSrcNumVec;

	 /* Push the root node.  No parent or child, path2Node = (0,0,0), depth 0 */
	 nodeVec.push_back( Node(-1,-1,Index3(0,0,0), 0) );
	 vecIdxs.push_back( vector<int>() );
	 /* vecIdxs[0] will store all indices put into curVecIdxs */
	 vector<int>& curVecIdxs = vecIdxs[0];
	 Point3 bbmin(ctr()-Point3(radius()));
	 Point3 bbmax(ctr()+Point3(radius()));
	 /* For all sample points, make sure that it is within the root's radius distance from the center of the tree */
	 for(int k=0; k<allSamPosVec.size(); k++) {
		Point3 tmp(allSamPosVec[k]);
		/* Assertion kills implementation at runtime of point not within prescribed box */
		pA( tmp>=bbmin && tmp<=bbmax );
		curVecIdxs.push_back(k);
	 }
	 /* the local number of sources vector stores the size of all sample points in the 0 position for the root */
	 lclSrcNumVec.push_back( curVecIdxs.size() );

	 /* Sample points max number of points per leaf box */
	 int PTSMAX = (int)ceil(ptsMax() * DENSITY);
	 
	 int level = 0;
	 int arrBeg = 0;
	 int arrEnd = 1;
	 while(arrBeg < arrEnd) {
		int arrCnt = arrEnd;
		/* In this for loop, we build each level as necessary
		* arrBeg will be set to the number of points already parsed after this for loop*/
		for(int k=arrBeg; k<arrEnd; k++) {
		  //--
		  /* lclSrcNumVec[0] is already set to the size of all of the sample positions within radius of the center
		  * So, lclSrcNumVec[k] at each level is set such that the maximum number of points is in the subsequent
		  * children boxes at a maximum level */
		  if( lclSrcNumVec[k]>PTSMAX && level<maxLevel()-1) {
			 /* the current node's child locator is set to arrCnt (the number of boxes already built)
			  * For example, when k = 0, arrCnt = 1.  The location of the root's lead child */
			 nodeVec[k].chd() = arrCnt;
			 /* Increment arrCnt by 8 */
			 arrCnt = arrCnt + pow2(dim());
			 /* Build all 8 of nodeVec[k]'s children */
			 for(int a=0; a<2; a++) {
				for(int b=0; b<2; b++) {
				  for(int c=0; c<2; c++) {
					 /* Create a new node with parent at location "k"
					  * child initially set to -1 (changed when new node is looked at)
					  * path set to 2*"k's path" + the binary index of the new node's location relative to k
					  * depth is set to "k's" depth + 1
					  */
		  			 nodeVec.push_back( Node(k,-1, 2*nodeVec[k].path2Node()+Index3(a,b,c), nodeVec[k].depth()+1) );
					 /* push a new vector of integers onto vecIdxs */
					 vecIdxs.push_back( vector<int>() );
					 /* Push an extra level onto lclSrcNumVec - set to zero for now and will be changed later */
					 lclSrcNumVec.push_back( 0 );
				  }
				}
			 }
			 /* children's vector of indices */ 
			 Point3 curctr( center(k) ); /* get center of current node */
			 /* vecIdxs[k] stores indices of sample points that node k has in it
			  * for vecIdxs[0], all of the sample points indices are stored */
			 for(vector<int>::iterator pi=vecIdxs[k].begin(); pi!=vecIdxs[k].end(); pi++) {
				/* Build a point out of the index *pi which is in the current node's vector of indices */
				Point3 tmp( allSamPosVec[*pi] );
				Index3 idx;
				/* Build an index based on the point tmp's location relative to the center of the node k */
				for(int j=0; j<dim(); j++)
				  idx(j) = (tmp(j) >= curctr(j));
				/* return the index of the child of node k located at idx */
				int chdgNodeIdx = child(k, idx);
				/* put the point references by *pi into the index of points that the node at chgNodeIdx has access to */
				vecIdxs[chdgNodeIdx].push_back(*pi);
			 }
			 vecIdxs[k].clear(); //VERY IMPORTANT
			 /* children's lclSrcNum */
			 for(int a=0; a<2; a++) {
				for(int b=0; b<2; b++) {
				  for(int c=0; c<2; c++) {
					 int chdgNodeIdx = child( k, Index3(a,b,c) );
					 /* the local number of sources for the node at chdgNodeIdx is the size of its vector of indices */
					 lclSrcNumVec[chdgNodeIdx] = vecIdxs[chdgNodeIdx].size();
				  }
				} 
			 }
		  } /* end of if */
		} /* end of for */
		/* the previous level has been built so increment the levels and update arr counters */
		level++;
		arrBeg = arrEnd;
		arrEnd = arrCnt;
	 } /* end of while */
	 
	 allSamPosVec.clear(); //SAVE SPACE
  } /* End of building the tree on processor zero */
  /* end of this scope.  everything in this scope, including lclSrcNumVec and vecIdxs no longer available */
  
  /* processor 0 sends the tree to all processors */
  {
	 int size = nodeVec.size();
	 /* Broadcast the size of NodeVec to all processes.  Then, brodcast
	  * the address of the beginning of nodeVec to all processor.
	  * More information at 
	  * http://www-unix.mcs.anl.gov/mpi/www/www3/MPI_Bcast.html
	  */
	 pC( MPI_Bcast(&size, 1, MPI_INT, 0, mpiComm()) );
	 /* Nust be resized for parallel operation, so that all nodeVec's of same size.  Does nothing for uniprocessor operation */
	 nodeVec.resize(size, Node(0,0,Index3(0,0,0),0));
	 pC( MPI_Bcast(&(nodeVec[0]), size*sizeof(Node), MPI_CHAR, 0, mpiComm()) );
  }

  /* Each processor puts its own srcdata into the tree  - new scope */
  /* local source number - number of source points in each box; different than above in differnent scope */
  vector<int> lclSrcNumVec;	 lclSrcNumVec.resize(nodeVec.size(), 0);
  /* new vecIdxs - previous one lost in previous scope.  Here vecIdxs[i] stores all indices for source positions
	* which are available to node i.  For example, vecIdxs[i].begin() to vecIdxs[i].end() store the indices of points
	* in node i or its descendants.  So, vecIdxs[0] stores the vector of indices of all points in the LET */
  vector< vector<int> > vecIdxs;	 vecIdxs.resize( nodeVec.size() );
  { 
	 /* get all source positions */
	 Vec pos = _srcPos;
	 double* posArr;	 pC( VecGetArray(pos, &posArr) );
	 /* get the number of positions for local porocessor */
	 int numLclPos = procLclNum(pos);
	 /* get the range of indices that this processor is responsible for */
	 int begLclPos, endLclPos;  procLclRan(pos, begLclPos, endLclPos);

	 /* push all local indices from this processor's range into current indices list */
	 vector<int>& curVecIdxs = vecIdxs[0];
	 Point3 bbmin(ctr()-Point3(radius()));
	 Point3 bbmax(ctr()+Point3(radius()));  
	 for(int k=begLclPos; k<endLclPos; k++) {
		Point3 tmp(posArr+(k-begLclPos)*dim());
		pA(tmp>=bbmin && tmp<=bbmax);	 /* Assert this point is within the desired radius of the center */
		curVecIdxs.push_back(k);
	 }
	 /* lclSrcNumVec[0] stores number of all indices of points available to root of tree for this processor
	  * Specifically, the number of points this processor is in charge of storing in LET */
	 lclSrcNumVec[0] = curVecIdxs.size();

	 /* ordVec - ordering of the boxes in down->up fashion */
	 vector<int> ordVec;	 pC( dwnOrderCollect(ordVec) );
	 
	 for(int i=0; i<ordVec.size(); i++) {
		int curgNodeIdx = ordVec[i]; /* current node index */
		/* store all indices put into curVecIdxs into the current Node's vector of indices */
		vector<int>& curVecIdxs = vecIdxs[curgNodeIdx];

		/* If the current node is NOT a leaf (i.e., has children in the LET)
		 * then go through current vector of indices at build its childrens'
		 * vector of indices */
		if(terminal(curgNodeIdx)==false) {
		  Point3 curctr( center(curgNodeIdx) );
		  for(vector<int>::iterator pi=curVecIdxs.begin(); pi!=curVecIdxs.end(); pi++) {
			 Point3 tmp(posArr+(*pi-begLclPos)*dim());
			 Index3 idx;
			 for(int j=0; j<dim(); j++) {
				idx(j) = (tmp(j)>=curctr(j));
			 }
			 int chdgNodeIdx = child(curgNodeIdx, idx);
			 vector<int>& chdVecIdxs = vecIdxs[chdgNodeIdx];
			 chdVecIdxs.push_back(*pi);
		  }
		  curVecIdxs.clear(); //VERY IMPORTANT
		  /* For all of the current Node's children, put the vector of indices
			* into the local source number vector */
		  for(int a=0; a<2; a++) {
			 for(int b=0; b<2; b++) {
				for(int c=0; c<2; c++) {
				  int chdgNodeIdx = child(curgNodeIdx, Index3(a,b,c));
				  lclSrcNumVec[chdgNodeIdx] = vecIdxs[chdgNodeIdx].size();
				}
			 }
		  }
		} /* end if */
	 } /* end for loop through downward ordering of the nodes in the LET */
	 pC( VecRestoreArray(pos, &posArr) );
  } /* At the end of this scope, we've built lclSrcNumVec and VecIdxs for this processor */
  
  /* decide owner of each node in LET by mpi_operator max */
  vector<int> ownVec;  ownVec.resize(nodeVec.size(), -1);
  vector<int> gsmVec;  gsmVec.resize(nodeVec.size(), 0);
  {
	 /* max to gsmVec */
	 vector<double> tmpVec;	 tmpVec.resize(nodeVec.size(), 0);
	 double extra = 1.0 - double(mpiRank()+1) / double(mpiSize());
	 /* If lclSrcNumVec for some node is zero, then set tmpVec[i] to zero
	  * Otherwise, set it to lclSrcNumVec there and add extra based on the processor rank */
	 for(int i=0; i<tmpVec.size(); i++) {
		if(lclSrcNumVec[i]==0) {
		  tmpVec[i] = 0;
		}
		else {
		  tmpVec[i] = lclSrcNumVec[i] + extra; /* extra is used to ? */
		}
	 }

	 /* maxVec will store the maximum value among processors for a specific node in the LET.
	  * That is, among N processors, if processor k owns more sources that will be allocated to
	  * Node m, then maxVec[m] = tempVec[m] = lclSrcNumVec[m] + extra for processor k */
	 vector<double> maxVec;  maxVec.resize(nodeVec.size(), 0);
	 /* Combine values from all processors and then broadcast it back out to all processors
	  * Here, take the max value among all processors from each tmpVec[i] and store it in maxVec[i]
	  * See http://www-unix.mcs.anl.gov/mpi/www/www3/MPI_Allreduce.html for more information */
	 pC( MPI_Allreduce( &(tmpVec[0]), &(maxVec[0]), nodeVec.size(), MPI_DOUBLE, MPI_MAX, mpiComm() ) );

	 /* Store in gsmVec the rank of the processor with the maximum number of local sources for a specific node */ 
	 for(int i=0; i<maxVec.size(); i++)
		if(maxVec[i]==tmpVec[i] && tmpVec[i]>0)
		  gsmVec[i] = mpiRank();
		else
		  gsmVec[i] = -1; /* gsmVec contains part of the owner information */

	 /* set owner - can be -1.  From gsmVec get which processors "owns" each node.  All processors have this information */
	 pC( MPI_Allreduce( &(gsmVec[0]), &(ownVec[0]), nodeVec.size(), MPI_INT, MPI_MAX, mpiComm() ) );
	 /* total number of each box.  Uses MPI_SUM as its operation handle.  All processors know the total number now */
	 pC( MPI_Allreduce( &(lclSrcNumVec[0]), &(gsmVec[0]), nodeVec.size(), MPI_INT, MPI_SUM, mpiComm() ) );

	 /* turn on the appropraite bits for each node */
	 for(int i=0; i<nodeVec.size(); i++)
		if(ownVec[i]==mpiRank()) /* This node is owned by this processor */
		  node(i).tag() |= LET_OWNRNODE;
	 for(int i=0; i<nodeVec.size(); i++)
		if(gsmVec[i]>0) /* This node has source points in it */
		  node(i).tag() |= LET_SRCENODE;
  } /* End of building gsmVec=number of sources for each node, ownVec=rank of node owner */


  
  /* based on the owner info, assign glbSrcNodeIdx, glbSrcExaBeg, glbSrcExaNum */
  {
	 /* numNodeVec will store how many nodes are owned by each processor */
	 vector<int> numNodeVec;	 numNodeVec.resize(mpiSize(), 0);
	 /* numNodeVec will store how many sources are owned by each processor */
	 vector<int> numExaVec;	 numExaVec.resize(mpiSize(), 0);
	 for(int i=0; i<nodeVec.size(); i++) {
		/* Assign owner to which rank owns nodeVec[i].  If none, then owner = -1 */
		int owner = ownVec[i];
		if(owner>=0) { /* Owner exists */
		  numNodeVec[owner] ++; /* Increment number of nodes processor "owner" has */
		  if(terminal(i)==true)
			 numExaVec[owner] += gsmVec[i]; /* Increment number of sources processor "owner" has */
		}
	 }	 

	 /* begNodeVec gives the beginning node index for each processor */
	 vector<int> begNodeVec;	 begNodeVec.resize(mpiSize(), 0);
	 /* begNodeVec gives the beginning source index for each processor */
	 vector<int> begExaVec;	 begExaVec.resize(mpiSize(), 0);
	 int disNode = 0; /* a displacement counter for nodes */
	 int disExa = 0; /* a displacement counter for sources */
	 /* For each processor set where that processor begins "ownership" in the indices */
	 for(int i=0; i<mpiSize(); i++) { 
		begNodeVec[i] = disNode;
		disNode += numNodeVec[i];
		begExaVec[i] = disExa;
		disExa += numExaVec[i];
		 }

	 /* For each node, the following sets a way to index the nodes and sources using parameters set in each node
	  * and available globally */
	 for(int i=0; i<nodeVec.size(); i++) {
		/* Find the owner of node(i) */
		int owner = ownVec[i];
		if(owner>=0) {
		  /* Set the global source node index for node(i) to begNodeVec[owner] and increment begNodeVec[owner] */
		  node(i).glbSrcNodeIdx() = begNodeVec[owner];
		  begNodeVec[owner] ++;
		  /* set "exact" source position variables for terminal/leaf nodes */
		  if(terminal(i)==true) {
			 /* beginning index of source positions */
			 node(i).glbSrcExaBeg() = begExaVec[owner];
			 /* exact number of sources */
			 node(i).glbSrcExaNum() = gsmVec[i]; 
			 begExaVec[owner] += gsmVec[i];
		  }	
		}
		/* If this node has no owner, then there are no sources in it */
		else { 
		  node(i).glbSrcNodeIdx() = -1; 
		}
	 }

	 /* _glbSrcNodeCnt is the total number of nodes "owned" by a processor (i.e., nodes with sources in them where ownVec != -1 */
	 _glbGlbSrcNodeCnt = 0;
	 for(int i=0; i<ownVec.size(); i++){
		if(ownVec[i]>=0)
		  _glbGlbSrcNodeCnt++;
	 }
	 /* global count of exact number of sources is just the size of _srcPos */
	 _glbGlbSrcExaCnt = procGlbNum(_srcPos);
	 /* local number of nodes from the global count is the the number of nodes for this processor */
	 _lclGlbSrcNodeCnt = numNodeVec[mpiRank()];
	 /* local number of sources from the global count is the the number of sources for this processor */
	 _lclGlbSrcExaCnt = numExaVec[mpiRank()];

	 /* ownVec and gsmVec no longer needed, so they are cleared */
	 ownVec.clear();
	 gsmVec.clear();
  }
  
  /* mpi_scan local number info to accumulative number info and write ctbSrcNodeIdx, ctbSrcExaBeg, ctbSrcExanum */
  {
	 vector<int> accVec;	 accVec.resize(nodeVec.size(), 0);
	 /* The following MPI_Scan puts in accVec[i] the total sum of lclSrcNumVec[i] ( i is a node in the LET )
	  * of all previous processors.  For example, if lclSrcNumVec[i] = 7 on processor 1,
	  *                                              lclSrcNumVec[i] = 3 on processor 2,
	  *                                              lclSrcNumVec[i] = 8 on processor 3,
	  * then, MPI_Scan computes the following:       accVec[i] = 7 on processor 1,
	  *                                              accVec[i] = 10 on processor 2,
	  *                                              accVec[i] = 18 on processor 3
	  * For more information on MPI_Scan:
	  * http://www-unix.mcs.anl.gov/mpi/www/www3/MPI_Scan.html
	  */
	 pC( MPI_Scan(&(lclSrcNumVec[0]), &(accVec[0]), nodeVec.size(), MPI_INT, MPI_SUM, mpiComm()) );
	 /* For each node on each processor, subtract lclSrcNumVec[i] from accVec[i], such that
	  * accVec[i] contains the number of local sources for node i only on processors
	  * of lower rank */
	 for(int i=0; i<accVec.size(); i++){
		accVec[i] -= lclSrcNumVec[i];
	 }

	 /* At end, will give total count of all nodes which have sources which contribute for this processor */
	 _ctbSrcNodeCnt = 0;
	 /* At end, will give total count of all source positions which contribute for this processor */
	 _ctbSrcExaCnt = 0; 
	 for(int i=0; i<nodeVec.size(); i++) {
		if(lclSrcNumVec[i]>0) {
		  /* Turn on Contributor bit - this processor contributes to this node */
		  node(i).tag() |= LET_CBTRNODE;
		  /* Number of additional nodes (ancestors) thus far which contribute to this node with non-zero source numbers */
		  node(i).ctbSrcNodeIdx() = _ctbSrcNodeCnt; 
		  _ctbSrcNodeCnt ++;
		  /* push map for contributor node to global source node index */
		  _ctb2GlbSrcNodeMap.push_back( node(i).glbSrcNodeIdx() );

		  /* For leaves set beginning index and exact number for contributor source positions */
		  if(terminal(i)==true) { 
			 node(i).ctbSrcExaBeg() = _ctbSrcExaCnt;
			 node(i).ctbSrcExaNum() = lclSrcNumVec[i];
			 /* Accumulate sum of the exact number of positions that contribute for this processor */
			 _ctbSrcExaCnt += lclSrcNumVec[i];
			 /* Set node(i)'s list of contributing source positions for this processor */
			 node(i).ctbSrcOwnVecIdxs() = vecIdxs[i];
			 /* Create a mapping from contributor source positions to global source position indices */
			 for(int k=0; k < vecIdxs[i].size(); k++){
				_ctb2GlbSrcExaMap.push_back( node(i).glbSrcExaBeg()+accVec[i]+k );
			 }
		  }
		}
	 }
  }  
  return(0);
}

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "Let3d_MPI::trgData"
int Let3d_MPI::trgData()
{
  /* each proc put its own target data into the tree */
  vector<Node>& nodeVec = this->nodeVec();
  vector<int> lclSrcNumVec;  lclSrcNumVec.resize(nodeVec.size(), 0);
  vector< vector<int> > vecIdxs;  vecIdxs.resize(nodeVec.size());
  {
	 Vec pos = _trgPos;
	 double* posArr;	 pC( VecGetArray(pos, &posArr) );
	 /* Get the number of positions for this processor */
	 int numLclPos = procLclNum(pos);
	 /* Get the beginning and ending range values this processor is responsible for */
	 int begLclPos, endLclPos;  procLclRan(pos, begLclPos, endLclPos);
	 
	 /* Create vector of indices of points for root of the LET */
	 vector<int>& curVecIdxs = vecIdxs[0];
	 Point3 bbmin(ctr()-Point3(radius()));
	 Point3 bbmax(ctr()+Point3(radius())); 
	 for(int k=begLclPos; k<endLclPos; k++) {
		Point3 tmp(posArr+(k-begLclPos)*dim());  
		pA(tmp>=bbmin && tmp<=bbmax);	 //LEXING: IMPORTANT : Asserts each point is within range of the center of the root
		curVecIdxs.push_back(k);
	 }
	 /* local number of sources for the root */
	 lclSrcNumVec[0] = curVecIdxs.size();

	 /* ordVec - ordering of the boxes in down->up fashion */
	 vector<int> ordVec;	 pC( dwnOrderCollect(ordVec) );
	 for(int i=0; i<ordVec.size(); i++) {
		int curgNodeIdx = ordVec[i]; /* current node index */
		vector<int>& curVecIdxs = vecIdxs[curgNodeIdx];

		/* Do the following for non-leaf nodes */
		if(terminal(curgNodeIdx)==false) {
		  Point3 curctr( center(curgNodeIdx) );
		  for(vector<int>::iterator pi=curVecIdxs.begin(); pi!=curVecIdxs.end(); pi++) {
			 /* Construct a point from the pointer as required for the current vector of indices
			  * Begins for root which is pre-built from above, and then all children indices of points will be built
			  * for use in as we traverse the ordVec list */
			 Point3 tmp(posArr+(*pi-begLclPos)*dim());
			 Index3 idx;
			 for(int j=0; j<dim(); j++) {
				idx(j) = (tmp(j)>=curctr(j));
			 }
			 /* Retrieve the child of the current node based on its location relative to the point tmp
			 * in order to find where the point *pi resides (which child node) and then push that point
			 * into the child's vector of indices
			 */
			 int chdgNodeIdx = child(curgNodeIdx, idx);
			 vector<int>& chdVecIdxs = vecIdxs[chdgNodeIdx];
			 chdVecIdxs.push_back(*pi);
		  } /* End of for loop for curVecIdxs */
		  curVecIdxs.clear(); /* VERY IMPORTANT */
		  /* For all 8 children of the current node, set the local number of sources based
			* on the size of its vector of indices */
		  for(int a=0; a<2; a++) {
			 for(int b=0; b<2; b++) {
				for(int c=0; c<2; c++) {
				  int chdgNodeIdx = child(curgNodeIdx, Index3(a,b,c));
				  lclSrcNumVec[chdgNodeIdx] = vecIdxs[chdgNodeIdx].size();
				}
			 }
		  }
		}
	 } /* End of for loop for ordVec traversal */
	 pC( VecRestoreArray(pos, &posArr) );
  } /* End of this scope:  vecIdxs and lclSrcNumVec built */
  
  /* write evaTrgNodeIdx, evaTrgExaBeg, evaTrgExaNum
	* Here, for all nodes, we look to see if the local processor actually stores
	* any local sources there such that we know if the local processor is responsible
	* for evaluations at these nodes for these points */
  {
	 _evaTrgNodeCnt = 0;
	 _evaTrgExaCnt = 0;
	 for(int i=0; i<nodeVec.size(); i++) {
		if(lclSrcNumVec[i]>0) {
		  node(i).tag() |= LET_EVTRNODE; //LEXING - Turn on evaluator built
		  node(i).evaTrgNodeIdx() = _evaTrgNodeCnt;
		  _evaTrgNodeCnt ++;
		  if(terminal(i)==true) {
			 node(i).evaTrgExaBeg() = _evaTrgExaCnt;
			 node(i).evaTrgExaNum() = lclSrcNumVec[i];
			 _evaTrgExaCnt += lclSrcNumVec[i];
			 node(i).evaTrgOwnVecIdxs() = vecIdxs[i];
		  }
		}
	 }
  }
  
  /* based on the owner info write usrSrcNodeIdx, usrSrcExaBeg, usrSrcExaNum
	* We will build the U,V,W, and X lists for each node in the LET, such that any box B
	* in these lists is used by the current processor.  That is, processor P is a "user" of B  */
  {
	 for(int i=0; i<nodeVec.size(); i++){
		if(lclSrcNumVec[i]>0) {
		  /* make sure that U,V,W, and X lists can be calculated/built properly.  See calGlbNodeLists(Node i) for more information */
		  pC( calGlbNodeLists(i) );
		  /* For all nodes in U,V,W, and X lists, turn on the bit to indicate that this processor
			* uses these nodes for computation */
		  for(vector<int>::iterator vi=node(i).Unodes().begin(); vi!=node(i).Unodes().end(); vi++)
			 node(*vi).tag() |= LET_USERNODE;
		  for(vector<int>::iterator vi=node(i).Vnodes().begin(); vi!=node(i).Vnodes().end(); vi++)
			 node(*vi).tag() |= LET_USERNODE;
		  for(vector<int>::iterator vi=node(i).Wnodes().begin(); vi!=node(i).Wnodes().end(); vi++)
			 node(*vi).tag() |= LET_USERNODE;
		  for(vector<int>::iterator vi=node(i).Xnodes().begin(); vi!=node(i).Xnodes().end(); vi++)
			 node(*vi).tag() |= LET_USERNODE;
		}
	 }

	 /* Count the number of nodes being used by the current processor
	  * and build a map between how these nodes are number and the global node indices.
	  * Do the same for the source positions available to this processor
	  */
	 _usrSrcNodeCnt = 0;
	 _usrSrcExaCnt = 0;
	 for(int i=0; i<nodeVec.size(); i++) {
		if(node(i).tag() & LET_USERNODE) {
		  node(i).usrSrcNodeIdx() = _usrSrcNodeCnt;
		  _usrSrcNodeCnt ++;
		  _usr2GlbSrcNodeMap.push_back( node(i).glbSrcNodeIdx() );
		  if(terminal(i)==true) {
			 node(i).usrSrcExaBeg() = _usrSrcExaCnt;
			 node(i).usrSrcExaNum() = node(i).glbSrcExaNum();
			 _usrSrcExaCnt += node(i).glbSrcExaNum();
			 for(int k=0; k<node(i).glbSrcExaNum(); k++){
				_usr2GlbSrcExaMap.push_back( node(i).glbSrcExaBeg()+k );
			 }
		  }
		}
	 }
  } 
  
  return(0);
} /* end of trgData()

// ----------------------------------------------------------------------
/* Build/Calculate the global node lists (U,V,W, and X lists) for a specific node */
#undef __FUNCT__
#define __FUNCT__ "Let3d_MPI::calGlbNodeLists"
int Let3d_MPI::calGlbNodeLists(int gNodeIdx)
{
  //begin
  /* We use sets here since the values will be unnecessary to keep associated with the keys.
	* Can also use a C++ map, but it is unnecessary here to maintain unique key/value pairs
	* See let3d_mpi.hpp for descriptions of what U,B,W, and X lsists are. */
  set<int> Uset, Vset, Wset, Xset;  
  int curgNodeIdx = gNodeIdx;
  /* Make sure current node is not the root as the root has no ancestors */
  if(root(curgNodeIdx)==false) {
	 /* get parent of node we're interested in */
	 int pargNodeIdx = parent(curgNodeIdx);
	 Index3 minIdx(0); /* = (0,0,0) */
	 Index3 maxIdx(pow2(depth(curgNodeIdx))); /* = (2^d, 2^d, 2^d) */
	 
	 /* Try several different paths.  Here we are looking for paths which
	  * will lead to children of neighbors of the parent of the current node
	  * in order to build U,V,W, and X lists */
	 for(int i=-2; i<4; i++) {
		for(int j=-2; j<4; j++) {
		  for(int k=-2; k<4; k++) {
			 /* New Path to try */
			 Index3 tryPath2Node( 2*path2Node(pargNodeIdx) + Index3(i,j,k) );
			 /* Verify new path does not exceed a maximu index and is greater than (0,0,0) and is not equal to current node's path */
			 if(tryPath2Node >= minIdx && tryPath2Node <  maxIdx && tryPath2Node != path2Node(curgNodeIdx)) {
				/* Look for the children of the neighbors of the current node's parent */
				int resgNodeIdx = findGlbNode(depth(curgNodeIdx), tryPath2Node);
				/* adj = true if nodes have at least one edge touching */
				bool adj = adjacent(resgNodeIdx, curgNodeIdx);
				/* When test node is higher in the tree than current node */
				if( depth(resgNodeIdx)<depth(curgNodeIdx) ) {
				  /* If adjacent and current node is a leaf, put test node in the Uset if the current node is a leaf */
				  if(adj){
					 if(terminal(curgNodeIdx)){ /* need to test if res a leaf?? HARPER */
						Uset.insert(resgNodeIdx);
					 }
					 else {;} /* non-leaves do not have U-lists */
				  }
				  else{
					 /* The nodes are not adjacent, but resgNodeIdx is still a child of the current node's parent's neighbors.
					  * Hence, the current node's parent is adjacent to resgNodeIdx.  IN general. the current node is in the
					  * W-list of resgNodeIdx */
					 Xset.insert(resgNodeIdx);
				  }
				} /* End of "if depth(res) < depth(current) " */
				
				/* Current node and test node at same depth */
				if( depth(resgNodeIdx)==depth(curgNodeIdx) ) {
				  /* Two nodes at same depth not adjacent */
				  if(!adj) {
					 Index3 bb(path2Node(resgNodeIdx)-path2Node(curgNodeIdx));
					 /* Verify that no single component of the two paths exceed a limit */
					 assert( bb.linfty()<=3 );
					 /* resgNodeIdx is a child of the neighbor's of the currentr node's parents and not-adjacent to current node.
					  * Hence, it is in the V-list */
					 Vset.insert(resgNodeIdx);
				  }
				  /* nodes at same depth and are adjacent, so resgNodeIdx could be in U OR W lists */
				  else {
					 if(terminal(curgNodeIdx)) {
						/* queue:  elements added/pushed to the back and removed/popped from the front */ 
						queue<int> rest;
						/* push resgNodeIdx into the queue */
						rest.push(resgNodeIdx);
						while(rest.empty()==false) {
						  int fntgNodeIdx = rest.front(); rest.pop(); /* Set front temp node index and pop the list */
						  /* If the current node and temp node are not adjacent, temp node is in W-list of current node */
						  if(adjacent(fntgNodeIdx, curgNodeIdx)==false) {
							 Wset.insert( fntgNodeIdx );
						  }
						  /* Current node and temp node ARE adjacent */
						  else {
							 /* If temp node is a leaf, it is in current node's U-list */
							 if(terminal(fntgNodeIdx)) {
								Uset.insert(fntgNodeIdx);
							 }
							 /* If temp node is not a leaf, then one of its descendants may be a leaf in the U or W lists of current node
							 * So, we push those into the queue
							 */
							 else { 
								for(int a=0; a<2; a++) {
								  for(int b=0; b<2; b++) {
									 for(int c=0; c<2; c++) {
										rest.push( child(fntgNodeIdx, Index3(a,b,c)) );
									 }
								  }
								}
							 }
						  }
						} /* End of while loop for rest queue */
					 } /* End of "if current node a leaf" */ 
				  } /* End of res and current nodes at same depth and adjacent */
				} /* End of "if depth(res) == depth(current)" */
			 } /* End of trypath */
		  }
		}
	 }
  } /* End of building sets for non-root node, curgNodeIdx */

  /* If the current node is a leaf, then it is in its own U-list */
  if(terminal(curgNodeIdx))
	 Uset.insert(curgNodeIdx);

  /* For all sets, check to make sure all of the nodes actually have sources in them.  If not
	* we do not need to put them in the lists.  If so, build the U,V,W, and X lists from respective sets
	*/
  for(set<int>::iterator si=Uset.begin(); si!=Uset.end(); si++)
	 if(node(*si).tag() & LET_SRCENODE) /* Check if sources bit set for this node's tag.  Same below */
		node(gNodeIdx).Unodes().push_back(*si);
  for(set<int>::iterator si=Vset.begin(); si!=Vset.end(); si++)
	 if(node(*si).tag() & LET_SRCENODE)
		node(gNodeIdx).Vnodes().push_back(*si);
  for(set<int>::iterator si=Wset.begin(); si!=Wset.end(); si++)
	 if(node(*si).tag() & LET_SRCENODE)
		node(gNodeIdx).Wnodes().push_back(*si);
  for(set<int>::iterator si=Xset.begin(); si!=Xset.end(); si++)
	 if(node(*si).tag() & LET_SRCENODE)
		node(gNodeIdx).Xnodes().push_back(*si);
  
  return(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "Let3d_MPI::dwnOrderCollect"
int Let3d_MPI::dwnOrderCollect(vector<int>& ordVec)
{
  /* ordVec - ordering of the boxes in downward fashion */
  ordVec.clear();
  for(int i=0; i<_nodeVec.size(); i++)
	 ordVec.push_back(i);
  return(0);
}
// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "Let3d_MPI::upwOrderCollect"
int Let3d_MPI::upwOrderCollect(vector<int>& ordVec)
{
  /* ordVec - ordering of the boxes in upward fashion */
  ordVec.clear();
  for(int i=_nodeVec.size()-1; i>=0; i--)
	 ordVec.push_back(i);
  return(0);
}
// ----------------------------------------------------------------------
/* Return one of the eight children of a node based on a binary index, idx.  For example,
 * if idx = (0,0,0), the first child is returned.
 * If idx = (0,0,1), the second child is returned.
 * If idx = (0,1,0), the third child is returned.
 * If idx = (0,1,1), the fourth child is returned.
 * ...
 * If idx = (1,1,1), the eighth child is returned.
 */
int Let3d_MPI::child(int gNodeIdx, const Index3& idx)
{
  assert(idx>=Index3(0) && idx<Index3(2));
  return node(gNodeIdx).chd() + (idx(0)*4+idx(1)*2+idx(2));
}
/* Construct the center for a specific node */
Point3 Let3d_MPI::center(int gNodeIdx) 
{
  Point3 ll( ctr() - Point3(radius()) );
  int tmp = pow2(depth(gNodeIdx));
  Index3 path2Node_Lcl(path2Node(gNodeIdx));
  Point3 res;
  for(int d=0; d<dim(); d++) {
	 res(d) = ll(d) + (2*radius()) * (path2Node_Lcl(d)+0.5) / double(tmp);
  }
  return res;
}
/* Radius of a node =
 * radius(root)/(2^depth(node)).
 * When radius of root is 1 (most often),
 * radius of a node is 2^(-d)
 */
double Let3d_MPI::radius(int gNodeIdx) //radius of a node
{
  return radius()/double(pow2(depth(gNodeIdx)));
}
// ----------------------------------------------------------------------
/* Find a node based on a depth and some path.
 * This function is used to try to find children of the neighbors of the parent of some
 * node at depth wntDepth by trying several paths
 */
int Let3d_MPI::findGlbNode(int wntDepth, const Index3& wntPath2Node)
{
  int cur = 0;  
  Index3 tmpPath2Node(wntPath2Node);
  /* Keep trying to find a node cur which has depth greater than/equal to the given depth
	* OR is a terminal/leaf */
  while(depth(cur)<wntDepth && terminal(cur)==false) {
	 /* Difference in the depths */
	 int dif = wntDepth-depth(cur);
	 /* 2^(dif-1):  2^(difference in depths > 0) */
	 int tmp = pow2(dif-1);
	 /* Returns a binary index because of tmp's size */
	 Index3 choice( tmpPath2Node/tmp );

	 /* Get new path to follow */
	 tmpPath2Node -= choice*tmp;
	 /* Choose new node - hopefully this node is at the same depth as our current node */
	 cur = child(cur, choice);	 
  }  //cerr<<endl;
  return cur;
}
// ----------------------------------------------------------------------
/* Returns true if two nodes have an edge or face which touch.  False otherwise */
bool Let3d_MPI::adjacent(int me, int you)
{
  int maxDepth = max(depth(me),depth(you)); /* Which node is at a greater depth */
  Index3 one(1); /* = (1,1,1) */
  /* Construct the center for both nodes */
  Index3 meCtr(  (2*path2Node(me)+one) * pow2(maxDepth - depth(me))  );
  Index3 youCtr(  (2*path2Node(you)+one) * pow2(maxDepth - depth(you))  );
  /* construct radii for both nodes */
  int meRad = pow2(maxDepth - depth(me));
  int youRad = pow2(maxDepth - depth(you));
  /* absolute value of difference of centers */
  Index3 dif( abs(meCtr-youCtr) );
  /* sum of radii */
  int radius  = meRad + youRad;
  /* if the abs. value of difference of centers is less than the sum of the radii
	* AND the infinity norm equals the sum of the radii, then the two nodes
	* are not too far away AND at least one edge touches */
  return
	 ( dif <= Index3(radius) ) && //not too far
	 ( dif.linfty() == radius ); //at least one edge touch
}

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "Let3d_MPI::print"
int Let3d_MPI::print()
{
  //begin  //pC( PetscPrintf(MPI_COMM_WORLD, "nodeVec %d\n", _nodeVec.size()) );
  pC( PetscPrintf(MPI_COMM_SELF, "%d: %d glbGlbSrc %d %d  lclGlbSrc %d %d ctbSrc %d %d evaTrg %d %d usrSrc %d %d\n",
						mpiRank(), _nodeVec.size(),
						_glbGlbSrcNodeCnt, _glbGlbSrcExaCnt,   _lclGlbSrcNodeCnt, _lclGlbSrcExaCnt,
						_ctbSrcNodeCnt, _ctbSrcExaCnt,	 _evaTrgNodeCnt, _evaTrgExaCnt,			_usrSrcNodeCnt, _usrSrcExaCnt) );
  return(0);
}


