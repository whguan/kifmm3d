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

/** \file
 * This file implements the check function for FMM3d_MPI.
 * Fmm3d_MPI::check allows for the checking of the relative error to make sure it meets the pre-specified deisred amount.
 */


/* ********************************************************************** */
/* Check takes the source densities, target values as solved for and a value
 * numChk, which is the point being checked.  The value rerr (relative error)
 * is returned.
 */
#undef __FUNCT__
#define __FUNCT__ "FMM3d_MPI::check"
int FMM3d_MPI::check(Vec srcDen, Vec trgVal, int numChk, double& rerr)
{
  //begin
  pA(procLclNum(_trgPos)>0); //have things to check
  //pC( PetscPrintf(mpiComm(), "checking.........\n") );
  int dim = this->dim();
  int srcDOF = this->srcDOF();
  int trgDOF = this->trgDOF();
  int slclnum = procLclNum(_srcPos);
  int tlclnum = procLclNum(_trgPos);   //int lclnum = this->lclnum();
  
  //1. select point to check
  int ttlChk = numChk*mpiSize();
  vector<int> chkVec(numChk);
  for(int k=0; k<numChk; k++) {
	 chkVec[k] = int( floor(drand48()*tlclnum) );
	 pA(chkVec[k]>=0 && chkVec[k]<tlclnum);
  }
  pC( MPI_Barrier(mpiComm()) );
  
  Vec chkPos; pC( VecCreateMPI(mpiComm(), numChk*dim,  ttlChk*dim,  &chkPos) ); //check position
  Vec chkVal; pC( VecCreateMPI(mpiComm(), numChk*trgDOF, ttlChk*trgDOF, &chkVal) ); //check value
  Vec chkDen; pC( VecCreateMPI(mpiComm(), numChk*trgDOF, ttlChk*trgDOF, &chkDen) ); //check density
  
  double* posArr; pC( VecGetArray(_trgPos, &posArr) );
  double* valArr; pC( VecGetArray( trgVal, &valArr) );
  double* chkPosArr; pC( VecGetArray(chkPos, &chkPosArr) );
  double* chkValArr; pC( VecGetArray(chkVal, &chkValArr) );
  for(int k=0; k<numChk; k++) {
	 for(int i=0; i<dim; i++)
		chkPosArr[ k*dim+i  ] = posArr[ chkVec[k]*dim+i  ];
	 for(int i=0; i<trgDOF;i++)
		chkValArr[ k*trgDOF+i ] = valArr[ chkVec[k]*trgDOF+i ];
  } 
  pC( VecRestoreArray(_trgPos, &posArr) );
  pC( VecRestoreArray( trgVal, &valArr) );
  pC( VecRestoreArray(chkPos, &chkPosArr) );
  pC( VecRestoreArray(chkVal, &chkValArr) );
  
  //2. compute and gather
  Vec allChkPos;
  {
	 VecScatter ctx;  pC( VecScatterCreateToAll(chkPos, &ctx, &allChkPos) );
	 pC( VecScatterBegin(chkPos, allChkPos, INSERT_VALUES, SCATTER_FORWARD, ctx) );
	 pC( VecScatterEnd(  chkPos, allChkPos, INSERT_VALUES, SCATTER_FORWARD, ctx) );
	 pC( VecScatterDestroy(ctx) );
  }
  Vec allChkDen; pC( VecCreateSeq(PETSC_COMM_SELF, ttlChk*trgDOF, &allChkDen) );
  Vec glbChkDen; pC( VecCreateSeq(PETSC_COMM_SELF, ttlChk*trgDOF, &glbChkDen) );
  
  double* lclSrcPosArr; pC( VecGetArray(_srcPos, &lclSrcPosArr) );
  double* lclSrcNorArr; pC( VecGetArray(_srcNor, &lclSrcNorArr) );
  double* lclSrcDenArr; pC( VecGetArray( srcDen, &lclSrcDenArr) );
  double* allChkPosArr; pC( VecGetArray(allChkPos, &allChkPosArr) );
  double* allChkDenArr; pC( VecGetArray(allChkDen, &allChkDenArr) );
  double* glbChkDenArr; pC( VecGetArray(glbChkDen, &glbChkDenArr) );
  
  DblNumMat lclSrcPosMat(dim, slclnum, false, lclSrcPosArr);
  DblNumMat lclSrcNorMat(dim, slclnum, false, lclSrcNorArr);
  DblNumVec lclSrcDenVec(srcDOF*slclnum, false, lclSrcDenArr);
  DblNumMat allChkPosMat(dim,  ttlChk, false, allChkPosArr);
  DblNumVec allChkDenVec(trgDOF* ttlChk, false, allChkDenArr);
  
  DblNumMat inter(trgDOF, slclnum*srcDOF);
  for(int i=0; i<ttlChk; i++) {
	 DblNumMat oneChkPosMat(dim, 1, false, allChkPosMat.clmdata(i));
	 DblNumVec oneChkDenVec(trgDOF, false, allChkDenVec.data()+i*trgDOF);
	 pC( _knl.buildKnlIntCtx(lclSrcPosMat, lclSrcNorMat, oneChkPosMat, inter) );
	 pC( dgemv(1.0, inter, lclSrcDenVec, 0.0, oneChkDenVec) );
  }
  
  pC( MPI_Barrier(mpiComm()) );
  pC( MPI_Allreduce(allChkDenArr, glbChkDenArr, ttlChk*trgDOF, MPI_DOUBLE, MPI_SUM, mpiComm()) );
  
  //3. distribute to individual
  double* lclChkDenArr; pC( VecGetArray(chkDen, &lclChkDenArr) );
  void* dst = (void*)lclChkDenArr;
  void* src = (void*)(glbChkDenArr + mpiRank()*numChk*trgDOF);  
  memcpy(dst, src, sizeof(double)*numChk*trgDOF);
  pC( VecRestoreArray(chkDen, &lclChkDenArr) );
  
  pC( VecRestoreArray(_srcPos, &lclSrcPosArr) );
  pC( VecRestoreArray(_srcNor, &lclSrcNorArr) );
  pC( VecRestoreArray( srcDen, &lclSrcDenArr) );
  pC( VecRestoreArray(allChkPos, &allChkPosArr) );
  pC( VecRestoreArray(allChkDen, &allChkDenArr) );
  pC( VecRestoreArray(glbChkDen, &glbChkDenArr) );
  
  pC( VecDestroy(allChkPos) );
  pC( VecDestroy(allChkDen) );
  pC( VecDestroy(glbChkDen) );
  
  PetscScalar mone = -1.0;
  pC( VecAXPY(chkVal, mone, chkDen) ); //chkVal now is the error - HARPER:  In correct order??
  double nerr, npot;
  pC( VecNorm(chkVal, NORM_2, &nerr) );
  pC( VecNorm(chkDen, NORM_2, &npot) );
  
  rerr = nerr/npot;//pC( PetscPrintf(mpiComm(), "Error %e, %e %e\n", nerr/npot, nerr, npot) );
  
  pC( VecDestroy(chkPos) );
  pC( VecDestroy(chkVal) );
  pC( VecDestroy(chkDen) );
  
  return(0);
}



