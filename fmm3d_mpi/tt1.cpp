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

using namespace std;

int main(int argc, char** argv)
{
  PetscInitialize(&argc,&argv,"options_0",NULL);  //PetscTruth flg = PETSC_FALSE;

  MPI_Comm comm; comm = PETSC_COMM_WORLD;
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  int dim = 3;
  srand48( (long)( time(NULL)+mpirank ) );
  
  PetscTruth flg = PETSC_FALSE;
  pC( PetscPrintf(MPI_COMM_WORLD, "mpisize %d\n", mpisize) );
  
  //1. random src and trg data
  int numsrc;  pC( PetscOptionsGetInt(PETSC_NULL, "-numsrc", &numsrc, &flg) ); pA(flg==PETSC_TRUE);
  int numtrg;  pC( PetscOptionsGetInt(PETSC_NULL, "-numtrg", &numtrg, &flg) ); pA(flg==PETSC_TRUE);
  int kt;  pC( PetscOptionsGetInt(PETSC_NULL, "-kt", &kt, &flg) ); pA(flg==PETSC_TRUE);
  
  vector<double> tmp(2);	 tmp[0] = 1;	 tmp[1] = 0.25; //coefs in the kernel, work for all examples
  Kernel3d_MPI knl(kt, tmp);
  
  int lclnumsrc = (numsrc+mpirank)/mpisize;
  int lclnumtrg = (numtrg+mpirank)/mpisize;
  
  //pick center
  double R = 0.85;
  double r = 1-R;
  double lclctr[3];
  for(int d=0; d<dim; d++)	 lclctr[d] = R * (2.0*drand48()-1.0);
  
  //generate src and trg data
  Vec srcPos;  pC( VecCreateMPI(comm, lclnumsrc*dim, PETSC_DETERMINE, &srcPos) );
  double* srcPosarr; pC( VecGetArray(srcPos, &srcPosarr) );
  for(int k=0; k<lclnumsrc; k++)
	 for(int d=0; d<dim; d++)
		srcPosarr[d+dim*k] = lclctr[d] + r*(2.0*drand48()-1.0);
  pC( VecRestoreArray(srcPos, &srcPosarr) );
  
  Vec trgPos;  pC( VecCreateMPI(comm, lclnumtrg*dim, PETSC_DETERMINE, &trgPos) );
  double* trgPosarr; pC( VecGetArray(trgPos, &trgPosarr) );
  for(int k=0; k<lclnumtrg; k++)
	 for(int d=0; d<dim; d++)
		trgPosarr[d+dim*k] = lclctr[d] + r*(2.0*drand48()-1.0);
  pC( VecRestoreArray(trgPos, &trgPosarr) );
  
  int srcDOF = knl.srcDOF();
  int trgDOF = knl.trgDOF();
  Vec srcDen;  pC( VecCreateMPI(comm, lclnumsrc*srcDOF, PETSC_DETERMINE, &srcDen) );
  double* srcDenarr; pC( VecGetArray(srcDen, &srcDenarr) );
  for(int k=0; k<lclnumsrc*srcDOF; k++)
	 srcDenarr[k] = drand48();
  pC( VecRestoreArray(srcDen, &srcDenarr) );
  
  Vec trgVal;  pC( VecCreateMPI(comm, lclnumtrg*trgDOF, PETSC_DETERMINE, &trgVal) );
  double zero = 0;  pC( VecSet(trgVal, zero) );
  
  //2. allocate fmm 
  clock_t ck0, ck1;
  
  FMM3d_MPI* fmm = new FMM3d_MPI("fmm3d_MPI_");
  fmm->srcPos()=srcPos;  fmm->srcNor()=srcPos;
  fmm->trgPos()=trgPos;
  fmm->ctr() = Point3(0,0,0); // CENTER OF THE TOPLEVEL BOX
  fmm->rootLevel() = 0;         // 2^(-rootlvl) is the RADIUS OF THE TOPLEVEL BOX
  fmm->knl() = knl;
  //setup
  ck0 = clock();
  pC( fmm->setup() );
  ck1 = clock();  
  pC( PetscPrintf(MPI_COMM_WORLD, "fmm setup used %e secs\n", double(ck1-ck0)/CLOCKS_PER_SEC) );
  
  //3. run fmm
  for(int i=0; i<3; i++) {
	 ck0 = clock();
	 pC( fmm->evaluate(srcDen, trgVal) );
	 //pC( MPI_Barrier(MPI_COMM_WORLD) );  //check Vincnt, Votcnt
	 ck1 = clock();
	 pC( PetscPrintf(MPI_COMM_WORLD, "fmm evaluate used %e secs\n", double(ck1-ck0)/CLOCKS_PER_SEC) );
  }
  
  //4. check
  ck0 = clock();
  double rerr;
  pC( fmm->check(srcDen, trgVal, 20, rerr) );
  pC( PetscPrintf(MPI_COMM_WORLD, "relative %e\n", rerr) );
  ck1 = clock();
  pC( PetscPrintf(MPI_COMM_WORLD, "fmm check used %e secs\n", double(ck1-ck0)/CLOCKS_PER_SEC) );
  //extra
  //{ PetscViewer viewer;  pC( PetscViewerASCIIOpen(comm, "srcPos.out", &viewer) );  pC( VecView(srcPos, viewer) );  pC( PetscViewerDestroy(viewer) ); }
  //{ PetscViewer viewer;  pC( PetscViewerASCIIOpen(comm, "trgPos.out", &viewer) );  pC( VecView(trgPos, viewer) );  pC( PetscViewerDestroy(viewer) ); }
  
  delete fmm;
  
  pC( VecDestroy(srcPos) );
  pC( VecDestroy(trgPos) );
  pC( VecDestroy(srcDen) );
  pC( VecDestroy(trgVal) );
  
  PetscFinalize();
  return 0;
}
