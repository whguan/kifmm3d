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
#ifndef _MATMGNT3D_HPP_
#define _MATMGNT3D_HPP_

#include "common/nummat.hpp"
#include "common/numtns.hpp"
#include "common/offtns.hpp"
#include "common/vec3t.hpp"
#include "kernel3d.hpp"
#include "comobject.hpp"

using std::map;
using std::pair;

//--------------------------------------
//unique identifier: equation
class MatMgnt3d
{
public:
  /*! UE = Upper Equivalent
	*  UC = Upper Check
	*  DE = Downward Equivalent
	*  DC = Downward Check
	*/
  enum {	 UE=0,	 UC=1,	 DE=2,	 DC=3,  };
public:
  //PARAMS(REQ) -- has to be set by parent class
  Kernel3d _knl; //the elq used by matmagnt (provide from fmm)
  int _np;
  bool _hom;
  vector<double> _degVec;
  //COMPONENTS

  /*! Upward Check To Upward Equivalent */
  map<int, DblNumMat> _upwChk2UpwEqu;
  /*! Upward Equivalent To Upward Check */
  map<int, NumTns<DblNumMat> > _upwEqu2UpwChk;
  /*! Downward Check To Downward Equivalent */
  map<int, DblNumMat> _dwnChk2DwnEqu;
  /*! Downward Equivalent To Downward Check */
  map<int, NumTns<DblNumMat> > _dwnEqu2DwnChk;
  /*! Upward Equivalent To Downward Check */
  map<int, OffTns<DblNumMat> > _upwEqu2DwnChk;
  /*! sample positions - different depending on whether UW, UC, DE or DC*/
  DblNumMat _samPos[4];
  /*! regular positions */
  DblNumMat _regPos; 
  
public:
  MatMgnt3d();
  ~MatMgnt3d();
  //MEMBER ACCESS
  Kernel3d& knl() { return _knl; }
  int& np() { return _np; }
  double alt(); //TODO: decide it based on np
  //...
  /*! src degree of freedom */
  int srcDOF() { return _knl.srcDOF(); }  //int tdof() { return eq().tdof(qt()); }
  /*! trg degree of freedom */
  int trgDOF() { return _knl.trgDOF(); }
  /*! dimension = 3 */
  int dim() { return 3; }
  //SETUP AND USE
  int setup();
  int report();
  /*! the size of plainn data */
  int plnDatSze(int tp); 
  /*! the size of eff data (data stored when using FFT) */
  int effDatSze(int tp); 
  
  int UpwChk2UpwEqu_dgemv (int level,             const DblNumVec&, DblNumVec&);
  int UpwEqu2UpwChk_dgemv (int level, Index3 ii,  const DblNumVec&, DblNumVec&);
  int DwnChk2DwnEqu_dgemv (int level,             const DblNumVec&, DblNumVec&);
  int DwnEqu2DwnChk_dgemv (int level, Index3 ii,  const DblNumVec&, DblNumVec&);
  int UpwEqu2DwnChk_dgemv (int level, Index3 ii,  const DblNumVec& effDen, DblNumVec& effVal);

  ///plain->regular->effective
  int plnDen2EffDen(int level, const DblNumVec&, DblNumVec&); 
  int samDen2RegDen(const DblNumVec&, DblNumVec&);
  ///effective->regular->plain
  int effVal2PlnVal(int level, const DblNumVec&, DblNumVec&); 
  int regVal2SamVal(const DblNumVec&, DblNumVec&);

  /*! Return sample positions of UE=0,UC=1,DE=2, or DC=3
	* where UE = upward equivalent,
	*       UC = upward check,
	*       DE = downward equivalent,
	*       DC = downward check
	*/
  const DblNumMat& samPos(int tp) { return _samPos[tp]; }
  /*! return regular positions */
  const DblNumMat& regPos()       { return _regPos; }
  /*! return local position */
  int localPos(int, Point3, double, DblNumMat&);
  /*! calculate various grids */
  int samPosCal(int n, double R, DblNumMat& ret);
  /*! calculate the fft regular position grid */
  int regPosCal(int n, double R, DblNumMat& ret); 
  int cptwvv(int, double, fftw_complex*, int, fftw_complex*, int, fftw_complex*, int);
  
public:
  static double _wsbuf[];
  static vector<MatMgnt3d> _mmvec;
public:
  static MatMgnt3d* getmmptr(Kernel3d, int);  //static void clearmmptrs();
};




/*	 int plnnum(ItlGrdType tp) { return _plnpos[tp].n(); }
	 int effnum(ItlGrdType tp) { return (2*_np+2)*(2*_np)*(2*_np); }
	 int regnum(ItlGrdType tp) { return _regpos[tp].n(); }  */
//int SE2TC_dgemv(const DblNumMat& srcpos, const DblNumMat& srcnor, const DblNumMat& trgpos, const DblNumVec& srcden, DblNumVec& trgval);
//int SE2UC_dgemv(const DblNumMat& srcpos, const DblNumMat& srcnor, Point3 trgctr, double trgrad, const DblNumVec& srcden, DblNumVec& trgval);
//int SE2DC_dgemv(const DblNumMat& srcpos, const DblNumMat& srcnor, Point3 trgctr, double trgrad, const DblNumVec& srcden, DblNumVec& trgval);
//int DE2TC_dgemv(Point3 srcctr, double srcrad, const DblNumMat& trgpos, const DblNumVec& srcden, DblNumVec& trgval);
//int UE2TC_dgemv(Point3 srcctr, double srcrad, const DblNumMat& trgpos, const DblNumVec& srcden, DblNumVec& trgval);




#endif


