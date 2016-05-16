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
#ifndef _KNLMAT3D_HPP_
#define _KNLMAT3D_HPP_

#include "common/vec3t.hpp"
#include "comobject.hpp"
#include "kernel3d.hpp"

//----------------------------------------------------------------------------------
class KnlMat3d: public ComObject
{
protected:
  //---PARAMS (REQ)
  /*!source positions */
  DblNumMat* _srcPos;
  /*!source normal */
  DblNumMat* _srcNor;
  /*!target position */
  DblNumMat* _trgPos;
  /*! Kernel.  See kernel3d.hpp and kernel3d.cpp */
  Kernel3d _knl;
public:
  KnlMat3d(const string& p):  ComObject(p), _srcPos(NULL), _srcNor(NULL), _trgPos(NULL) {;}
  virtual ~KnlMat3d() { }
  //MEMBER ACESS
  /*! return source positions */
  DblNumMat*& srcPos() { return _srcPos; }
  /*! return source normals */
  DblNumMat*& srcNor() { return _srcNor; }
  /*! return target positions */
  DblNumMat*& trgPos() { return _trgPos; }
  /*! return kernel */
  Kernel3d& knl()    { return _knl; }
  //SETUP and USE
  /*! Virtual setup function */
  virtual int setup(map<string, string>&)=0;
  /*! Virtual evaluation function */
  virtual int evaluate(const DblNumVec& srcDen, DblNumVec& trgVal) = 0;
  //OTHER ACCESS
  /*! return dimension of the kernel matrix/fmm, set to 3 */  
  int dim() { return 3; }
  /*! return source degrees of freedom */
  int srcDOF() { return _knl.srcDOF(); }
  /*! return target degrees of freedom */
  int trgDOF() { return _knl.trgDOF(); }
};

#endif
