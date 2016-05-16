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
#ifndef _DENSE3D_HPP_
#define _DENSE3D_HPP_

#include "knlmat3d.hpp"

using std::vector;

//! Dense3d implements the virtual setup and eval functions of KnlMat3d */
class Dense3d: public KnlMat3d
{
public:
  /*! Dense3d constructor takes string to declare type */
  Dense3d(const string& p);
  /*! Dense3d destructor */
  ~Dense3d();

  /* Setup the environment for dense3d.  Very little is done for this */
  int setup(map<string,string>&);
  /* Evaluate target data.  This is done directly, using a dense solver */
  int evaluate(const DblNumVec& srcDen, DblNumVec& trgVal);
};


#endif
