/*!
 *  @file ijkscalar.h
 *  @brief Combine or modify scalar and gradient grid data
 *  - Version 0.6.0
 */


/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2012-2025 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 3.0 of the License, or any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _IJKSCALAR_
#define _IJKSCALAR_

#include <cmath>
#include <iostream>
#include <string>

#include "ijkscalar_grid.tpp"
#include "ijkvector_grid.tpp"

/// Types and data structures for ijkscalar
namespace IJKSCALAR {

  // **************************************************
  // TYPES
  // **************************************************

  const int GRID_CUBEV_BITSET_SIZE(8);
  
  typedef float SCALAR_TYPE;
  typedef int AXIS_SIZE_TYPE;
  typedef int LENGTH_TYPE;
  typedef int VERTEX_INDEX;
  typedef int GRID_COORD_TYPE;
  typedef float COORD_TYPE;
  typedef float GRADIENT_COORD_TYPE;
  typedef GRADIENT_COORD_TYPE ANGLE_TYPE;
  typedef int NUM_TYPE;
  typedef IJK::GRID_PLUS
  <GRID_CUBEV_BITSET_SIZE,NUM_TYPE,AXIS_SIZE_TYPE,VERTEX_INDEX,NUM_TYPE> 
  GRID_PLUS;
  typedef IJK::GRID_SPACING<COORD_TYPE, GRID_PLUS> GRID_BASE;
  typedef IJK::SCALAR_GRID_BASE<GRID_BASE, SCALAR_TYPE> SGRID_BASE;
  typedef IJK::SCALAR_GRID<GRID_BASE, SCALAR_TYPE> SGRID;
  typedef IJK::VECTOR_GRID<GRID_BASE, LENGTH_TYPE, GRADIENT_COORD_TYPE> 
  GRADIENT_GRID;

};

#endif
