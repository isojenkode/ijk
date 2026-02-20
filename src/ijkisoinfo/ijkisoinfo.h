/*!
 *  @file ijkisoinfo.h
 *  @brief Get info about a "generic" isosurface in a scalar grid.
 *  - Version 0.6.0
 */


/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2025 Rephael Wenger

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

#ifndef _IJKISOINFO_
#define _IJKISOINFO_

#include "ijk.tpp"
#include "ijkgrid.tpp"
#include "ijkgrid_nrrd.tpp"

namespace IJKISOINFO {

  // *****************************************************************
  // TYPE DEFINITIONS
  // *****************************************************************

  const int GRID_CUBEV_BITSET_SIZE(8);
  
  typedef float SCALAR_TYPE;
  typedef int AXIS_SIZE_TYPE;
  typedef int VERTEX_INDEX;
  typedef int GRID_COORD_TYPE;
  typedef IJK::GRID_SIZE_TYPE NUM_TYPE;
  typedef IJK::GRID_PLUS
  <GRID_CUBEV_BITSET_SIZE,NUM_TYPE, AXIS_SIZE_TYPE, VERTEX_INDEX, NUM_TYPE> 
    INFO_GRID;
  typedef IJK::SCALAR_GRID<INFO_GRID, SCALAR_TYPE> SGRID;
  typedef IJK::GRID_NRRD_IN<int,AXIS_SIZE_TYPE> SGRID_NRRD_IN;
  
  /// @brief Command line options.
  typedef enum {
    SUBSAMPLE_OPT
  } COMMAND_LINE_OPTION_TYPE;
  

  // *****************************************************************
  // Class INPUT_ARG
  // *****************************************************************

  /// @brief Input arguments.
  class INPUT_ARG {

  protected:
    /// @brief Initialize data structure.
    void Init();
    
  public:    
    std::string nrrd_filename;
    std::vector<std::string> isovalue_string;
    std::vector<SCALAR_TYPE> isovalue;

    /// @brief If true, subsample the data set.
    bool flag_subsample;

    /// @brief Subsample resolution.
    int subsample_resolution;

  public:
    /// @brief Constructor.
    INPUT_ARG() { Init(); }
  };


  // *****************************************************************
  // Class NRRD_INFO
  // *****************************************************************

  class NRRD_INFO {
  public:
    /// @brief Indicates scalar type in nrrd input file.
    int scalar_type;
  };
    
};


#endif
