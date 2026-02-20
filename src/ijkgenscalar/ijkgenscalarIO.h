/*!
 *  @file ijkgenscalarIO.h
 *  @brief Read/write/prompt for ijkgenscalar
 *  - Version 0.6.0
 */

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2011-2025 Rephael Wenger

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


#ifndef _IJKGENSCALARIO_
#define _IJKGENSCALARIO_

#include "ijkgenscalar.h"

#include<string>
#include<vector>

namespace IJKGENSCALAR {


  // **************************************************
  // ENUMERATED TYPES
  // **************************************************

  typedef enum { NRRD_RAW, NRRD_ASCII, NRRD_GZIP }
    NRRD_ENCODING_TYPE;

  
  // **************************************************
  // WRITE NRRD FILE
  // **************************************************

  void write_scalar_grid
    (const char * output_filename, const SCALAR_GRID & grid, 
     const FIELD_PARAM & field_param,
     const std::vector<FIELD_INFO> & field_info,
     const bool flag_gzip);

  void write_gradient_grid
    (const std::string & output_filename, const GRADIENT_GRID & grid,
     const FIELD_PARAM & field_param, 
     const std::vector<FIELD_INFO> & field_info,
     const bool flag_gzip);

  void write_scalar_grid
    (const char * output_filename, const SCALAR_GRID & grid, 
     const FIELD_PARAM & field_param,
     const std::vector<FIELD_INFO> & field_info,
     const NRRD_ENCODING_TYPE nrrd_encoding_type);

};

#endif

