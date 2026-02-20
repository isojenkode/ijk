/*!
 *  \file ijkdual_position.tpp
 *  @brief Position dual isosurface vertices.
 *  - Version 0.6.0
 */


/*
  IJK: Isosurface Jeneration Kode
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

#ifndef IJKDUAL_POSITION_TPP_
#define IJKDUAL_POSITION_TPP_

#include <bitset>
#include <unordered_map>
#include <vector>

#ifdef INCLUDE_RANDOM
#include "ijkrandom.tpp"
#endif

#include "ijkdual_datastruct.h"
#include "ijkdual_types.h"

#include "ijk.tpp"
#include "ijkcoord.tpp"
#include "ijkcube.tpp"
#include "ijkhash.tpp"
#include "ijkinterpolate.tpp"
#include "ijkisocoord.tpp"
#include "ijkisopoly.tpp"
#include "ijkscalar_grid.tpp"
#include "ijksort.tpp"


// *** DEBUG ***
#include "ijkprint.tpp"


namespace IJKDUAL {

  // *****************************************************************
  //! @name Position single isosurface vertex in each grid cube.
  // *****************************************************************

  //@{

  /*!
   *  @brief Position dual isosurface vertices in cube centers.
   *  - One isosurface vertex per active cube.
   *  @param active_cube_list[kw] Cube containing isosurface vertex kw.
   *    - One isosurface vertex per cube.
   *  @param[out] coord[] Isosurface vertex coordinates.
   *    - coord[i*dimension+d] = d'th coordinate
   *      of isosurface vertex in cube active_cube_list[i].
   *    @pre Array coord[] is preallocated to size at least
   *    active_cube_list.size() * grid.Dimension().
   */
  template <typename GRID_TYPE, typename ITYPE,
            typename CTYPE>
  void position_all_isov_cube_center_single
  (const GRID_TYPE & grid,
   const std::vector<ITYPE> & active_cube_list, CTYPE * coord)
  {
    typedef typename std::vector<ITYPE>::size_type SIZE_TYPE;

    const int dimension = grid.Dimension();

    for (SIZE_TYPE i = 0; i < active_cube_list.size(); i++) {
      const ITYPE icube = active_cube_list[i];

      grid.ComputeCubeCenterCoord(icube, coord+i*dimension);
    }
  }


  /// @brief Position dual isosurface vertices in cube centers.
  /// - C++ STL vector format for array coord[].
  template <typename GRID_TYPE, typename ITYPE,
            typename CTYPE>
  void position_all_isov_cube_center_single
  (const GRID_TYPE & grid,
   const std::vector<ITYPE> & active_cube_list,
   std::vector<CTYPE> & coord)
  {
    const int dimension = grid.Dimension();

    coord.resize(active_cube_list.size()*dimension);
    position_all_isov_cube_center_single
      (grid, active_cube_list, IJK::vector2pointerNC(coord));
  }


  /*!
   *  @brief Position dual isosurface vertex in centroid
   *  of isosurface-cube edge intersections.
   *  - One isosurface vertex per active cube.
   *  @param icube Index of cube containing isosurface vertex.
   *  @param[out] coord[] Isosurface vertex coordinates.
   *    - coord[i*dimension+d] = d'th coordinate
   *      of isosurface vertex in cube icube
   *    @pre Array coord[] is preallocated to size at least
   *      grid.Dimension().
   */
  template <typename GRID_TYPE, typename STYPE, typename ITYPE,
            typename CTYPEV, typename CTYPET>
  void position_isov_in_cube_centroid_single
  (const GRID_TYPE & scalar_grid, const STYPE isovalue,
   const ITYPE icube,
   CTYPEV * vcoord, CTYPET * temp_coord)
  {
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;

    const int dimension = scalar_grid.Dimension();

    NTYPE num_intersected_edges = 0;
    IJK::set_coord(dimension, 0.0, vcoord);

    for (int edge_dir = 0; edge_dir < dimension; edge_dir++)
      for (NTYPE k = 0; k < scalar_grid.NumCubeFacetVertices(); k++) {
        const VERTEX_INDEX_TYPE iend0 =
          scalar_grid.FacetVertex(icube, edge_dir, k);
        const VERTEX_INDEX_TYPE iend1 =
          scalar_grid.NextVertex(iend0, edge_dir);

        const STYPE s0 = scalar_grid.Scalar(iend0);
        bool is_end0_positive = true;
        if (s0 < isovalue)
          { is_end0_positive = false; };

        const STYPE s1 = scalar_grid.Scalar(iend1);
        bool is_end1_positive = true;
        if (s1 < isovalue)
          { is_end1_positive = false; };

        if (is_end0_positive != is_end1_positive) {
          IJK::compute_isosurface_grid_edge_intersection_linear
            (scalar_grid, isovalue, iend0, edge_dir, s0, s1, temp_coord);
            
          IJK::add_coord
            (dimension, vcoord, temp_coord, vcoord);

          num_intersected_edges++;
        }

      }

    if (num_intersected_edges > 0) {
      IJK::multiply_coord
        (dimension, 1.0/num_intersected_edges, vcoord, vcoord);
    }
    else {
      scalar_grid.ComputeCubeCenterCoord(icube, vcoord);
    }
  }


  /*!
   *  @brief Position all dual isosurface vertices in centroid
   *  of isosurface-edge intersections. (Simple version.)
   *  - One isosurface vertex per active cube.
   *  - Simple version that does not precompute intersections
   *    of isosurface and grid edges.
   *  @param active_cube_list[kw] Cube containing isosurface vertex kw.
   *    - One isosurface vertex per cube.
   */
  template <typename GRID_TYPE, typename STYPE, typename ITYPE,
            typename CTYPE>
  void position_all_isov_centroid_single_simple
  (const GRID_TYPE & scalar_grid,
   const STYPE isovalue,
   const std::vector<ITYPE> & active_cube_list,
   CTYPE * coord)
  {
    typedef typename std::vector<ITYPE>::size_type SIZE_TYPE;

    const int dimension = scalar_grid.Dimension();
    IJK::ARRAY<CTYPE> temp_coord(dimension);

    for (SIZE_TYPE i = 0; i < active_cube_list.size(); i++) {
      const ITYPE icube = active_cube_list[i];

      position_isov_in_cube_centroid_single
        (scalar_grid, isovalue, icube, coord+i*dimension,
         temp_coord.Ptr());
    }
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices in centroid
   *    of isosurface-edge intersections. (Simple version. C++ vector.)
   *  - One isosurface vertex per active cube.
   *  - Simple version that does not precompute intersections
   *    of isosurface and grid edges.
   *  - C++ STL vector format for array coord[].
   */
  template <typename GRID_TYPE, typename STYPE, typename ITYPE,
            typename CTYPE>
  void position_all_isov_centroid_single_simple
  (const GRID_TYPE & scalar_grid,
   const STYPE isovalue,
   const std::vector<ITYPE> & active_cube_list,
   std::vector<CTYPE> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(active_cube_list.size()*dimension);
    position_all_isov_centroid_single_simple
      (scalar_grid, isovalue, active_cube_list,
       IJK::vector2pointerNC(coord));
  }


  /*!
   *  @brief Position dual isosurface vertex in centroid
   *  of isosurface-cube edge intersections with offset.
   *  - One isosurface vertex per active cube.
   *  - Offset isosurface-(grid edge) intersections so that all
   *    centroid is in cube interior.
   *  @param offset Offset value.
   *    @pre Value of offset must be in range [0,0.5].
   *  @param[out] coord[] Isosurface vertex coordinates.
   *    - coord[i*dimension+d] = d'th coordinate
   *      of isosurface vertex in cube icube
   *    @pre Array coord[] is preallocated to size at least
   *      grid.Dimension().
   */
  template <typename GRID_TYPE, typename STYPE, typename ITYPE,
            typename OFFSET_TYPE, typename CTYPEV, typename CTYPET>
  void position_isov_in_cube_centroid_offset_single
  (const GRID_TYPE & scalar_grid, const STYPE isovalue,
   const ITYPE icube, const OFFSET_TYPE offset,
   CTYPEV * vcoord, CTYPET * temp_coord)
  {
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;

    const int dimension = scalar_grid.Dimension();

    NTYPE num_intersected_edges = 0;
    IJK::set_coord(dimension, 0.0, vcoord);

    for (int edge_dir = 0; edge_dir < dimension; edge_dir++)
      for (NTYPE k = 0; k < scalar_grid.NumCubeFacetVertices(); k++) {
        const VERTEX_INDEX_TYPE iend0 =
          scalar_grid.FacetVertex(icube, edge_dir, k);
        const VERTEX_INDEX_TYPE iend1 =
          scalar_grid.NextVertex(iend0, edge_dir);

        const STYPE s0 = scalar_grid.Scalar(iend0);
        bool is_end0_positive = true;
        if (s0 < isovalue)
          { is_end0_positive = false; };

        const STYPE s1 = scalar_grid.Scalar(iend1);
        bool is_end1_positive = true;
        if (s1 < isovalue)
          { is_end1_positive = false; };

        if (is_end0_positive != is_end1_positive) {

          IJK::compute_isosurface_grid_edge_intersection_linear_offset
            (scalar_grid, isovalue, iend0, edge_dir, s0, s1,
             offset, temp_coord);

          IJK::add_coord(dimension, vcoord, temp_coord, vcoord);

          num_intersected_edges++;
        }

      }

    if (num_intersected_edges > 0) {
      IJK::multiply_coord
        (dimension, 1.0/num_intersected_edges, vcoord, vcoord);
    }
    else {
      scalar_grid.ComputeCubeCenterCoord(icube, vcoord);
    }
  }


  /*!
   *  @brief Position all dual isosurface vertices in centroid
   *  of isosurface-edge intersections. Simple version with offset.
   *  - One isosurface vertex per active cube.
   *  - Simple version that does not precompute intersections
   *    of isosurface and grid edges.
   *  - Offset isosurface-(grid edge) intersections so that all
   *    centroid is in cube interior.
   *  @param active_cube_list[kw] Cube containing isosurface vertex kw.
   *    - One isosurface vertex per cube.
   *  @param offset Offset value.
   *    @pre Value of offset must be in range [0,0.5].
   */
  template <typename GRID_TYPE, typename STYPE, typename ITYPE,
            typename OFFSET_TYPE, typename CTYPE>
  void position_all_isov_centroid_offset_single_simple
  (const GRID_TYPE & scalar_grid,
   const STYPE isovalue,
   const std::vector<ITYPE> & active_cube_list,
   const OFFSET_TYPE offset,
   CTYPE * coord)
  {
    typedef typename std::vector<ITYPE>::size_type SIZE_TYPE;

    const int dimension = scalar_grid.Dimension();
    IJK::ARRAY<CTYPE> temp_coord(dimension);

    for (SIZE_TYPE i = 0; i < active_cube_list.size(); i++) {
      const ITYPE icube = active_cube_list[i];

      position_isov_in_cube_centroid_offset_single
        (scalar_grid, isovalue, icube, offset,
         coord+i*dimension, temp_coord.Ptr());
    }
  }
  
  /*!
   *  @overload
   *  @brief Position dual isosurface vertices in centroid
   *    of isosurface-edge intersections. Simple version with offset.
   *    (C++ vector.)
   *  - One isosurface vertex per active cube.
   *  - Simple version that does not precompute intersections
   *    of isosurface and grid edges.
   *  - Offset isosurface-(grid edge) intersections so that all
   *    centroid is in cube interior.
   *  - C++ STL vector format for array coord[].
   */
  template <typename GRID_TYPE, typename STYPE, typename ITYPE,
            typename OFFSET_TYPE, typename CTYPE>
  void position_all_isov_centroid_offset_single_simple
  (const GRID_TYPE & scalar_grid,
   const STYPE isovalue,
   const std::vector<ITYPE> & active_cube_list,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(active_cube_list.size()*dimension);
    position_all_isov_centroid_offset_single_simple
      (scalar_grid, isovalue, active_cube_list, offset,
       IJK::vector2pointerNC(coord));
  }

  
  /*!
   *  @brief Position dual isosurface vertices in centroid
   *  of isosurface-edge intersections.
   *  @param active_cube_list[kw] Cube containing isosurface vertex kw.
   *    - One isosurface vertex per cube.
   *  @param[out] coord[] Isosurface vertex coordinates.
   *    - coord[i*dimension+d] = d'th coordinate
   *      of isosurface vertex in cube active_cube_list[i].
   *    @pre Array coord[] is preallocated to size at least
   *    active_cube_list.size() * grid.Dimension().
   */
  template <typename GRID_TYPE, typename STYPE, typename ITYPE,
            typename ISOV_INDEX_TYPE, typename ISOPOLY_INFO_TYPE,
            typename CTYPE>
  void position_all_isov_centroid_single
  (const GRID_TYPE & scalar_grid,
   const STYPE isovalue,
   const std::vector<ITYPE> & active_cube_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   CTYPE * isov_coord)
  {
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;
    typedef typename std::vector<ITYPE>::size_type SIZE_TYPE;

    const int dimension = scalar_grid.Dimension();
    const int num_vertices_per_isopoly =
      scalar_grid.NumFacetVertices();
    const NTYPE num_isov = active_cube_list.size();
    const NTYPE num_isopoly = isopoly_info.size();
    IJK::ARRAY<CTYPE> temp_coord(dimension);
    IJK::PROCEDURE_ERROR
      error("position_all_isov_centroid_single");

    if (!IJK::check_array_allocated
        (isov_coord, "isov_coord", error)) { throw error; }

    if (!check_isopoly_vert_size
        (isopoly_vert, num_vertices_per_isopoly, isopoly_info.size(),
         error)) { throw error; }

    // num_incident_isopoly[i] =
    //   Number of isosurface polytopes incident on isov_list[i].
    std::vector<ISO_VERTEX_DEGREE>
      num_incident_isopoly(num_isov, 0);

    // Set all isov_coord entries to 0.0.
    for (NTYPE i = 0; i < num_isov*dimension; i++)
      { isov_coord[i] = 0.0; }

    SIZE_TYPE k = 0;
    for (NTYPE ipoly = 0; ipoly < num_isopoly; ipoly++) {
      const VTYPE iend0 = isopoly_info[ipoly].GridEdgeEndpoint0();
      const int edge_direction = isopoly_info[ipoly].GridEdgeDirection();
      
      IJK::compute_isosurface_grid_edge_intersection_linear
        (scalar_grid, isovalue, iend0, edge_direction,
         temp_coord.Ptr());

      for (int j = 0; j < num_vertices_per_isopoly; j++) {
        const ISOV_INDEX_TYPE kw = isopoly_vert[k];
        CTYPE * isov_coord_ptr = isov_coord + dimension*kw;

        IJK::add_coord
          (dimension, isov_coord_ptr, temp_coord.Ptr(), isov_coord_ptr);
        num_incident_isopoly[kw]++;
        k++;
      }
    }

    if (k != isopoly_vert.size()) {
      error.AddMessage
        ("Programming error. Processed incorrect number of isopoly vertices.");
      error.AddMessage
        ("  Isosurface polytopes are incident on ",
         isopoly_vert.size(), " vertices.");
      error.AddMessage("  Processed ", k, " vertices.");
      throw error;
    }

    for (NTYPE kw = 0; kw < num_isov; kw++) {
      const ITYPE icube = active_cube_list[kw];
      CTYPE * isov_coord_ptr = isov_coord + dimension*kw;

      if (scalar_grid.IsCubeOnGridBoundary(icube)) {
        // Isosurface vertex is in cube on grid boundary.
        // Some dual facets may be missing from isopoly_info.
        // Use slower method to compute isosurface vertex position.
        position_isov_in_cube_centroid_single
          (scalar_grid, isovalue, icube, isov_coord_ptr,
           temp_coord.Ptr());
      }
      else {
        IJK::multiply_coord
          (dimension, 1.0/num_incident_isopoly[kw],
           isov_coord_ptr, isov_coord_ptr);
      }
    }
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices in centroid
   *  of isosurface-edge intersections. (C++ vector.)
   *  - C++ STL vector format for array coord[].
   *  @param[out] coord[] Isosurface vertex coordinates.
   *    - coord[i*dimension+d] = d'th coordinate
   *      of isosurface vertex in cube active_cube_list[i].
   */
  template <typename GRID_TYPE, typename STYPE, typename ITYPE,
            typename ISOV_INDEX_TYPE, typename ISOPOLY_INFO_TYPE,
            typename CTYPE>
  void position_all_isov_centroid_single
  (const GRID_TYPE & scalar_grid,
   const STYPE isovalue,
   const std::vector<ITYPE> & active_cube_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   std::vector<CTYPE> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(active_cube_list.size()*dimension);
    position_all_isov_centroid_single
      (scalar_grid, isovalue, active_cube_list,
       isopoly_vert, isopoly_info, IJK::vector2pointerNC(coord));
  }


  /*!
   *  @brief Position dual isosurface vertices in centroid
   *  of isosurface-edge intersections with offset.
   *  - Offset isosurface-(grid edge) intersections so that all
   *    centroid is in cube interior.
   *  @param active_cube_list[kw] Cube containing isosurface vertex kw.
   *    - One isosurface vertex per cube.
   *  @param offset Offset value.
   *    @pre Value of offset must be in range [0,0.5].
   *  @param[out] coord[] Isosurface vertex coordinates.
   *    - coord[i*dimension+d] = d'th coordinate
   *      of isosurface vertex in cube active_cube_list[i].
   *    @pre Array coord[] is preallocated to size at least
   *    active_cube_list.size() * grid.Dimension().
   */
  template <typename GRID_TYPE, typename STYPE, typename ITYPE,
            typename ISOV_INDEX_TYPE, typename ISOPOLY_INFO_TYPE,
            typename OFFSET_TYPE, typename CTYPE>
  void position_all_isov_centroid_offset_single
  (const GRID_TYPE & scalar_grid,
   const STYPE isovalue,
   const std::vector<ITYPE> & active_cube_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   const OFFSET_TYPE offset,
   CTYPE * isov_coord)
  {
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;
    typedef typename std::vector<ITYPE>::size_type SIZE_TYPE;

    const int dimension = scalar_grid.Dimension();
    const int num_vertices_per_isopoly =
      scalar_grid.NumFacetVertices();
    const NTYPE num_isov = active_cube_list.size();
    const NTYPE num_isopoly = isopoly_info.size();
    IJK::ARRAY<CTYPE> temp_coord(dimension);
    IJK::PROCEDURE_ERROR
      error("position_all_isov_centroid_offset_single");

    if (!IJK::check_array_allocated
        (isov_coord, "isov_coord", error)) { throw error; }

    if (!check_isopoly_vert_size
        (isopoly_vert, num_vertices_per_isopoly, isopoly_info.size(),
         error)) { throw error; }

    // num_incident_isopoly[i] =
    //   Number of isosurface polytopes incident on isov_list[i].
    std::vector<ISO_VERTEX_DEGREE>
      num_incident_isopoly(num_isov, 0);

    // Set all isov_coord entries to 0.0.
    for (NTYPE i = 0; i < num_isov*dimension; i++)
      { isov_coord[i] = 0.0; }

    SIZE_TYPE k = 0;
    for (NTYPE ipoly = 0; ipoly < num_isopoly; ipoly++) {
      const VTYPE iend0 = isopoly_info[ipoly].GridEdgeEndpoint0();
      const int edge_direction = isopoly_info[ipoly].GridEdgeDirection();

      IJK::compute_isosurface_grid_edge_intersection_linear_offset
        (scalar_grid, isovalue, iend0, edge_direction,
         offset, temp_coord.Ptr());

      for (int j = 0; j < num_vertices_per_isopoly; j++) {
        const ISOV_INDEX_TYPE kw = isopoly_vert[k];
        CTYPE * isov_coord_ptr = isov_coord + dimension*kw;

        IJK::add_coord
          (dimension, isov_coord_ptr, temp_coord.Ptr(), isov_coord_ptr);
        num_incident_isopoly[kw]++;
        k++;
      }
    }

    if (k != isopoly_vert.size()) {
      error.AddMessage
        ("Programming error. Processed incorrect number of isopoly vertices.");
      error.AddMessage
        ("  Isosurface polytopes are incident on ",
         isopoly_vert.size(), " vertices.");
      error.AddMessage("  Processed ", k, " vertices.");
      throw error;
    }

    for (NTYPE kw = 0; kw < num_isov; kw++) {
      const ITYPE icube = active_cube_list[kw];
      CTYPE * isov_coord_ptr = isov_coord + dimension*kw;

      if (scalar_grid.IsCubeOnGridBoundary(icube)) {
        // Isosurface vertex is in cube on grid boundary.
        // Some dual facets may be missing from isopoly_info.
        // Use slower method to compute isosurface vertex position.
        position_isov_in_cube_centroid_offset_single
          (scalar_grid, isovalue, icube, offset,
           isov_coord_ptr, temp_coord.Ptr());
      }
      else {
        IJK::multiply_coord
          (dimension, 1.0/num_incident_isopoly[kw],
           isov_coord_ptr, isov_coord_ptr);
      }
    }
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices in centroid
   *  of isosurface-edge intersections. (C++ vector.)
   *  - C++ STL vector format for array coord[].
   *  @param[out] coord[] Isosurface vertex coordinates.
   *    - coord[i*dimension+d] = d'th coordinate
   *      of isosurface vertex in cube active_cube_list[i].
   */
  template <typename GRID_TYPE, typename STYPE, typename ITYPE,
            typename ISOV_INDEX_TYPE, typename ISOPOLY_INFO_TYPE,
            typename OFFSET_TYPE, typename CTYPE>
  void position_all_isov_centroid_offset_single
  (const GRID_TYPE & scalar_grid,
   const STYPE isovalue,
   const std::vector<ITYPE> & active_cube_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(active_cube_list.size()*dimension);
    position_all_isov_centroid_offset_single
      (scalar_grid, isovalue, active_cube_list,
       isopoly_vert, isopoly_info, offset,
       IJK::vector2pointerNC(coord));
  }


  //@}


  // *****************************************************************
  //! @name Position multiple isosurface vertices in each grid cube.
  // *****************************************************************

  //@{

  /*!
   *  @brief Position dual isosurface vertex in centroid
   *    of isosurface patch-cube edge intersections.
   *  - Cube may contain multiple isosurface vertices.
   *  - Only use bipolar cube edges associated with the
   *    isosurface vertex, i.e., that intersect the
   *    isosurface patch containing the isosurface vertex.
   *  @pre scalar_grid.Dimension() == cube.Dimension().
   *  @param[out] vcoord[d] d'th coordinate of isosurface vertex.
   *    @pre Array vcoord[] is preallocated to size at least
   *      scalar_grid.Dimension().
   *  @param temp_coord[] Temporary coordinate array.
   *    @pre Pre-allocated to size at least dimension.
   */
  template <typename GRID_TYPE,
            typename ISODUAL_TABLE_TYPE,
            typename STYPE,
            typename DUAL_ISOV_TYPE,
            typename CUBE_TYPE,
            typename CTYPE, typename CTYPE2>
  void position_isov_in_cube_centroid_multi
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue, const DUAL_ISOV_TYPE & isov_info,
   const CUBE_TYPE & cube,
   CTYPE * vcoord, CTYPE2 * temp_coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;

    const DTYPE dimension = scalar_grid.Dimension();
    const VTYPE icube = isov_info.cube_index;
    const NTYPE ipatch = isov_info.patch_index;
    const TABLE_INDEX table_index = isov_info.table_index;
    const NTYPE num_incident_isopoly =
      isodual_table.NumIncidentIsoPoly(table_index, ipatch);

    IJK::set_coord(dimension, 0.0, vcoord);

    for (NTYPE ie = 0; ie < cube.NumEdges(); ie++) {
      if (isodual_table.IsBipolar(table_index, ie)) {
        if (isodual_table.IncidentIsoVertex(table_index, ie) == ipatch) {
          const NTYPE k0 = cube.EdgeEndpoint(ie, 0);
          const VTYPE iend0 = scalar_grid.CubeVertex(icube, k0);
          const int edge_direction = cube.EdgeDir(ie);
          
          compute_isosurface_grid_edge_intersection_linear
            (scalar_grid, isovalue, iend0, edge_direction,
             temp_coord);

          IJK::add_coord(dimension, vcoord, temp_coord, vcoord);
        }
      }
    }

    if (num_incident_isopoly > 0) {
      IJK::multiply_coord
        (dimension, 1.0/num_incident_isopoly, vcoord, vcoord);
    }
    else {
      scalar_grid.ComputeCubeCenterCoord(icube, vcoord);
    }
  }


  /*!
   *  @brief Position dual isosurface vertices in centroid
   *    of isosurface-edge intersections. (Check intersection.)
   *  - Check iso grid edge intersection.
   *  - Cube may contain multiple isosurface vertices.
   *  - If isosurface does not intersect grid edge,
   *    use grid edge midpoint in computing centroid.
   *  - Used in computing positions in cubes on grid ridges,
   *    where isosurface lookup table index may be changed
   *    to position two isosurfac vertices in the cube.
   *    - Table index is changed to avoid non-manifold conditions.
   *    - Some (boundary) edges may be identified as bipolar,
   *      even though the scalar values at the endpoints are both
   *      above or both below the isovalue.
   *  @pre scalar_grid.Dimension() == cube.Dimension().
   *  @param[out] vcoord[d] d'th coordinate of isosurface vertex.
   *    @pre Array vcoord[] is preallocated to size at least
   *      scalar_grid.Dimension().
   *  @param temp_coord[] Temporary coordinate array.
   *    @pre Pre-allocated to size at least dimension.
   */
  template <typename GRID_TYPE,
            typename ISODUAL_TABLE_TYPE,
            typename STYPE,
            typename DUAL_ISOV_TYPE,
            typename CUBE_TYPE,
            typename CTYPE, typename CTYPE2>
  void position_isov_in_cube_centroid_multi_checkI
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue, const DUAL_ISOV_TYPE & isov_info,
   const CUBE_TYPE & cube,
   CTYPE * vcoord, CTYPE2 * temp_coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;

    const DTYPE dimension = scalar_grid.Dimension();
    const VTYPE icube = isov_info.cube_index;
    const NTYPE ipatch = isov_info.patch_index;
    const TABLE_INDEX table_index = isov_info.table_index;
    const NTYPE num_incident_isopoly =
      isodual_table.NumIncidentIsoPoly(table_index, ipatch);

    IJK::set_coord(dimension, 0.0, vcoord);

    for (NTYPE ie = 0; ie < cube.NumEdges(); ie++) {
      if (isodual_table.IsBipolar(table_index, ie)) {
        if (isodual_table.IncidentIsoVertex(table_index, ie) == ipatch) {
          const NTYPE k0 = cube.EdgeEndpoint(ie, 0);
          const VTYPE iend0 = scalar_grid.CubeVertex(icube, k0);
          const int edge_direction = cube.EdgeDir(ie);
          const VTYPE iend1 =
            scalar_grid.NextVertex(iend0, edge_direction);

          const STYPE s0 = scalar_grid.Scalar(iend0);
          bool is_end0_positive = true;
          if (s0 < isovalue)
            { is_end0_positive = false; };

          const STYPE s1 = scalar_grid.Scalar(iend1);
          bool is_end1_positive = true;
          if (s1 < isovalue)
            { is_end1_positive = false; };

          if (is_end0_positive == is_end1_positive) {
            scalar_grid.ComputeEdgeMidpoint
              (iend0, edge_direction, temp_coord);
          }
          else {
            compute_isosurface_grid_edge_intersection_linear
              (scalar_grid, isovalue, iend0, edge_direction,
               s0, s1, temp_coord);
          }

          IJK::add_coord(dimension, vcoord, temp_coord, vcoord);
        }
      }
    }

    if (num_incident_isopoly > 0) {
      IJK::multiply_coord
        (dimension, 1.0/num_incident_isopoly, vcoord, vcoord);
    }
    else {
      scalar_grid.ComputeCubeCenterCoord(icube, vcoord);
    }
  }


  /*!
   *  @brief Position all dual isosurface vertex in centroid
   *    of isosurface-edge intersections. (Simple version.)
   *  - A single cube may contain multiple isosurface vertices.
   */
  template <typename GRID_TYPE,
            typename ISODUAL_TABLE_TYPE,
            typename STYPE,
            typename DUAL_ISOV_TYPE, typename CTYPE>
  void position_all_isov_centroid_multi_simple
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   CTYPE * coord)
  {
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;

    const int dimension = scalar_grid.Dimension();
    IJK::ARRAY<CTYPE> temp_coord(dimension);
    IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);

    for (SIZE_TYPE kw = 0; kw < isov_list.size(); kw++) {
      position_isov_in_cube_centroid_multi_checkI
        (scalar_grid, isodual_table, isovalue, isov_list[kw], cube,
         coord+kw*dimension, temp_coord.Ptr());
    }
  }


  /*!
   *  @overload
   *  @brief Position all dual isosurface vertex in centroid
   *    of isosurface-edge intersections. (Simple version. C++ vector.)
   *  - Version using C++ STL vector for array coord[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename DUAL_ISOV_TYPE, typename CTYPE>
  void position_all_isov_centroid_multi_simple
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   std::vector<CTYPE> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(isov_list.size()*dimension);
    position_all_isov_centroid_multi_simple
      (scalar_grid, isodual_table,isovalue, isov_list,
       IJK::vector2pointerNC(coord));
  }


  /*!
   *  @brief Position dual isosurface vertex in centroid
   *    of isosurface patch-cube edge intersections with offsets.
   *  - Cube may contain multiple isosurface vertices.
   *  - Only use bipolar cube edges associated with the
   *    isosurface vertex, i.e., that intersect the
   *    isosurface patch containing the isosurface vertex.
   *  - Offset isosurface-(grid edge) intersections so that all
   *    centroid is in cube interior.
   *  - If isosurface does not intersect grid edge,
   *    use grid edge midpoint in computing centroid.
   *  - Used in computing positions in cubes on grid ridges,
   *    where isosurface lookup table index may be changed
   *    to position two isosurfac vertices in the cube.
   *    - Table index is changed to avoid non-manifold conditions.
   *    - Some (boundary) edges may be identified as bipolar,
   *      even though the scalar values at the endpoints are both
   *      above or both below the isovalue.   
   *  @pre scalar_grid.Dimension() == cube.Dimension().
   *  @param offset Offset value.
   *    @pre Value of offset must be in range [0,0.5].
   *  @param[out] vcoord[d] d'th coordinate of isosurface vertex.
   *    @pre Array vcoord[] is preallocated to size at least
   *      scalar_grid.Dimension().
   *  @param temp_coord[] Temporary coordinate array.
   *    @pre Pre-allocated to size at least dimension.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE, typename CUBE_TYPE,
            typename CTYPE, typename CTYPE2>
  void position_isov_in_cube_centroid_offset_multi
  (const GRID_TYPE & scalar_grid, const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue, const DUAL_ISOV_TYPE & isov_info,
   const OFFSET_TYPE offset, const CUBE_TYPE & cube,
   CTYPE * vcoord, CTYPE2 * temp_coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;

    const DTYPE dimension = scalar_grid.Dimension();
    const VTYPE icube = isov_info.cube_index;
    const NTYPE ipatch = isov_info.patch_index;
    const TABLE_INDEX table_index = isov_info.table_index;
    const NTYPE num_incident_isopoly =
      isodual_table.NumIncidentIsoPoly(table_index, ipatch);

    IJK::set_coord(dimension, 0.0, vcoord);

    for (NTYPE ie = 0; ie < cube.NumEdges(); ie++) {
      if (isodual_table.IsBipolar(table_index, ie)) {
        if (isodual_table.IncidentIsoVertex(table_index, ie) == ipatch) {
          const NTYPE k0 = cube.EdgeEndpoint(ie, 0);
          const VTYPE iend0 = scalar_grid.CubeVertex(icube, k0);
          const int edge_direction = cube.EdgeDir(ie);

          compute_isosurface_grid_edge_intersection_linear_offset
            (scalar_grid, isovalue, iend0, edge_direction, offset,
             temp_coord);

          IJK::add_coord(dimension, vcoord, temp_coord, vcoord);
        }
      }
    }

    if (num_incident_isopoly > 0) {
      IJK::multiply_coord
        (dimension, 1.0/num_incident_isopoly, vcoord, vcoord);
    }
    else {
      scalar_grid.ComputeCubeCenterCoord(icube, vcoord);
    }
  }


  /*!
   *  @brief Position dual isosurface vertex in centroid
   *    of isosurface patch-cube edge intersections with offsets.
   *    (Check intersection.)
   *  - Check iso grid edge intersection.
   *  - Cube may contain multiple isosurface vertices.
   *  - Only use bipolar cube edges associated with the
   *    isosurface vertex, i.e., that intersect the
   *    isosurface patch containing the isosurface vertex.
   *  - Offset isosurface-(grid edge) intersections so that all
   *    centroid is in cube interior.
   *  @pre scalar_grid.Dimension() == cube.Dimension().
   *  @param offset Offset value.
   *    @pre Value of offset must be in range [0,0.5].
   *  @param[out] vcoord[d] d'th coordinate of isosurface vertex.
   *    @pre Array vcoord[] is preallocated to size at least
   *      scalar_grid.Dimension().
   *  @param temp_coord[] Temporary coordinate array.
   *    @pre Pre-allocated to size at least dimension.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE, typename CUBE_TYPE,
            typename CTYPE, typename CTYPE2>
  void position_isov_in_cube_centroid_offset_multi_checkI
  (const GRID_TYPE & scalar_grid, const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue, const DUAL_ISOV_TYPE & isov_info,
   const OFFSET_TYPE offset, const CUBE_TYPE & cube,
   CTYPE * vcoord, CTYPE2 * temp_coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;

    const DTYPE dimension = scalar_grid.Dimension();
    const VTYPE icube = isov_info.cube_index;
    const NTYPE ipatch = isov_info.patch_index;
    const TABLE_INDEX table_index = isov_info.table_index;
    const NTYPE num_incident_isopoly =
      isodual_table.NumIncidentIsoPoly(table_index, ipatch);

    IJK::set_coord(dimension, 0.0, vcoord);

    for (NTYPE ie = 0; ie < cube.NumEdges(); ie++) {
      if (isodual_table.IsBipolar(table_index, ie)) {
        if (isodual_table.IncidentIsoVertex(table_index, ie) == ipatch) {
          const NTYPE k0 = cube.EdgeEndpoint(ie, 0);
          const VTYPE iend0 = scalar_grid.CubeVertex(icube, k0);
          const int edge_direction = cube.EdgeDir(ie);
          const VTYPE iend1 =
            scalar_grid.NextVertex(iend0, edge_direction);

          const STYPE s0 = scalar_grid.Scalar(iend0);
          bool is_end0_positive = true;
          if (s0 < isovalue)
            { is_end0_positive = false; };

          const STYPE s1 = scalar_grid.Scalar(iend1);
          bool is_end1_positive = true;
          if (s1 < isovalue)
            { is_end1_positive = false; };

          if (is_end0_positive == is_end1_positive) {
            scalar_grid.ComputeEdgeMidpoint
              (iend0, edge_direction, temp_coord);
          }
          else {
            compute_isosurface_grid_edge_intersection_linear_offset
              (scalar_grid, isovalue, iend0, edge_direction,
               s0, s1, offset, temp_coord);
          }

          IJK::add_coord(dimension, vcoord, temp_coord, vcoord);
        }
      }
    }

    if (num_incident_isopoly > 0) {
      IJK::multiply_coord
        (dimension, 1.0/num_incident_isopoly, vcoord, vcoord);
    }
    else {
      scalar_grid.ComputeCubeCenterCoord(icube, vcoord);
    }
  }


  /*!
   *  @brief Position dual isosurface vertices in centroid
   *    of isosurface-edge intersections with offsets. Simple version.
   */
  template <typename GRID_TYPE,typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE,
            typename CTYPE>
  void position_all_isov_centroid_offset_multi_simple
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const OFFSET_TYPE offset,
   CTYPE * coord)
  {
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;

    const int dimension = scalar_grid.Dimension();
    IJK::ARRAY<CTYPE> temp_coord(dimension);
    IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);

    for (SIZE_TYPE kw = 0; kw < isov_list.size(); kw++) {
      position_isov_in_cube_centroid_offset_multi_checkI
        (scalar_grid, isodual_table, isovalue, isov_list[kw],
         offset, cube, coord+kw*dimension, temp_coord.Ptr());
    }
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices in centroid
   *    of isosurface-edge intersections with offsets. Simple version.
   *    (C++ vector.)
   *  - Version using C++ STL vector for array coord[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE, typename CTYPE>
  void position_all_isov_centroid_offset_multi_simple
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(isov_list.size()*dimension);
    position_all_isov_centroid_offset_multi_simple
      (scalar_grid, isodual_table,isovalue, isov_list, offset,
       IJK::vector2pointerNC(coord));
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices in centroid
   *    of isosurface-edge intersections.
   *    Simple version with offset argument.
   *    - Version with argument ISO_GRID_EDGE_INTERSECTION_OFFSET_METHOD.
   */
  template <typename GRID_TYPE,typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE,
            typename CTYPE>
  void position_all_isov_centroid_offset_multi_simple
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const ISO_GRID_EDGE_INTERSECTION_OFFSET_METHOD offset_method,
   const OFFSET_TYPE offset,
   CTYPE * coord)
  {
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;

    if (offset_method == NO_OFFSET) {
      position_all_isov_centroid_multi_simple
        (scalar_grid, isodual_table, isovalue, isov_list, coord);
    }
    else if (offset_method == OFFSET_ONLY_MULTI_IN_CUBE) {
      const int dimension = scalar_grid.Dimension();
      IJK::ARRAY<CTYPE> temp_coord(dimension);
      IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);

      for (SIZE_TYPE kw = 0; kw < isov_list.size(); kw++) {
        const TABLE_INDEX table_index = isov_list[kw].table_index;

        if (isodual_table.NumIsoVertices(table_index) > 1) {
          // Apply offset.
          position_isov_in_cube_centroid_offset_multi_checkI
            (scalar_grid, isodual_table, isovalue, isov_list[kw],
             offset, cube, coord+kw*dimension, temp_coord.Ptr());
        }
        else {
          // No offset.
          position_isov_in_cube_centroid_multi_checkI
            (scalar_grid, isodual_table, isovalue, isov_list[kw],
             cube, coord+kw*dimension, temp_coord.Ptr());          
        }
      }
    }
    else {
      // offset_method == OFFSET_ALL
      position_all_isov_centroid_offset_multi_simple
        (scalar_grid, isodual_table, isovalue, isov_list, offset,
         coord);
    }
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices in centroid
   *    of isosurface-edge intersections.
   *    Simple version with offset argument. (C++ vector.)
   *  - Version with argument ISO_GRID_EDGE_INTERSECTION_OFFSET_METHOD.
   *  - Version using C++ STL vector for array coord[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE, typename CTYPE>
  void position_all_isov_centroid_offset_multi_simple
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const ISO_GRID_EDGE_INTERSECTION_OFFSET_METHOD offset_method,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(isov_list.size()*dimension);
    position_all_isov_centroid_offset_multi_simple
      (scalar_grid, isodual_table,isovalue, isov_list,
       offset_method, offset, IJK::vector2pointerNC(coord));
  }


  /*!
   *  @brief Position dual isosurface vertices in centroid
   *    of isosurface-edge intersections.    
   *  @param isov_list[] List of isosurface vertices
   *    including containing cube and associated isosurface patch.
   *  @param isopoly_vert[i*numv_per_isopoly+k]
   *    Index of k'th vertex of isosurface polytope i.
   *  @param isopoly_info[] Isososurface polygon information,
   *    including dual grid edge.
   *  @param[out] coord[] Array of vertex coordinates.
   *    - coord[i*dimension+d] = d'th coordinate of isosurface vertex i.
   *    @pre coord[] is pre-allocated to size at least
   *      isov_list.size()*scalar_grid.Dimension().
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename ISOV_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE, typename DUAL_ISOV_TYPE,
            typename CTYPE>
  void position_all_isov_centroid_multi
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   CTYPE * coord)
  {
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;

    const int dimension = scalar_grid.Dimension();
    const int num_vertices_per_isopoly =
      scalar_grid.NumFacetVertices();
    IJK::ARRAY<CTYPE> temp_coord(dimension);
    IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);
    IJK::PROCEDURE_ERROR
      error("position_all_isov_centroid_multi");

    if (!check_isopoly_vert_size
        (isopoly_vert, num_vertices_per_isopoly, isopoly_info.size(),
         error)) { throw error; }

    // num_incident_isopoly[i] =
    //   Number of isosurface polytopes incident on isov_list[i].
    std::vector<ISO_VERTEX_DEGREE>
      num_incident_isopoly(isov_list.size(), 0);

    // Set all coord entries to 0.0.
    for (SIZE_TYPE i = 0; i < isov_list.size()*dimension; i++)
      { coord[i] = 0.0; }

    SIZE_TYPE k = 0;
    for (SIZE_TYPE ipoly = 0; ipoly < isopoly_info.size(); ipoly++) {
      const VTYPE iend0 = isopoly_info[ipoly].GridEdgeEndpoint0();
      const int edge_direction = isopoly_info[ipoly].GridEdgeDirection();
      
      IJK::compute_isosurface_grid_edge_intersection_linear
        (scalar_grid, isovalue, iend0, edge_direction,
         temp_coord.Ptr());
      
      for (int j = 0; j < num_vertices_per_isopoly; j++) {
        const ISOV_INDEX_TYPE kw = isopoly_vert[k];
        CTYPE * coord_ptr = coord + dimension*kw;

        IJK::add_coord
          (dimension, coord_ptr, temp_coord.Ptr(), coord_ptr);
        num_incident_isopoly[kw]++;
        k++;
      }
    }

    if (k != isopoly_vert.size()) {
      error.AddMessage
        ("Programming error. Processed incorrect number of isopoly vertices.");
      error.AddMessage
        ("  Isosurface polytopes are incident on ",
         isopoly_vert.size(), " vertices.");
      error.AddMessage("  Processed ", k, " vertices.");
      throw error;
    }

    for (SIZE_TYPE kw = 0; kw < isov_list.size(); kw++) {
      const int ipatch = isov_list[kw].patch_index;
      const TABLE_INDEX table_index = isov_list[kw].table_index;
      const ISO_VERTEX_DEGREE num_incident_isopolyB =
        isodual_table.NumIncidentIsoPoly(table_index, ipatch);

      CTYPE * coord_ptr = coord + dimension*kw;

      if (num_incident_isopoly[kw] == num_incident_isopolyB) {
        IJK::multiply_coord
          (dimension, 1.0/num_incident_isopoly[kw], coord_ptr, coord_ptr);
      }
      else {
        // Isosurface vertex is in grid cube with some bipolar edges
        //   on the grid boundary.
        // Dual facets are missing from isopoly_info.
        // Use slower method to compute isosurface vertex position.
        position_isov_in_cube_centroid_multi_checkI
          (scalar_grid, isodual_table, isovalue, isov_list[kw], cube,
           coord_ptr, temp_coord.Ptr());
      }
    }
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices in centroid
   *    of isosurface-edge intersections. (C++ vector.)
   *  - Version using C++ STL vector for array coord[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename ISOV_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE, typename DUAL_ISOV_TYPE,
            typename CTYPE>
  void position_all_isov_centroid_multi
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   std::vector<CTYPE> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(isov_list.size()*dimension);
    position_all_isov_centroid_multi
      (scalar_grid, isodual_table,isovalue,
       isov_list, isopoly_vert, isopoly_info,
       IJK::vector2pointerNC(coord));

  }


  /*!
   *  @brief Position dual isosurface vertices in centroid
   *    of isosurface-edge intersections with offsets.
   *  - Offset isosurface-(grid edge) intersections so that all
   *    centroids are in cube interiors.
   *  @param isov_list[] List of isosurface vertices
   *    including containing cube and associated isosurface patch.
   *  @param isopoly_vert[i*numv_per_isopoly+k]
   *    Index of k'th vertex of isosurface polytope i.
   *  @param isopoly_info[] Isososurface polygon information,
   *    including dual grid edge.
   *  @param[out] coord[] Array of vertex coordinates.
   *    - coord[i*dimension+d] = d'th coordinate of isosurface vertex i.
   *    @pre coord[] is pre-allocated to size at least
   *      isov_list.size()*scalar_grid.Dimension().
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename ISOV_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE, typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE, typename CTYPE>
  void position_all_isov_centroid_offset_multi
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   const OFFSET_TYPE offset,
   CTYPE * coord)
  {
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;

    const int dimension = scalar_grid.Dimension();
    const int num_vertices_per_isopoly =
      scalar_grid.NumFacetVertices();
    IJK::ARRAY<CTYPE> temp_coord(dimension);
    IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);
    IJK::PROCEDURE_ERROR
      error("position_all_isov_centroid_offset_multi");

    if (!check_isopoly_vert_size
        (isopoly_vert, num_vertices_per_isopoly, isopoly_info.size(),
         error)) { throw error; }

    // num_incident_isopoly[i] =
    //   Number of isosurface polytopes incident on isov_list[i].
    std::vector<ISO_VERTEX_DEGREE>
      num_incident_isopoly(isov_list.size(), 0);

    // Set all coord entries to 0.0.
    for (SIZE_TYPE i = 0; i < isov_list.size()*dimension; i++)
      { coord[i] = 0.0; }

    SIZE_TYPE k = 0;
    for (SIZE_TYPE ipoly = 0; ipoly < isopoly_info.size(); ipoly++) {
      const VTYPE iend0 = isopoly_info[ipoly].GridEdgeEndpoint0();
      const int edge_direction = isopoly_info[ipoly].GridEdgeDirection();
      
      IJK::compute_isosurface_grid_edge_intersection_linear_offset
        (scalar_grid, isovalue, iend0, edge_direction, offset,
         temp_coord.Ptr());
      
      for (int j = 0; j < num_vertices_per_isopoly; j++) {
        const ISOV_INDEX_TYPE kw = isopoly_vert[k];
        CTYPE * coord_ptr = coord + dimension*kw;

        IJK::add_coord
          (dimension, coord_ptr, temp_coord.Ptr(), coord_ptr);
        num_incident_isopoly[kw]++;
        k++;
      }
    }

    if (k != isopoly_vert.size()) {
      error.AddMessage
        ("Programming error. Processed incorrect number of isopoly vertices.");
      error.AddMessage
        ("  Isosurface polytopes are incident on ",
         isopoly_vert.size(), " vertices.");
      error.AddMessage("  Processed ", k, " vertices.");
      throw error;
    }

    for (SIZE_TYPE kw = 0; kw < isov_list.size(); kw++) {
      const int ipatch = isov_list[kw].patch_index;
      const TABLE_INDEX table_index = isov_list[kw].table_index;
      const ISO_VERTEX_DEGREE num_incident_isopolyB =
        isodual_table.NumIncidentIsoPoly(table_index, ipatch);

      CTYPE * coord_ptr = coord + dimension*kw;

      if (num_incident_isopoly[kw] == num_incident_isopolyB) {
        IJK::multiply_coord
          (dimension, 1.0/num_incident_isopoly[kw], coord_ptr, coord_ptr);
      }
      else {
        // Isosurface vertex is in grid cube with some bipolar edges
        //   on the grid boundary.
        // Dual facets are missing from isopoly_info.
        // Use slower method to compute isosurface vertex position.
        position_isov_in_cube_centroid_offset_multi_checkI
          (scalar_grid, isodual_table, isovalue, isov_list[kw],
           offset, cube, coord_ptr, temp_coord.Ptr());
      }
    }
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices in centroid
   *    of isosurface-edge intersections with offsets. (C++ vector.)
   *  - Version using C++ STL vector for array coord[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename ISOV_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE, typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE, typename CTYPE>
  void position_all_isov_centroid_offset_multi
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(isov_list.size()*dimension);
    position_all_isov_centroid_offset_multi
      (scalar_grid, isodual_table,isovalue,
       isov_list, isopoly_vert, isopoly_info, offset,
       IJK::vector2pointerNC(coord));

  }


  /*!
   *  @brief Position dual isosurface vertices in centroid
   *    of isosurface-edge intersections. Offset in cubes with multi isov.
   *  - Offset isosurface-(grid edge) intersections in cubes
   *    containing multiple isosurface vertices.
   *  @param isov_list[] List of isosurface vertices
   *    including containing cube and associated isosurface patch.
   *  @param isopoly_vert[i*numv_per_isopoly+k]
   *    Index of k'th vertex of isosurface polytope i.
   *  @param isopoly_info[] Isososurface polygon information,
   *    including dual grid edge.
   *  @param[out] coord[] Array of vertex coordinates.
   *    - coord[i*dimension+d] = d'th coordinate of isosurface vertex i.
   *    @pre coord[] is pre-allocated to size at least
   *      isov_list.size()*scalar_grid.Dimension().
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename ISOV_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE, typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE, typename CTYPE>
  void position_all_isov_centroid_offset_only_multi
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   const OFFSET_TYPE offset,
   CTYPE * coord)
  {
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;

    const int dimension = scalar_grid.Dimension();
    const int num_vertices_per_isopoly =
      scalar_grid.NumFacetVertices();
    IJK::ARRAY<CTYPE> temp_coord(dimension);
    IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);
    IJK::PROCEDURE_ERROR
      error("position_all_isov_centroid_offset_only_multi");

    if (!check_isopoly_vert_size
        (isopoly_vert, num_vertices_per_isopoly, isopoly_info.size(),
         error)) { throw error; }

    // num_incident_isopoly[i] =
    //   Number of isosurface polytopes incident on isov_list[i].
    std::vector<ISO_VERTEX_DEGREE>
      num_incident_isopoly(isov_list.size(), 0);

    // Set all coord entries to 0.0.
    for (SIZE_TYPE i = 0; i < isov_list.size()*dimension; i++)
      { coord[i] = 0.0; }

    SIZE_TYPE k = 0;
    for (SIZE_TYPE ipoly = 0; ipoly < isopoly_info.size(); ipoly++) {
      const VTYPE iend0 = isopoly_info[ipoly].GridEdgeEndpoint0();
      const int edge_direction = isopoly_info[ipoly].GridEdgeDirection();

      // Default: No offset.
      IJK::compute_isosurface_grid_edge_intersection_linear
        (scalar_grid, isovalue, iend0, edge_direction,
         temp_coord.Ptr());
      
      for (int j = 0; j < num_vertices_per_isopoly; j++) {
        const ISOV_INDEX_TYPE kw = isopoly_vert[k];
        CTYPE * coord_ptr = coord + dimension*kw;

        IJK::add_coord
          (dimension, coord_ptr, temp_coord.Ptr(), coord_ptr);
        num_incident_isopoly[kw]++;
        k++;
      }
    }

    if (k != isopoly_vert.size()) {
      error.AddMessage
        ("Programming error. Processed incorrect number of isopoly vertices.");
      error.AddMessage
        ("  Isosurface polytopes are incident on ",
         isopoly_vert.size(), " vertices.");
      error.AddMessage("  Processed ", k, " vertices.");
      throw error;
    }

    for (SIZE_TYPE kw = 0; kw < isov_list.size(); kw++) {
      const int ipatch = isov_list[kw].patch_index;
      const TABLE_INDEX table_index = isov_list[kw].table_index;
      const ISO_VERTEX_DEGREE num_incident_isopolyB =
        isodual_table.NumIncidentIsoPoly(table_index, ipatch);

      CTYPE * coord_ptr = coord + dimension*kw;

      if (num_incident_isopoly[kw] == num_incident_isopolyB) {
        if (isodual_table.NumIsoVertices(table_index) == 1) {
          IJK::multiply_coord
            (dimension, 1.0/num_incident_isopoly[kw], coord_ptr, coord_ptr);
        }
        else {
          // Isosurface is in cube containing multiple iso vertices. Apply offset.
          position_isov_in_cube_centroid_offset_multi
          (scalar_grid, isodual_table, isovalue, isov_list[kw],
           offset, cube, coord_ptr, temp_coord.Ptr());
        }
      }
      else {
        // Isosurface vertex is in grid cube with some bipolar edges
        //   on the grid boundary.
        //   Dual facets are missing from isopoly_info.
        // Use slower method to compute isosurface vertex position.
        
        if (isodual_table.NumIsoVertices(table_index) == 1) {
          // No offset.
          position_isov_in_cube_centroid_multi
            (scalar_grid, isodual_table, isovalue, isov_list[kw],
             cube, coord_ptr, temp_coord.Ptr());          
        }
        else {
          // Apply offset.
          position_isov_in_cube_centroid_offset_multi
            (scalar_grid, isodual_table, isovalue, isov_list[kw],
             offset, cube, coord_ptr, temp_coord.Ptr());
        }
      }
    }
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices in centroid
   *    of isosurface-edge intersections. (Offset argument.)
   *    - Version with argument ISO_GRID_EDGE_INTERSECTION_OFFSET_METHOD.
   */
  template <typename GRID_TYPE,typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename DUAL_ISOV_TYPE,
            typename ISOV_INDEX_TYPE, typename ISOPOLY_INFO_TYPE, 
            typename OFFSET_TYPE,
            typename CTYPE>
  void position_all_isov_centroid_offset_multi
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   const ISO_GRID_EDGE_INTERSECTION_OFFSET_METHOD offset_method,
   const OFFSET_TYPE offset,
   CTYPE * coord)
  {
    if (offset_method == NO_OFFSET) {
      position_all_isov_centroid_multi
        (scalar_grid, isodual_table, isovalue, isov_list,
         isopoly_vert, isopoly_info, coord);
    }
    else if (offset_method == OFFSET_ONLY_MULTI_IN_CUBE) {
      position_all_isov_centroid_offset_only_multi
        (scalar_grid, isodual_table, isovalue, isov_list,
         isopoly_vert, isopoly_info, offset, coord);
    }
    else {
      // offset_method == OFFSET_ALL
      position_all_isov_centroid_offset_multi
        (scalar_grid, isodual_table, isovalue, isov_list,
         isopoly_vert, isopoly_info, offset,         coord);
    }
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices in centroid
   *    of isosurface-edge intersections. (Offset argument. C++ vector.)
   *  - Version with argument ISO_GRID_EDGE_INTERSECTION_OFFSET_METHOD.
   *  - Version using C++ STL vector for coord.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename DUAL_ISOV_TYPE,
            typename ISOV_INDEX_TYPE, typename ISOPOLY_INFO_TYPE, 
            typename OFFSET_TYPE, typename CTYPE>
  void position_all_isov_centroid_offset_multi
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   const ISO_GRID_EDGE_INTERSECTION_OFFSET_METHOD offset_method,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(isov_list.size()*dimension);
    position_all_isov_centroid_offset_multi
      (scalar_grid, isodual_table, isovalue, isov_list,
       isopoly_vert, isopoly_info, offset_method, offset,
       IJK::vector2pointerNC(coord));
  }


  /*!
   *  @brief Position dual isosurface vertices near cube centers.
   *  - More than one vertex can be in a cube.
   *  - If cube contains multiple isosurface then vertices are positioned
   *    near but not on cube center.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename DUAL_ISOV_TYPE, typename CTYPE0, typename CTYPE1>
  void position_all_isov_near_cube_center_multi
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const CTYPE0 offset,
   CTYPE1 * coord)
  {
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;

    const int dimension = scalar_grid.Dimension();

    IJK::ARRAY<CTYPE1> vcoord(dimension);

    IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);
    IJK::UNIT_CUBE<int,int,int> unit_cube(dimension);
    bool intersects_facet[2*dimension];

    for (SIZE_TYPE i = 0; i < isov_list.size(); i++) {
      const VTYPE icube = isov_list[i].cube_index;
      const NTYPE ipatch = isov_list[i].patch_index;
      const TABLE_INDEX it = isov_list[i].table_index;

      if (isodual_table.NumIsoVertices(it) == 1) {
        scalar_grid.ComputeCubeCenterCoord(icube, coord+i*dimension);
      }
      else {

        scalar_grid.ComputeCubeCenterCoord(icube, vcoord.Ptr());

        // Set intersects facet to false.
        for (int d = 0; d < dimension; d++) {
          intersects_facet[2*d] = false;
          intersects_facet[2*d+1] = false;
        }

        for (NTYPE ie = 0; ie < cube.NumEdges(); ie++) {
          if (isodual_table.IsBipolar(it, ie)) {
            if (isodual_table.IncidentIsoVertex(it, ie) == ipatch) {
              NTYPE k0 = cube.EdgeEndpoint(ie, 0);
              NTYPE k1 = cube.EdgeEndpoint(ie, 1);

              for (int d = 0; d < dimension; d++) {
                NTYPE c = unit_cube.VertexCoord(k0, d);
                if (c == unit_cube.VertexCoord(k1, d)) {
                  if (c > 0)
                    { intersects_facet[2*d+1] = true; }
                  else
                    { intersects_facet[2*d] = true; }
                }
              }
            }
          }
        }

        for (int d = 0; d < dimension; d++) {

          if (intersects_facet[2*d] != intersects_facet[2*d+1]) {
            if (intersects_facet[2*d])
              { vcoord[d] -= offset; }
            else
              { vcoord[d] += offset; }
          }
        }
        IJK::copy_coord(dimension, vcoord.PtrConst(), coord+i*dimension);
      }
    }
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices near cube centers.
   *    (C++ vector)
   *  - More than one vertex can be in a cube.
   *  - C++ STL vector format for array coord[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename DUAL_ISOV_TYPE, typename CTYPE0, typename CTYPE1>
  void position_all_isov_near_cube_center_multi
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const CTYPE0 offset,
   std::vector<CTYPE1> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(isov_list.size()*dimension);
    position_all_isov_near_cube_center_multi
      (scalar_grid, isodual_table, isov_list, offset,
       IJK::vector2pointerNC(coord));
  }

  //@}


  // *****************************************************************
  //! @name Position lifted interval volume vertices.
  // *****************************************************************

  //@{

  /*!
   *  @brief Position interval volume vertices which have been lifted
   *    to one higher dimension.
   *  @param isovalue0 Lower isovalue.
   *  @param isovalue1 Upper isovalue.
   *  @pre isovalue0 < isovalue1.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISOVAL0_TYPE, typename ISOVAL1_TYPE,
            typename DUAL_ISOV_TYPE, typename CTYPE>
  void position_all_isov_ivol_lifted
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISOVAL0_TYPE isovalue0,
   const ISOVAL1_TYPE isovalue1,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   CTYPE * coord)
  {
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const int dimension = scalar_grid.Dimension();

    if (dimension < 1) { return; }

    const VTYPE numv_in_grid_facet_maxd =
      scalar_grid.AxisIncrement(dimension-1);
    IJK::ARRAY<CTYPE> temp_coord(dimension);
    IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);

    for (SIZE_TYPE isov = 0; isov < isov_list.size(); isov++) {

      const VTYPE icube = isov_list[isov].cube_index;

      if (icube < numv_in_grid_facet_maxd) {
        // Apply isovalue1 to (*,*,*,0) grid cubes.
        position_isov_in_cube_centroid_multi
          (scalar_grid, isodual_table, isovalue1, isov_list[isov], cube,
           coord+isov*dimension, temp_coord.Ptr());
      }
      else {
        // Apply isovalue0 to (*,*,*,1) grid cubes.
        position_isov_in_cube_centroid_multi
          (scalar_grid, isodual_table, isovalue0, isov_list[isov], cube,
           coord+isov*dimension, temp_coord.Ptr());

      }
    }
  }


  /*!
   *  @overload
   *  @brief Position interval volume vertices which have been lifted
   *    to one higher dimension. (C++ vector.)
   *  - Version using C++ STL vector for array coord[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE0, typename STYPE1,
            typename DUAL_ISOV_TYPE, typename CTYPE>
  void position_all_isov_ivol_lifted
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE0 isovalue0,
   const STYPE1 isovalue1,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   std::vector<CTYPE> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(isov_list.size()*dimension);
    position_all_isov_ivol_lifted
      (scalar_grid, isodual_table, isovalue0, isovalue1, isov_list,
       IJK::vector2pointerNC(coord));
  }

  //@}


  // *****************************************************************
  //! @name Reposition isosurface vertices away from cube facets.
  // *****************************************************************


  /*!
   *  @brief Reposition isosurface vertices away from cube facets.
   *  - Note: Forces isosurface vertex to lie in interior
   *    of associated cube.
   *  @param offset Offset. Must be in range [0,0.5].
   */
  template <typename GRID_TYPE,
            typename ISOV_INDEX_TYPE, typename CUBE_INDEX_TYPE,
            typename OFFSET_TYPE, typename CTYPE,
            typename NTYPE>
  void reposition_isov_away_from_cube_facets
  (const GRID_TYPE & grid,
   const ISOV_INDEX_TYPE isov,
   const CUBE_INDEX_TYPE icube,
   const OFFSET_TYPE offset,
   CTYPE * coord,
   NTYPE & num_repositions)
  {
    const int dimension = grid.Dimension();
    IJK::ARRAY<CTYPE> cube_coord(dimension);

    // Initialize
    num_repositions = 0;

    grid.ComputeCoord(icube, cube_coord.Ptr());
    CTYPE * isov_coord = coord + dimension*isov;

    for (int d = 0; d < dimension; d++) {
      const CTYPE minc = cube_coord[d] + offset;
      const CTYPE maxc = cube_coord[d] + 1 - offset;
      if (isov_coord[d] < minc) {
        isov_coord[d] = minc;
        num_repositions++;
      }
      else if (isov_coord[d] > maxc) {
        isov_coord[d] = maxc;
        num_repositions++;
      }
    }
  }


  /*!
   *  @brief Reposition all isosurface vertices away from cube facets.
   *  - Note: Forces each isosurface vertices to lie in associated cube.
   *  @param offset Offset. Must be in range [0,0.5].
   *  @param num_moved[out] Number of isosurface vertices
   *    moved away from some cube facet.
   *    - Note: If an isosurface vertex is moved away from more
   *      than one cube facet, it is still only counted once in num_moved.
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename NTYPE, typename OFFSET_TYPE,
            typename CTYPE, typename NTYPEM>
  void reposition_all_isov_away_from_cube_facets
  (const GRID_TYPE & grid,
   const DUAL_ISOV_TYPE isov_list[],
   const NTYPE num_isov,
   const OFFSET_TYPE offset,
   CTYPE * coord,
   NTYPEM & num_moved)
  {
    typedef typename DUAL_ISOV_TYPE::CUBE_INDEX_TYPE CUBE_INDEX_TYPE;

    // Initialize
    num_moved = 0;

    for (NTYPE isov = 0; isov < num_isov; isov++) {
      const CUBE_INDEX_TYPE icube = isov_list[isov].cube_index;
      NTYPEM num_repositions_i;
      reposition_isov_away_from_cube_facets
        (grid, isov, icube, offset, coord, num_repositions_i);
      if (num_repositions_i > 0)
        { num_moved++; }
    }
  }


  /*!
   *  @overload
   *  @brief Reposition all isosurface vertices away from cube facets.
   *    (Specialization of isov_list[] to int.)
   *  - Specialization for isov_list[] being a list of cube indices.
   *  - Note: Forces each isosurface vertices to lie in associated cube.
   *  @param isov_list[isov] Index of cube containing isosurface vertex isov.
   *  @param num_moved[out] Number of isosurface vertices
   *    moved away from some cube facet.
   *    - Note: If an isosurface vertex is moved away from more
   *      than one cube facet, it is still only counted once in num_moved.
   *  @param offset Offset. Must be in range [0,0.5].
   */
  template <typename GRID_TYPE, typename NTYPE, typename OFFSET_TYPE,
            typename CTYPE, typename NTYPEM>
  void reposition_all_isov_away_from_cube_facets
  (const GRID_TYPE & grid,
   const int isov_list[],
   const NTYPE num_isov,
   const OFFSET_TYPE offset,
   CTYPE * coord,
   NTYPEM & num_moved)
  {
    // Initialize
    num_moved = 0;

    for (NTYPE isov = 0; isov < num_isov; isov++) {
      NTYPEM num_repositions_i;
      reposition_isov_away_from_cube_facets
        (grid, isov, isov_list[isov], offset, coord, num_repositions_i);
      if (num_repositions_i > 0)
        { num_moved++; }
    }
  }


  /*!
   *  @overload
   *  @brief Reposition all isosurface vertices away from cube facets.
   *    (C++ STL vector for isov_list[].)
   *  - Version with C++ STL vector format for array isov_list[].
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE, typename CTYPE, typename NTYPE>
  void reposition_all_isov_away_from_cube_facets
  (const GRID_TYPE & grid,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const OFFSET_TYPE offset,
   CTYPE * coord,
   NTYPE & num_moved)
  {
    reposition_all_isov_away_from_cube_facets
      (grid, IJK::vector2pointer(isov_list), isov_list.size(),
       offset, coord, num_moved);
  }


  /*!
   *  @overload
   *  @brief Reposition all isosurface vertices away from cube facets.
   *    (C++ STL vector for isov_list[] and coord[].)
   *  - Version with C++ STL vector format for array isov_list[].
   *  - Version with C++ STL vector format for array coord[].
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE, typename CTYPE, typename NTYPE>
  void reposition_all_isov_away_from_cube_facets
  (const GRID_TYPE & grid,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord,
   NTYPE & num_moved)
  {
    reposition_all_isov_away_from_cube_facets
      (grid, isov_list, offset, IJK::vector2pointerNC(coord),
       num_moved);
  }


  // ******************************************************************
  //! @name Reposition isosurface vertices away from cube ridges.
  // ******************************************************************

  /*!
   *  @brief Reposition isosurface vertices away from cube ridges.
   *  - Return true if isosurface vertex isov is moved.
   *  - Note: Forces isosurface vertex to lie in interior
   *    of associated cube.
   *  - In 3D, cube ridges are cube edges.
   *  @param offset Offset. Must be in range [0,0.5].
   */
  template <typename GRID_TYPE,
            typename ISOV_INDEX_TYPE, typename CUBE_INDEX_TYPE,
            typename OFFSET_TYPE, typename CTYPE>
  bool reposition_isov_away_from_cube_ridges
  (const GRID_TYPE & grid,
   const ISOV_INDEX_TYPE isov,
   const CUBE_INDEX_TYPE icube,
   const OFFSET_TYPE offset,
   CTYPE * coord)
  {
    const int dimension = grid.Dimension();
    const int TWO(2);
    IJK::ARRAY<CTYPE> cube_coord(dimension);

    grid.ComputeCoord(icube, cube_coord.Ptr());
    CTYPE * isov_coord = coord + dimension*isov;

    int num_close_facets = 0;
    
    // Count number of facets close to isosurface vertex.
    for (int d = 0; d < dimension; d++) {
      const CTYPE minc = cube_coord[d] + offset;
      const CTYPE maxc = cube_coord[d] + 1 - offset;
      if (isov_coord[d] < minc)
        { num_close_facets++; }
      else if (isov_coord[d] > maxc) {
        { num_close_facets++; }
      }
    }

    if (num_close_facets >= TWO) {
      int num_repositions = 0;

      reposition_isov_away_from_cube_facets
        (grid, isov, icube, offset, coord, num_repositions);
      
      if (num_repositions > 0) {
        return true;
      }
    }

    return false;
  }


  /*!
   *  @brief Reposition all isosurface vertices away from cube ridges.
   *  - Note: Forces each isosurface vertices to lie in associated cube.
   *  @param offset Offset. Must be in range [0,0.5].
   *  - In 3D, cube ridges are cube edges.
   *  @param num_moved[out] Number of isosurface vertices
   *    moved away from some cube ridge.
   *    - Note: If an isosurface vertex is moved away from more
   *      than one cube ridge, it is still only counted once in num_moved.
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename NTYPE, typename OFFSET_TYPE,
            typename CTYPE, typename NTYPEM>
  void reposition_all_isov_away_from_cube_ridges
  (const GRID_TYPE & grid,
   const DUAL_ISOV_TYPE isov_list[],
   const NTYPE num_isov,
   const OFFSET_TYPE offset,
   CTYPE * coord,
   NTYPEM & num_moved)
  {
    typedef typename DUAL_ISOV_TYPE::CUBE_INDEX_TYPE CUBE_INDEX_TYPE;

    // Initialize
    num_moved = 0;

    for (NTYPE isov = 0; isov < num_isov; isov++) {
      const CUBE_INDEX_TYPE icube = isov_list[isov].cube_index;
      if (reposition_isov_away_from_cube_ridges
          (grid, isov, icube, offset, coord))
        { num_moved++; }
    }
  }


  /*!
   *  @overload
   *  @brief Reposition all isosurface vertices away from cube ridges.
   *    (Specialization of isov_list[] to int.)
   *  - Specialization for isov_list[] being a list of cube indices.
   *  - Note: Forces each isosurface vertices to lie in associated cube.
   *  - In 3D, cube ridges are cube edges.
   *  @param isov_list[isov] Index of cube containing isosurface vertex isov.
   *  @param offset Offset. Must be in range [0,0.5].
   */
  template <typename GRID_TYPE, typename NTYPE, typename OFFSET_TYPE,
            typename CTYPE, typename NTYPEM>
  void reposition_all_isov_away_from_cube_ridges
  (const GRID_TYPE & grid,
   const int isov_list[],
   const NTYPE num_isov,
   const OFFSET_TYPE offset,
   CTYPE * coord,
   NTYPEM & num_moved)
  {
    // Initialize
    num_moved = 0;

    for (NTYPE isov = 0; isov < num_isov; isov++) {
      if (reposition_isov_away_from_cube_ridges
          (grid, isov, isov_list[isov], offset, coord))
        { num_moved++; }
    }
  }


  /*!
   *  @overload
   *  @brief Reposition all isosurface vertices away from cube ridges.
   *    (C++ STL vector for isov_list[].)
   *  - Version with C++ STL vector format for array isov_list[].
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE, typename CTYPE, typename NTYPE>
  void reposition_all_isov_away_from_cube_ridges
  (const GRID_TYPE & grid,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const OFFSET_TYPE offset,
   CTYPE * coord,
   NTYPE & num_moved)
  {
    reposition_all_isov_away_from_cube_ridges
      (grid, IJK::vector2pointer(isov_list), isov_list.size(),
       offset, coord, num_moved);
  }


  /*!
   *  @overload
   *  @brief Reposition all isosurface vertices away from cube ridges.
   *    (C++ STL vector for isov_list[] and coord[].)
   *  - Version with C++ STL vector format for array isov_list[].
   *  - Version with C++ STL vector format for array coord[].
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE, typename CTYPE, typename NTYPE>
  void reposition_all_isov_away_from_cube_ridges
  (const GRID_TYPE & grid,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord,
   NTYPE & num_moved)
  {
    reposition_all_isov_away_from_cube_ridges
      (grid, isov_list, offset, IJK::vector2pointerNC(coord),
       num_moved);
  }
  

  // ******************************************************************
  //! @name Separate isosurface vertices near cube facets.
  // ******************************************************************

  /*!
   *  @brief Store information about grid face, isosurface vertex pair.
   *  - Note: Face, NOT facet.
   *  @tparam BITSET_TYPE Bitset used to indicate position
   *    of facet with respect to cube containing isosurface vertex.
   *    @pre BITSET_TYPE::size() >= Grid dimension.
   */
  template <typename BITSET_TYPE,
            typename FACEV_TYPE, typename ISOV_TYPE, typename DIST_TYPE>
  class FACE_ISOV_PAIR_BASE {

  public:

    // Bitset type indicating positition of d'th facet containing the face
    //   relative to the cube containing the isosurface vertex.
    typedef BITSET_TYPE FACET_POS_BITSET_TYPE;


  protected:

    /// @brief Set face_vertex, isosurface vertex, and distance.
    void _Set
    (const FACEV_TYPE face_vertex,
     const ISOV_TYPE isov, 
     const DIST_TYPE distance)
    {
      this->face_vertex = face_vertex;
      this->isov = isov;
      this->distance = distance;
    }
    
    /// @brief Constructor that does not initialize bitset.
    /// - Only for use by derived classes.
    FACE_ISOV_PAIR_BASE
    (const FACEV_TYPE face_vertex, const ISOV_TYPE isov,
     const DIST_TYPE distance)
    {
      _Set(face_vertex, isov, distance);
    }


  public:

    /// @brief Index of lower/leftmost vertex in face.
    FACEV_TYPE face_vertex;
    
    /*!
     *  @brief isov Isosurface vertex.
     *  - Note: isov_list[isov] contains index of cube containing isov.
     */
    ISOV_TYPE isov;

    /// @brief Distance from isov to facet
    DIST_TYPE distance;

    /*!
     *  @brief Positition of d'th facet containing the face
     *     relative to the cube containing the isosurface vertex.
     *  - Facet is orthogonal to direction d.
     *  - facet_pos[d] = 0, if facet is below/left of cube.
     *  - facet_pos[d] = 1, if facet is above/right of cube.
     */
    FACET_POS_BITSET_TYPE facet_pos;
    

    void Set
    (const FACEV_TYPE face_vertex,
     const ISOV_TYPE isov, const FACET_POS_BITSET_TYPE flag_pos_bitset,
     const DIST_TYPE distance)
    {
      _Set(face_vertex, isov, distance);
      facet_pos = flag_pos_bitset;
    }

    /// @brief Constructor.
    FACE_ISOV_PAIR_BASE
    (const FACEV_TYPE face_vertex, const ISOV_TYPE isov,
     const FACET_POS_BITSET_TYPE facet_pos_bitset,const DIST_TYPE distance)
    {
      Set(face_vertex, isov, facet_pos_bitset, distance);
    }

  };


  /*!
   *  @overload
   *  @brief Store information about grid face, isosurface vertex pair.
   *  - Version defined by bitset size, not bitset type.
   *  @tparam BITSET_SIZE Number of bits in bitset.
   *    - Usually 8.
   *    @pre BITSET_SIZE >= Grid dimension.
   */
  template <const int BITSET_SIZE,
            typename FACEV_TYPE, typename ISOV_TYPE, typename DIST_TYPE>
  using FACE_ISOV_PAIR =
    FACE_ISOV_PAIR_BASE<std::bitset<BITSET_SIZE>, FACEV_TYPE, ISOV_TYPE, DIST_TYPE>;
  

  /*!
   *  @brief Store information about grid facet, isosurface vertex pair.
   *  - Note: Facet, NOT face.
   *  - Note: Bitset size in FACE_ISOV_PAIR<> is 8,
   *      since only 1 bit is actually needed to store above/below.
   */
  template <typename FACETV_TYPE, typename ISOV_TYPE, typename DIST_TYPE>
  class FACET_ISOV_PAIR:
    public FACE_ISOV_PAIR<8, FACETV_TYPE, ISOV_TYPE, DIST_TYPE> {
    
  public:

    void Set
    (const FACETV_TYPE facet_vertex,
     const ISOV_TYPE isov, const int iside,
     const DIST_TYPE distance)
    {
      this->_Set(facet_vertex, isov, distance);
      this->facet_pos[0] = iside;
    }

    /// Constructor.
    FACET_ISOV_PAIR
    (const FACETV_TYPE facet_vertex, const ISOV_TYPE isov,
     const int iside, const DIST_TYPE distance):
      FACE_ISOV_PAIR<8,FACETV_TYPE,ISOV_TYPE,DIST_TYPE>(facet_vertex,isov,distance)
    {
      Set(facet_vertex, isov, iside, distance);
    }

    /// Return side of cube containing facet.
    int CubeSideContainingFacet() const
    { return this->facet_pos[0]; }
  };


  // local
  namespace {

    // If flag_moved[i] is false, then set to true
    //   and increment num_moved.
    void record_moved(const int i, std::vector<bool> & flag_moved,
                      int & num_moved)
    {
      if (!flag_moved[i]) {
        flag_moved[i] = true;
        num_moved++;
      }
    }


    /*!
     *  @brief Add facet isov pair.
     *  @param facet_index Facet index. Index of lowest/leftmost facet vertex.
     *    - Facet orthogonal direction is implicit.
     *  @param iside Side (0 or 1) of cube containing facet.
     */
    template <typename FACET_INDEX_TYPE, typename ISOV_TYPE,
              typename CTYPE,
              typename FACET_ISOV_PAIR_TYPE,
              typename VNEAR_BITSET_TYPE>
    void add_facet_isov_pair
    (const FACET_INDEX_TYPE facet_index, const ISOV_TYPE isov,
     const CTYPE cdiff, const int iside,
     std::vector<FACET_ISOV_PAIR_TYPE> & facet_isov_list,
     std::unordered_map<FACET_INDEX_TYPE, VNEAR_BITSET_TYPE> & vnear_facet)
    {
      const FACET_ISOV_PAIR_TYPE
        facet_isov_pair(facet_index, isov, iside, cdiff);
      facet_isov_list.push_back(facet_isov_pair);

      IJK::insert_or_update_hash_table_bitset
        (facet_index, iside, 1, vnear_facet);
    }

    
    // Throw error if facet not found.
    template <typename ITER_TYPE, typename MAP_TYPE, typename VTYPE>
    void throw_error_if_facet_not_found
    (const ITER_TYPE facet_iter, const MAP_TYPE & vnear_facet,
     const int orth_dir, const VTYPE facet_vertex,
     IJK::ERROR & error)
    {
      if (facet_iter == vnear_facet.end()) {
        error.AddMessage
          ("Programming error. Facet does not appear in map vnear_facet[].");
        error.AddMessage("  Facet orthogonal direction: ", orth_dir);
        error.AddMessage("  Lower/leftmost facet vertex: ", facet_vertex);
        throw error;
      }
    }
    
  };

  
  /*!
   *  @brief Separate isosurface vertices that are near shared grid cube facets.
   *  - Separate isosurface vertices w0, w1 that are from two different
   *    grid cubes that intersect at a facet.
   *  - Move both vertices away from the facet.
   *  @param min_distance Move vertices if both cubes have
   *    isosurface vertices within distance (Linf) \a min_distance
   *    of facet f.
   *  @param offset Move vertices so that they are distance
   *    \a offset from the facet f.
   *  @param[out] num_repositions Number of times isosurface vertices
   *    are repositioned away from some grid cube facet.
   *    - Note: An isosurface vertex that is repositioned away
   *      from multiple grid cube facets will contribute
   *      multiple times to this count.
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename MIN_DIST_TYPE, typename OFFSET_TYPE,
            typename CTYPE, typename NTYPEV, typename NTYPEM>
  void separate_isov_near_shared_grid_cube_facets
  (const GRID_TYPE & grid,
   const DUAL_ISOV_TYPE isov_list[],
   const NTYPEV num_isov,
   const MIN_DIST_TYPE min_distance,
   const OFFSET_TYPE offset,
   CTYPE * coord,
   NTYPEM & num_moved)
  {
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE CUBE_INDEX_TYPE;
    typedef FACET_ISOV_PAIR<VERTEX_INDEX_TYPE,NTYPEV,CTYPE>
      FACET_ISOV_PAIR_TYPE;
    typedef typename std::vector<FACET_ISOV_PAIR_TYPE>::size_type SIZE_TYPE;

    // Bitset type indicating vertices near facet.
    typedef std::bitset<8> VNEAR_BITSET_TYPE;

    typedef std::unordered_map<VERTEX_INDEX_TYPE, VNEAR_BITSET_TYPE>
      VNEAR_FACET_MAP;

    const int dimension = grid.Dimension();
    VNEAR_FACET_MAP vnear_facet;
    std::vector<FACET_ISOV_PAIR_TYPE> facet_isov_list;
    IJK::PROCEDURE_ERROR error
      ("separate_isov_near_shared_grid_cube_facets");


    // Array used to track which isosurface vertices have moved.
    // Really only needed to get an accurate count of num_moved, but...
    std::vector<bool> flag_moved(num_isov, false);
    
    // Initialize
    num_moved = 0;

    for (int d = 0; d < dimension; d++ ) {

      vnear_facet.clear();
      facet_isov_list.clear();

      for (NTYPEV isov = 0; isov < num_isov; isov++) {
        const CUBE_INDEX_TYPE icube =
          get_isovert_grid_cube_index(isov_list[isov]);

        const CTYPE isov_coord_d = coord[dimension*isov+d];
        const CTYPE cube_coord_d = grid.CoordD(icube, d);
        const CTYPE cdiff0 = isov_coord_d - cube_coord_d;
        const CTYPE cdiff1 = (cube_coord_d+1) - isov_coord_d;

        if (cdiff0 < min_distance) {
          add_facet_isov_pair
            (icube, isov, cdiff0, 0,
             facet_isov_list, vnear_facet);
        }

        if (cdiff1 < min_distance) {
          const VERTEX_INDEX_TYPE iv_primary =
            grid.NextVertex(icube, d);

          add_facet_isov_pair
            (iv_primary, isov, cdiff1, 1,
             facet_isov_list, vnear_facet);          
        }
      }


      const VNEAR_BITSET_TYPE BITSET_11 = 3;
      for (SIZE_TYPE i = 0; i < facet_isov_list.size(); i++) {
        const VERTEX_INDEX_TYPE facet_vertex =
          facet_isov_list[i].face_vertex;
        const NTYPEV isov = facet_isov_list[i].isov;
        const int iside =
          facet_isov_list[i].CubeSideContainingFacet();

        const CTYPE facet_coord_d = grid.CoordD(facet_vertex, d);

        auto facet_iter = vnear_facet.find(facet_vertex);
        throw_error_if_facet_not_found
          (facet_iter, vnear_facet, d, facet_vertex, error);

        const VNEAR_BITSET_TYPE vnear_bitset = facet_iter->second;
        
        // Move isov away from facet.
        if (vnear_bitset == BITSET_11) {

          if (iside == 0) {
            // Facet is below isov and above some other close iso vertex.
            coord[dimension*isov+d] = facet_coord_d + offset;
          }
          else {
            // Facet is above isov and below some other close iso vertex.
            coord[dimension*isov+d] = facet_coord_d - offset;
            record_moved(isov, flag_moved, num_moved);
          }

          record_moved(isov, flag_moved, num_moved);
        }
      }
    }
  }


  /*!
   *  @overload
   *  @brief Separate isosurface vertices that are near shared grid cube facets.
   *    (C++ vector.)
   *  - Version using C++ STL vector for isov_list[] and coord[].
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename MIN_DIST_TYPE, typename OFFSET_TYPE,
            typename CTYPE, typename NTYPEM>
  void separate_isov_near_shared_grid_cube_facets
  (const GRID_TYPE & grid,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const MIN_DIST_TYPE min_distance,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord,
   NTYPEM & num_moved)
  {
    separate_isov_near_shared_grid_cube_facets
    (grid, IJK::vector2pointer(isov_list), isov_list.size(),
     min_distance, offset, IJK::vector2pointerNC(coord), num_moved);
  }


  // ******************************************************************
  //! @name Separate isosurface vertices that are near grid ridges.
  // ******************************************************************
  
  // local
  namespace {

    template <typename CTYPE, typename CUBE_CTYPE,
              typename MIN_DISTANCE_TYPE>
    bool is_isov_close_to_facet
    (const CTYPE isov_coord[],
     const CUBE_CTYPE cube_coord[],
     const int orth_dir,
     const int iside,
     const MIN_DISTANCE_TYPE min_distance)
    {
      if (iside == 0) {
        if (isov_coord[orth_dir] < cube_coord[orth_dir] + min_distance)
          { return true; }
      }
      else {
        if (isov_coord[orth_dir] > cube_coord[orth_dir] + 1 - min_distance)
          { return true; }
      }

      return false;
    }
    
    
    template <typename CTYPE, typename CUBE_CTYPE,
              typename MIN_DISTANCE_TYPE>
    bool is_isov_close_to_ridge
    (const CTYPE isov_coord[],
     const CUBE_CTYPE cube_coord[],
     const int orth_dir0, const int orth_dir1,
     const int iside0, const int iside1,
     const MIN_DISTANCE_TYPE min_distance)
    {
      if (!is_isov_close_to_facet
          (isov_coord, cube_coord, orth_dir0, iside0, min_distance))
        { return false; }

      if (!is_isov_close_to_facet
          (isov_coord, cube_coord, orth_dir1, iside1, min_distance))
        { return false; }      

      return true;
    }


    template <typename CTYPE, typename CUBE_CTYPE>
    CTYPE distance_to_facet
    (const CTYPE isov_coord[],
     const CUBE_CTYPE cube_coord[],
     const int orth_dir, const int iside)
    {
      if (iside == 0) {
        return (isov_coord[orth_dir] - cube_coord[orth_dir]);
      }
      else {
        return (cube_coord[orth_dir]+1 - isov_coord[orth_dir]);
      }
    }


    // L-infinity distance from isov to closest point on ridge.
    template <typename CTYPE, typename CUBE_CTYPE>
    CTYPE distance_to_ridge
    (const CTYPE isov_coord[],
     const CUBE_CTYPE cube_coord[],
     const int orth_dir0, const int orth_dir1,
     const int iside0, const int iside1)
    {
      const CTYPE cdiff0 =
        distance_to_facet(isov_coord, cube_coord, orth_dir0, iside0);
      const CTYPE cdiff1 =
        distance_to_facet(isov_coord, cube_coord, orth_dir1, iside1);

      const CTYPE cdiff = std::min(cdiff0, cdiff1);
      
      return cdiff;
    }
      
    
    
    template <typename GRID_TYPE, typename ITYPE>
    ITYPE get_ridge_primary_vertex
    (const GRID_TYPE & grid, const ITYPE & icube,
     const int orth_dir0, const int orth_dir1,
     const int iside0, const int iside1)
    {
      ITYPE iv_primary = icube;
      
      if (iside0 == 1) {
        iv_primary = grid.NextVertex(iv_primary, orth_dir0);
      }

      if (iside1 == 1) {
        iv_primary = grid.NextVertex(iv_primary, orth_dir1);
      }

      return iv_primary;
    }

    
    // Add ridge-isov pair.
    template <typename RIDGE_INDEX_TYPE, typename ISOV_TYPE,
              typename CTYPE,
              typename FACE_ISOV_PAIR_TYPE,
              typename VNEAR_BITSET_TYPE>
    void add_ridge_isov_pair
    (const RIDGE_INDEX_TYPE facet_index, const ISOV_TYPE isov,
     const CTYPE cdiff, const int iside0, const int iside1,
     std::vector<FACE_ISOV_PAIR_TYPE> & ridge_isov_list,
     std::unordered_map<RIDGE_INDEX_TYPE, VNEAR_BITSET_TYPE> & vnear_ridge)
    {
      typedef typename FACE_ISOV_PAIR_TYPE::FACET_POS_BITSET_TYPE
        FACET_POS_BITSET_TYPE;
      FACET_POS_BITSET_TYPE facet_pos_bitset;

      facet_pos_bitset[0] = iside0;
      facet_pos_bitset[1] = iside1;
      
      const FACE_ISOV_PAIR_TYPE
        ridge_isov_pair(facet_index, isov, facet_pos_bitset, cdiff);

      ridge_isov_list.push_back(ridge_isov_pair);

      IJK::insert_or_update_hash_table_bitset
        (facet_index, iside0, 1, vnear_ridge);
      IJK::insert_or_update_hash_table_bitset
        (facet_index, 2+iside1, 1, vnear_ridge);      
    }

    
    template <typename BITSET_TYPE>
    bool are_both_bits1
    (const BITSET_TYPE & bitset, const int i0, const int i1)
    {
      if ((bitset[i0] == 1) && (bitset[i1] == 1))
        { return true; }
      else
        { return false; }
    }


    template <typename CTYPEV, typename OFFSET_TYPE,
              typename CTYPE>
    void move_coord_away_from_facet
    (const CTYPEV gridv_coord[], const int orth_dir, const int iside,
     const OFFSET_TYPE offset, CTYPE * isov_coord)
    {
      if (iside == 0) {
        isov_coord[orth_dir] = gridv_coord[orth_dir] + offset;
      }
      else {
        isov_coord[orth_dir] = gridv_coord[orth_dir] - offset;
      }
    }

    
    // Throw error if ridge not found.
    template <typename ITER_TYPE, typename MAP_TYPE, typename VTYPE>
    void throw_error_if_ridge_not_found
    (const ITER_TYPE facet_iter, const MAP_TYPE & vnear_facet,
     const int orth_dir0, const int orth_dir1,
     const VTYPE facet_vertex, IJK::ERROR & error)
    {
      if (facet_iter == vnear_facet.end()) {
        error.AddMessage
          ("Programming error. Facet does not appear in map vnear_facet[].");
        error.AddMessage
          ("  Ridge orthogonal directions: (", orth_dir0, ",", orth_dir1, ")");
        error.AddMessage("  Lower/leftmost facet vertex: ", facet_vertex);
        throw error;
      }
    }
    
  };

  
  /*!
   *  @brief Separate isosurface vertices that are near shared grid cube ridges.
   *  - Separate isosurface vertices w0, w1 that are from two different
   *    grid cubes that intersect at a ridge.
   *  - Move both vertices away from the ridge.
   *  @param min_distance Move vertices if both cubes have
   *    isosurface vertices within distance (Linf) \a min_distance
   *    of ridge r.
   *  @param offset Move vertices so that they are distance
   *    \a offset from ridge r.
   *  @param[out] num_repositions Number of times isosurface vertices
   *    are repositioned away from some grid cube ridge.
   *    - Note: An isosurface vertex that is repositioned away
   *      from multiple grid cube ridges will contribute
   *      multiple times to this count.
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename MIN_DIST_TYPE, typename OFFSET_TYPE,
            typename CTYPE, typename NTYPEV, typename NTYPEM>
  void separate_isov_near_shared_grid_cube_ridges
  (const GRID_TYPE & grid,
   const DUAL_ISOV_TYPE isov_list[],
   const NTYPEV num_isov,
   const MIN_DIST_TYPE min_distance,
   const OFFSET_TYPE offset,
   CTYPE * coord,
   NTYPEM & num_moved)
  {
    constexpr int BITSET_SIZE = 8;
    
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE CUBE_INDEX_TYPE;
    typedef FACE_ISOV_PAIR<BITSET_SIZE,VERTEX_INDEX_TYPE,NTYPEV,CTYPE>
      FACE_ISOV_PAIR_TYPE;
    typedef typename FACE_ISOV_PAIR_TYPE::FACET_POS_BITSET_TYPE
      FACET_POS_BITSET_TYPE;
    typedef typename std::vector<FACE_ISOV_PAIR_TYPE>::size_type SIZE_TYPE;

    // Bitset type indicating vertices near ridge.
    typedef std::bitset<BITSET_SIZE> VNEAR_BITSET_TYPE;

    typedef std::unordered_map<VERTEX_INDEX_TYPE, VNEAR_BITSET_TYPE>
      VNEAR_RIDGE_MAP;

    const int dimension = grid.Dimension();
    VNEAR_RIDGE_MAP vnear_ridge;
    std::vector<FACE_ISOV_PAIR_TYPE> ridge_isov_list;
    IJK::ARRAY<CTYPE> cube_coord(dimension);
    IJK::ARRAY<CTYPE> gridv_coord(dimension);
    IJK::PROCEDURE_ERROR error
      ("separate_isov_near_shared_grid_cube_ridges");


    // Array used to track which isosurface vertices have moved.
    // Really only needed to get an accurate count of num_moved, but...
    std::vector<bool> flag_moved(num_isov, false);
    
    // Initialize
    num_moved = 0;
    
    // For every combination of orthogonal directions.
    for (int orth_dir0 = 0; orth_dir0 < dimension; orth_dir0++ ) {
      for (int orth_dir1 = orth_dir0+1; orth_dir1 < dimension;
           orth_dir1++) {

        vnear_ridge.clear();
        ridge_isov_list.clear();

        for (NTYPEV isov = 0; isov < num_isov; isov++) {
          const CTYPE * isov_coord = coord+dimension*isov;
          const CUBE_INDEX_TYPE icube =
            get_isovert_grid_cube_index(isov_list[isov]);
          grid.ComputeCoord(icube, cube_coord.Ptr());


          for (int iside0 = 0; iside0 < 2; iside0++) {
            for (int iside1 = 0; iside1 < 2; iside1++) {

              if (is_isov_close_to_ridge
                  (isov_coord, cube_coord.PtrConst(),
                   orth_dir0, orth_dir1, iside0, iside1,
                   min_distance)) {
                    
                const CTYPE cdiff =
                  distance_to_ridge
                  (isov_coord, cube_coord.PtrConst(),
                   orth_dir0, orth_dir1, iside0, iside1);
                

                const VERTEX_INDEX iv_primary =
                  get_ridge_primary_vertex
                  (grid, icube, orth_dir0, orth_dir1, iside0, iside1);
                
                add_ridge_isov_pair
                  (iv_primary, isov, cdiff, iside0, iside1,
                   ridge_isov_list, vnear_ridge);
              }
            }
          }

        }

        const VNEAR_BITSET_TYPE BITSET_1111 = 15;
        
        // Reset bitset to 1111 whenever there are two vertices
        //   on both sides of the one of the two planes
        //   through the ridge.
        for (auto & pair : vnear_ridge) {
          const VNEAR_BITSET_TYPE vnear_bitset = pair.second;

          if (are_both_bits1(vnear_bitset, 0, 1))
            { (pair.second) = BITSET_1111; }
          else if (are_both_bits1(vnear_bitset, 2, 3))
            { (pair.second) = BITSET_1111; }
        }

        for (SIZE_TYPE i = 0; i < ridge_isov_list.size(); i++) {
          const VERTEX_INDEX_TYPE ridge_vertex =
            ridge_isov_list[i].face_vertex;
          const NTYPEV isov = ridge_isov_list[i].isov;
          const FACET_POS_BITSET_TYPE facet_pos =
            ridge_isov_list[i].facet_pos;
        
          CTYPE * isov_coord = coord+dimension*isov;
          grid.ComputeCoord(ridge_vertex, gridv_coord.Ptr());
        
          auto ridge_iter = vnear_ridge.find(ridge_vertex);
          throw_error_if_ridge_not_found
            (ridge_iter, vnear_ridge, orth_dir0, orth_dir1,
             ridge_vertex, error);

          const VNEAR_BITSET_TYPE vnear_bitset = ridge_iter->second;
        
          // Move isov away from ridge.
          if (vnear_bitset == BITSET_1111) {

            const int iside0 = facet_pos[0];
            const int iside1 = facet_pos[1];
            move_coord_away_from_facet
              (gridv_coord.PtrConst(), orth_dir0, iside0, offset, isov_coord);
            move_coord_away_from_facet
              (gridv_coord.PtrConst(), orth_dir1, iside1, offset, isov_coord);
            
            record_moved(isov, flag_moved, num_moved);
          }
        }
      }
    }
  }


  /*!
   *  @overload
   *  @brief Separate isosurface vertices that are near shared grid cube ridges.
   *    (C++ vector.)
   *  - Version using C++ STL vector for isov_list[] and coord[].
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename MIN_DIST_TYPE, typename OFFSET_TYPE,
            typename CTYPE, typename NTYPEM>
  void separate_isov_near_shared_grid_cube_ridges
  (const GRID_TYPE & grid,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const MIN_DIST_TYPE min_distance,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord,
   NTYPEM & num_moved)
  {
    separate_isov_near_shared_grid_cube_ridges
    (grid, IJK::vector2pointer(isov_list), isov_list.size(),
     min_distance, offset, IJK::vector2pointerNC(coord), num_moved);
  }


  // ******************************************************************
  //! @name Separate isosurface vertices that are near grid vertices.
  // ******************************************************************
  
  // local
  namespace {

    template <typename CTYPE, typename CUBE_CTYPE,
              typename GRIDV_POS_BITSET_TYPE, typename MIN_DISTANCE_TYPE>
    bool is_isov_close_to_grid_vertex
    (const int dimension,
     const CTYPE isov_coord[],
     const CUBE_CTYPE cube_coord[],
     const GRIDV_POS_BITSET_TYPE gridv_pos_bitset,
     const MIN_DISTANCE_TYPE min_distance)
    {
      for (int d = 0; d < dimension; d++) {
        if (gridv_pos_bitset[d] == 0) {
          // Grid vertex is on lower/leftmost side of cube.
          if (isov_coord[d] >= cube_coord[d] + min_distance)
            { return false; }
        }
        else {
          if (isov_coord[d] + min_distance <= cube_coord[d]+1)
            { return false; }
        }
      }

      return true;
    }


    // L-infinity distance from isov to grid vertex.
    template <typename CTYPE, typename CUBE_CTYPE,
              typename GRIDV_POS_BITSET_TYPE>
    CTYPE distance_to_grid_vertex
    (const int dimension,
     const CTYPE isov_coord[],
     const CUBE_CTYPE cube_coord[],
     const GRIDV_POS_BITSET_TYPE gridv_pos_bitset)
    {
      CTYPE cdiff =
        distance_to_facet
        (isov_coord, cube_coord, 0, gridv_pos_bitset[0]);
      
      for (int d = 1; d < dimension; d++) {
        if (gridv_pos_bitset[d] == 0) {
          const CTYPE cdiff_d = isov_coord[d] - cube_coord[d];
          cdiff = std::min(cdiff, cdiff_d);
        }
        else {
          const CTYPE cdiff_d = cube_coord[d] + 1 - isov_coord[d];
          cdiff = std::min(cdiff, cdiff_d);
        }
      }

      return cdiff;
    }
      
    
    // Add (grid vertex)-isov pair.
    template <typename VTYPE, typename ISOV_TYPE,
              typename GRIDV_POS_BITSET_TYPE,
              typename FACE_ISOV_PAIR_TYPE,
              typename CTYPE,
              typename VNEAR_BITSET_TYPE>
    void add_gridv_isov_pair
    (const int dimension, const VTYPE iv, const ISOV_TYPE isov,
     const CTYPE cdiff, 
     const GRIDV_POS_BITSET_TYPE & gridv_pos_bitset,
     std::vector<FACE_ISOV_PAIR_TYPE> & gridv_isov_list,
     std::unordered_map<VTYPE, VNEAR_BITSET_TYPE> & vnear_gridv)
    {
      const FACE_ISOV_PAIR_TYPE
        gridv_isov_pair(iv, isov, gridv_pos_bitset, cdiff);

      gridv_isov_list.push_back(gridv_isov_pair);

      for (int d = 0; d < dimension; d++) {
        IJK::insert_or_update_hash_table_bitset
          (iv, 2*d+gridv_pos_bitset[d], 1, vnear_gridv);
      }
    }


    template <typename BITSET_TYPE>
    void set_bits01_to_1(BITSET_TYPE & bitset)
    {
      bitset[0] = 1;
      bitset[1] = 1;
    }

    
    template <typename CTYPEV, typename VPOS_BITSET_TYPE,
              typename OFFSET_TYPE, typename CTYPE>      
    void move_coord_away_from_gridv
    (const int dimension, const CTYPEV vcoord[],
     const VPOS_BITSET_TYPE vpos_bitset,
     const OFFSET_TYPE offset, CTYPE * isov_coord)
    {
      for (int d = 0; d < dimension; d++) {
        move_coord_away_from_facet
          (vcoord, d, vpos_bitset[d], offset, isov_coord);
      }
    }

    
    // Throw error if grid vertex not found.
    template <typename ITER_TYPE, typename MAP_TYPE, typename VTYPE>
    void throw_error_if_gridv_not_found
    (const ITER_TYPE gridv_iter, const MAP_TYPE & vnear_gridv,
     const VTYPE iv, IJK::ERROR & error)
    {
      if (gridv_iter == vnear_gridv.end()) {
        error.AddMessage
          ("Programming error. Facet does not appear in map vnear_gridv[].");
        error.AddMessage("  Grid vertex: ", iv);
        throw error;
      }
    }
    
  };

  
  /*!
   *  @brief Separate isosurface vertices that are near shared grid vertices.
   *  - Separate isosurface vertices w0, w1 that are from two different
   *    grid cubes that intersect at a grid vertex.
   *  - Move both vertices away from the vertex.
   *  @param min_distance Move vertices if both cubes have
   *    isosurface vertices within distance (Linf) \a min_distance
   *    of grid vertex v.
   *  @param offset Move vertices so that they are distance
   *    \a offset from grid vertex v.
   *  @param[out] num_repositions Number of times isosurface vertices
   *    are repositioned away from some grid vertex.
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename MIN_DIST_TYPE, typename OFFSET_TYPE,
            typename CTYPE, typename NTYPEV, typename NTYPEM>
  void separate_isov_near_shared_grid_vertices
  (const GRID_TYPE & grid,
   const DUAL_ISOV_TYPE isov_list[],
   const NTYPEV num_isov,
   const MIN_DIST_TYPE min_distance,
   const OFFSET_TYPE offset,
   CTYPE * coord,
   NTYPEM & num_moved)
  {
    constexpr int BITSET_SIZE = 8;
    
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE CUBE_INDEX_TYPE;
    typedef typename GRID_TYPE::CUBEV_BITSET_TYPE CUBEV_BITSET_TYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef FACE_ISOV_PAIR_BASE
      <CUBEV_BITSET_TYPE,VERTEX_INDEX_TYPE,NTYPEV,CTYPE> GRIDV_ISOV_PAIR_TYPE;
    typedef typename GRIDV_ISOV_PAIR_TYPE::FACET_POS_BITSET_TYPE
      GRIDV_POS_BITSET_TYPE;
    typedef typename std::vector<GRIDV_ISOV_PAIR_TYPE>::size_type SIZE_TYPE;

    // Bitset type indicating isosurface vertices near grid vertex.
    typedef std::bitset<BITSET_SIZE> VNEAR_BITSET_TYPE;

    typedef std::unordered_map<VERTEX_INDEX_TYPE, VNEAR_BITSET_TYPE>
      VNEAR_GRIDV_MAP;

    const int dimension = grid.Dimension();
    VNEAR_GRIDV_MAP vnear_gridv;
    std::vector<GRIDV_ISOV_PAIR_TYPE> gridv_isov_list;
    IJK::ARRAY<CTYPE> cube_coord(dimension);
    IJK::ARRAY<CTYPE> gridv_coord(dimension);
    IJK::PROCEDURE_ERROR error
      ("separate_isov_near_shared_grid_vertices");

    if (2*dimension > BITSET_SIZE) {
      error.AddMessage("Programming error. Bitset size too small.");
      error.AddMessage
        ("  Bitset size for bitset VPOS_BITSET_TYPE::flag_below: ",
         BITSET_SIZE);
      error.AddMessage
        ("  Bitset size required to be at least ", 2*dimension,
         " for grid with dimension ", dimension, ".");
      throw error;
    }

    // Array used to track which isosurface vertices have moved.
    // Really only needed to get an accurate count of num_moved, but...
    std::vector<bool> flag_moved(num_isov, false);
    
    // Initialize
    num_moved = 0;

    for (NTYPEV isov = 0; isov < num_isov; isov++) {
      const CTYPE * isov_coord = coord+dimension*isov;
      const CUBE_INDEX_TYPE icube =
        get_isovert_grid_cube_index(isov_list[isov]);
      grid.ComputeCoord(icube, cube_coord.Ptr());

      for (NUMBER_TYPE i = 0; i < grid.NumCubeVertices(); i++) {
        const VERTEX_INDEX_TYPE iv = grid.CubeVertex(icube, i);
        const CUBEV_BITSET_TYPE gridv_pos_bitset = grid.CubeVertexBitset(i);
        
        if (is_isov_close_to_grid_vertex
            (dimension, isov_coord, cube_coord.Ptr(), gridv_pos_bitset,
             min_distance)) {

          const CTYPE cdiff =
            distance_to_grid_vertex
            (dimension, isov_coord, cube_coord.Ptr(), gridv_pos_bitset);

          add_gridv_isov_pair
            (dimension, iv, isov, cdiff, gridv_pos_bitset,
             gridv_isov_list, vnear_gridv);
        }
      }
    }

    // Bitset used to identify grid vertices with two isosurface vertices
    //   on opposite side of some plane containing the grid vertex.
    constexpr VNEAR_BITSET_TYPE BITSET_11(3);
        
    // Reset vnear_bitset to BITSET_11 whenever
    //   there are two isosurface vertices on both sides
    //   of some plane containing the grid vertex.
    for (auto & pair : vnear_gridv) {
      const VNEAR_BITSET_TYPE vnear_bitset = pair.second;

      for (int d = 0; d < dimension; d++) {
        if (are_both_bits1(vnear_bitset, 2*d, 2*d+1)) {
          (pair.second) = BITSET_11;
          break;
        }
      }
    }

    for (SIZE_TYPE i = 0; i < gridv_isov_list.size(); i++) {
      const VERTEX_INDEX_TYPE iv =
            gridv_isov_list[i].face_vertex;
      const NTYPEV isov = gridv_isov_list[i].isov;
      const GRIDV_POS_BITSET_TYPE facet_pos =
        gridv_isov_list[i].facet_pos;
        
      CTYPE * isov_coord = coord+dimension*isov;
      grid.ComputeCoord(iv, gridv_coord.Ptr());
        
      auto gridv_iter = vnear_gridv.find(iv);
      throw_error_if_gridv_not_found
        (gridv_iter, vnear_gridv, iv, error);

      const VNEAR_BITSET_TYPE vnear_bitset = gridv_iter->second;

      // Move isov away from ridge.
      if (vnear_bitset == BITSET_11) {
        move_coord_away_from_gridv
          (dimension, gridv_coord.PtrConst(), facet_pos, offset,
           isov_coord);
        
        record_moved(isov, flag_moved, num_moved);
      }
    }

  }


  /*!
   *  @overload
   *  @brief Separate isosurface vertices that are near shared grid vertices.
   *    (C++ vector.)
   *  - Version using C++ STL vector for isov_list[] and coord[].
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename MIN_DIST_TYPE, typename OFFSET_TYPE,
            typename CTYPE, typename NTYPEM>
  void separate_isov_near_shared_grid_vertices
  (const GRID_TYPE & grid,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const MIN_DIST_TYPE min_distance,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord,
   NTYPEM & num_moved)
  {
    separate_isov_near_shared_grid_vertices
    (grid, IJK::vector2pointer(isov_list), isov_list.size(),
     min_distance, offset, IJK::vector2pointerNC(coord), num_moved);
  }

  
  // ******************************************************************
  //! @name Separate isosurface polytopes that are near grid vertices.
  // ******************************************************************

  /*!
   *  @brief Store information about isosurface polytopes that are near grid vertices.
   */
  template <typename ISOP_TYPE, typename GRIDV_TYPE>
  class ISOPOLY_NEAR_GRIDV {
  public:

    /// @brief Index of isosurface polytope.
    ISOP_TYPE isopoly_index;

    /// @brief Grid vertex.
    GRIDV_TYPE gridv_index;

    
    /// @brief True if vertex is below/left of isosurface polytope.
    bool flag_below;

    void Set
    (const ISOP_TYPE isopoly_index, const GRIDV_TYPE gridv_index,
     const bool flag_below)
    {
      this->isopoly_index = isopoly_index;
      this->gridv_index = gridv_index;
      this->flag_below = flag_below;
    }

    /// Constructor.
    ISOPOLY_NEAR_GRIDV
    (const ISOP_TYPE isopoly_index, const GRIDV_TYPE gridv_index,
     const bool flag_below)
    {
      Set(isopoly_index, gridv_index, flag_below);
    }

  };

  
  /*!
   *  @brief Return true if any isopoly vertex is near the plane
   *     orthogonal to orth_dir and with coordinate c.
   */
  template <typename GRID_TYPE, typename ISOV_TYPE,
            typename ISOV_CTYPE, typename CTYPEP, typename DIST_TYPE>
  bool are_isopoly_vert_near_plane
  (const GRID_TYPE & grid,
   const ISOV_TYPE isopoly_vert[],
   const int numv_per_isopoly,
   const ISOV_CTYPE coord[],
   const int orth_dir, 
   const CTYPEP c,
   const bool flag_below,
   const DIST_TYPE min_distance)
  {
    const int dimension = grid.Dimension();
    
    for (int i = 0; i < numv_per_isopoly; i++) {
      const ISOV_TYPE isov = isopoly_vert[i];
      const ISOV_CTYPE isov_coord_d = coord[dimension*isov+orth_dir];

      ISOV_CTYPE cdiff;
      if (flag_below) { cdiff = isov_coord_d - c; }
      else { cdiff = c - isov_coord_d; }

      if (cdiff < min_distance)
        { return true; }
    }

    return false;
  }


  /*!
   *  @brief Separate isosurface polytopes that are near grid vertices
   *     that are shared by their dual edges.
   *  @param isopoly_vert[i*numv_per_isopoly+k]
   *    Index of k'th vertex of isosurface polytope i.
   *  @param isopoly_info[jpoly] Information on isososurface polygon jpoly,
   *    including dual grid edge.
   */
  template <typename GRID_TYPE, typename ISOV_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE,
            typename DUAL_ISOV_TYPE,
            typename MIN_DIST_TYPE, typename OFFSET_TYPE,
            typename CTYPE, typename NTYPEV, typename NTYPEM>
  void separate_iso_poly_near_shared_grid_vertices
  (const GRID_TYPE & grid,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   const DUAL_ISOV_TYPE isov_list[],
   const NTYPEV num_isov,
   const MIN_DIST_TYPE min_distance,
   const OFFSET_TYPE offset,
   CTYPE coord[],
   NTYPEM & num_moved)
  {
    typedef typename GRID_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    // Bitset type indicating isopoly near grid vertex.
    typedef std::bitset<8> PNEAR_BITSET_TYPE;
    typedef ISOPOLY_NEAR_GRIDV<NTYPEV,VTYPE> ISOPOLY_NEAR_GRIDV_TYPE;
    typedef typename std::vector<ISOPOLY_NEAR_GRIDV_TYPE>::size_type SIZE_TYPE;
    
    typedef std::unordered_map<VTYPE, PNEAR_BITSET_TYPE>
      VNEAR_ISOPOLY_MAP;

    
    const int dimension = grid.Dimension();
    const int numv_per_isopoly = grid.NumCubeFacetVertices();
    const NUMBER_TYPE num_isopoly = isopoly_info.size();
    VNEAR_ISOPOLY_MAP vnear_isopoly;
    std::vector<ISOPOLY_NEAR_GRIDV_TYPE> isopoly_near_gridv_list;
    IJK::PROCEDURE_ERROR error
      ("separate_isopoly_near_shared_grid_vertices");

    // Array used to track which isosurface vertices have moved.
    // Really only needed to get an accurate count of num_moved, but...
    std::vector<bool> flag_moved(num_isov, false);

    // Initialize
    num_moved = 0;

    for (int d = 0; d < dimension; d++ ) {

      vnear_isopoly.clear();
      isopoly_near_gridv_list.clear();

      for (NUMBER_TYPE ipoly = 0; ipoly < num_isopoly; ipoly++) {

        if (isopoly_info[ipoly].GridEdgeDirection() != d) { continue; }
        
        const VTYPE iend0 = isopoly_info[ipoly].GridEdgeEndpoint0();
        const VTYPE iend1 = grid.NextVertex(iend0, d);
        const VTYPE iend[2] = { iend0, iend1 };

        const ISOV_INDEX_TYPE * isopoly_vert_i =
          (isopoly_vert.data() + ipoly*numv_per_isopoly);

        bool flag_below = true;
        int iloc = 0;
        for (int j = 0; j < 2; j++) {

          const COORD_TYPE gridv_coord = grid.CoordD(iend[j], d);
    
          if (are_isopoly_vert_near_plane
              (grid, isopoly_vert_i, numv_per_isopoly,
               coord, d, gridv_coord, flag_below, min_distance)) {

            const ISOPOLY_NEAR_GRIDV_TYPE
              isopoly_near_gridv(ipoly, iend[j], flag_below);
            isopoly_near_gridv_list.push_back(isopoly_near_gridv);

            IJK::insert_or_update_hash_table_bitset
              (iend[j], iloc, 1, vnear_isopoly);
          }

          flag_below = !flag_below;
          iloc++;
        }
      }

        
      for (SIZE_TYPE i = 0;
           i < isopoly_near_gridv_list.size(); i++) {

        const NUMBER_TYPE ipoly =
          isopoly_near_gridv_list[i].isopoly_index;
        const VERTEX_INDEX gridv_index =
          isopoly_near_gridv_list[i].gridv_index;
        const bool flag_below = isopoly_near_gridv_list[i].flag_below;

        auto vnear_iter = vnear_isopoly.find(gridv_index);

        if (vnear_iter == vnear_isopoly.end()) {
          error.AddMessage
            ("Programming error. Grid vertex does not apper in map vnear_isopoly[].");
          error.AddMessage("  Isopoly dual grid edge direction: ", d);
          error.AddMessage("  Grid vertex: ", gridv_index);
          throw error;
        }

        const PNEAR_BITSET_TYPE bitset = vnear_iter->second;

        if (flag_below) {
          if (bitset[1] == 1) {
            const ISOV_INDEX_TYPE * isopoly_vert_i =
              (isopoly_vert.data() + ipoly*numv_per_isopoly);
            const COORD_TYPE gridv_coord_d = grid.CoordD(gridv_index, d);
          
            for (int j = 0; j < numv_per_isopoly; j++) {
              const ISOV_INDEX_TYPE isov = isopoly_vert_i[j];
              const COORD_TYPE isov_coord_d = coord[dimension*isov+d];
              const COORD_TYPE cdiff = isov_coord_d - gridv_coord_d;
              if (cdiff < min_distance) {
                coord[dimension*isov+d] = gridv_coord_d + offset;

                if (!flag_moved[isov]) {
                  flag_moved[isov] = true;
                  num_moved++;
                }
              }
            }
          }
        }
        else {
          if (bitset[0] == 1) {
            const ISOV_INDEX_TYPE * isopoly_vert_i =
              (isopoly_vert.data() + ipoly*numv_per_isopoly);
            const COORD_TYPE gridv_coord_d = grid.CoordD(gridv_index, d);
          
            for (int j = 0; j < numv_per_isopoly; j++) {
              const ISOV_INDEX_TYPE isov = isopoly_vert_i[j];
              const COORD_TYPE isov_coord_d = coord[dimension*isov+d];
              const COORD_TYPE cdiff = gridv_coord_d - isov_coord_d;
              if (cdiff < min_distance) {
                coord[dimension*isov+d] = gridv_coord_d - offset;

                if (!flag_moved[isov]) {
                  flag_moved[isov] = true;
                  num_moved++;
                }
              }
            }
          }
        }
      }
    }

  }


  /*!
   *  @overload
   *  @brief Separate isosurface polytopes that are near grid vertices
   *     that are shared by their dual edges. (C++ vector.)
   *  - Version using C++ STL vector for isov_list[] and coord[].
   */
  template <typename GRID_TYPE, typename ISOV_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE, typename DUAL_ISOV_TYPE,
            typename MIN_DIST_TYPE, typename OFFSET_TYPE,
            typename CTYPE, typename NTYPEM>
  void separate_iso_poly_near_shared_grid_vertices
  (const GRID_TYPE & grid,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const MIN_DIST_TYPE min_distance,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord,
   NTYPEM & num_moved)
  {
    separate_iso_poly_near_shared_grid_vertices    
      (grid, isopoly_vert, isopoly_info, IJK::vector2pointer(isov_list),
       isov_list.size(), min_distance, offset,
       IJK::vector2pointerNC(coord), num_moved);
  }


  // ******************************************************************
  //! @name Reposition "doubly connected" isosurface vertices.
  // ******************************************************************

  //@{

  /*!
   *  @brief Return c0 clamped to range
   *    [min(c1A,c1B)+offset, max(c1A,c1B)-offset].
   *  - If |c1A-c1B| <= 2*offset, return (c1A+c1B)/2.0.
   */
  template <typename CTYPE0, typename CTYPE1, typename OFFSET_TYPE>
  CTYPE0 clamp_between_offset
  (const CTYPE0 c0, const CTYPE1 c1A, const CTYPE1 c1B,
   const OFFSET_TYPE offset)
  {
    const CTYPE1 c1min = std::min(c1A,c1B);
    const CTYPE1 c1max = std::max(c1A,c1B);
    CTYPE0 result;

    if ((c1max-c1min) > 2*offset) {
      result =
        IJK::clamp_coord_to_range(c0, c1min+offset, c1max-offset);
    }
    else {
      result = (c1min + c1max)/2.0;
    }
    
    return result;
  }


  /*!
   *  @brief Move coordinates d1 and d2 of isosurface vertex isov0
   *    between coordinates d1 and d2 of isov1 and isov2.
   *  - Apply offset.
   *  - NOTE: Only for dimension 3.
   *  @pre Isosurface vertices isov1 and isov2 are appropriately
   *    separated by at least 2*offset in directions d1 and d2.
   */
  template <typename ISOV_INDEX_TYPE, typename OFFSET_TYPE,
            typename CTYPE>
  void move_isov_coord_between_3D
  (const ISOV_INDEX_TYPE isov0,
   const ISOV_INDEX_TYPE isov1,
   const ISOV_INDEX_TYPE isov2,
   const int d1,
   const int d2,
   const OFFSET_TYPE offset,
   CTYPE * coord)
  {
    CTYPE * isov0_coord = coord + isov0*DIM3;
    const CTYPE * isov1_coord = coord + isov1*DIM3;
    const CTYPE * isov2_coord = coord + isov2*DIM3;

    // Move isov0_coord[d1] to between
    //   isov1_coord[d1] and isov2_coord[d1] (+/- offset).
    // Move isov0_coord[d2] to between
    //   isov1_coord[d2] and isov2_coord[d2] (+/- offset).
    isov0_coord[d1] =
      clamp_between_offset
      (isov0_coord[d1], isov1_coord[d1], isov2_coord[d1], offset);
    isov0_coord[d2] =
      clamp_between_offset
      (isov0_coord[d2], isov1_coord[d2], isov2_coord[d2], offset);
  }

  /*!
   *  Reposition all doubly connected isosurface vertices.
   *  - Clamp all doubly connected vertex coordinates to coordinates
   *    of two adjacent vertices.
   *  - An isosurface vertex is doubly connected to cube c if it is adjacent
   *    to two (or more) isosurface vertices in cube c.
   *  - NOTE: Only for dimension 3.
   *  @pre Isosurface vertices in adjacent cube are already
   *    appropriately separated by at least 2*offset in each direction.
   */
  template <typename GRID_TYPE,
            typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_CUBE_DATA_TYPE, typename NTYPEC,
            typename OFFSET_TYPE, typename CTYPE, typename NTYPEM>
  void reposition_all_doubly_connected_3D_clamp_all
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const GRID_CUBE_DATA_TYPE active_cube_list[],
   const NTYPEC num_active_cubes,
   const OFFSET_TYPE offset,
   CTYPE * coord, NTYPEM & num_moved)
  {
    typedef typename GRID_CUBE_DATA_TYPE::ISOV_TYPE ISOV_TYPE;
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;

    const int DIM3(3);
    std::vector<int> isotable_vlist;
    IJK::PROCEDURE_ERROR
      error("reposition_all_doubly_connected_3D_clamp_all");

    // Initialize
    num_moved = 0;
    
    if (grid.Dimension() != DIM3) {
      error.AddMessage("Programming error. Dimension not equal to 3.");
      error.AddMessage("  grid.Dimension(): ", grid.Dimension());
      throw error;
    }
    
    for (NTYPEC i0 = 0; i0 < num_active_cubes; i0++) {
      const TABLE_INDEX table_index0 =
        active_cube_list[i0].table_index;
      
      if (isodual_table.NumAmbiguousFacets(table_index0) != 1)
        { continue; }

      if (isodual_table.NumIsoVertices(table_index0) != 1)
        { continue; }

      if (!active_cube_list[i0].is_adjacent_cube_set) {
        // Ambiguous facet is a boundary facet.
        continue;
      }

      const ISOV_TYPE isov0 = active_cube_list[i0].first_isov;
      const int ifacet_ambiguous =
        isodual_table.IndexOfSomeAmbiguousFacet(table_index0);
      const int facet_orth_dir =
        IJK::cube_facet_orth_dir(DIM3, ifacet_ambiguous);
      const int iopposite_facet =
        IJK::opposite_cube_facet(DIM3, ifacet_ambiguous);
      const int d1 = (facet_orth_dir+1)%DIM3;
      const int d2 = (facet_orth_dir+2)%DIM3;
      CTYPE * isov0_coord = coord + isov0*DIM3;

      const NTYPEC i1 =
        active_cube_list[i0].list_loc_of_adjacent_cube;

      const TABLE_INDEX table_index1 =
        active_cube_list[i1].table_index;

      const CTYPE c0d1 = isov0_coord[d1];
      const CTYPE c0d2 = isov0_coord[d2];

      isodual_table_vinfo.GetVerticesIncidentOnEdgesCrossingFacet
        (table_index1, iopposite_facet, isotable_vlist);

      if (isotable_vlist.size() != 2) {
        error.AddMessage
          ("Programming error. Incorrect number of vertices returned by function.");
        error.AddMessage
          ("  GetVerticesIncidentOnEdgesCrossingFacet() returned ",
           isotable_vlist.size(), " vertices.");
        error.AddMessage("  instead of expected two vertices.");
        error.AddMessage("  Table index: ", table_index1);
        error.AddMessage("  Facet: ", iopposite_facet);
        throw error;
      }

      const ISOV_TYPE isov1A =
        active_cube_list[i1].first_isov+isotable_vlist[0];
      const ISOV_TYPE isov1B =
        active_cube_list[i1].first_isov+isotable_vlist[1];

      move_isov_coord_between_3D
        (isov0, isov1A, isov1B, d1, d2, offset, coord);

      if ((c0d1 != isov0_coord[d1]) || (c0d2 != isov0_coord[d2])) 
        { num_moved++; }
    }
  }


  /*!
   *  @overload
   *  @brief Reposition all doubly connected isosurface vertices. (C++ vector.)
   *  - Version using C++ vectors for arrays active_cube_list[] and coord[].
   */
  template <typename GRID_TYPE,
            typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_CUBE_DATA_TYPE,
            typename OFFSET_TYPE, typename CTYPE, typename NTYPEM>
  void reposition_all_doubly_connected_3D_clamp_all
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord, NTYPEM & num_moved)
  {
    reposition_all_doubly_connected_3D_clamp_all
      (grid, isodual_table, isodual_table_vinfo,
       IJK::vector2pointer(active_cube_list), active_cube_list.size(),
       offset, IJK::vector2pointerNC(coord), num_moved);
  }

  
  /*!
   *  Reposition all doubly connected isosurface vertices.
   *  - For isosurface vertices doubly connected to cubes containing 3
   *    or more isosurface vertices, align with cube center.
   *  - An isosurface vertex is doubly connected to cube c if it is adjacent
   *    to two (or more) isosurface vertices in cube c.
   *  - NOTE: Only for dimension 3.
   *  @pre Isosurface vertices in adjacent cube are already
   *    appropriately separated by at least 2*offset in each direction.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename GRID_CUBE_DATA_TYPE, typename NTYPEC,
            typename OFFSET_TYPE, typename CTYPE, typename NTYPEM>
  void reposition_all_doubly_connected_3D_centerIII
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const GRID_CUBE_DATA_TYPE active_cube_list[],
   const NTYPEC num_active_cubes,
   const OFFSET_TYPE offset,
   CTYPE * coord, NTYPEM & num_moved)
  {
    typedef typename GRID_CUBE_DATA_TYPE::ISOV_TYPE ISOV_TYPE;
    typedef typename
      GRID_CUBE_DATA_TYPE::CUBE_INDEX_TYPE CUBE_INDEX_TYPE;
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;

    const int DIM3(3);
    CTYPE cube_center_coord[DIM3];
    IJK::PROCEDURE_ERROR
      error("reposition_all_doubly_connected_3D_centerIII");

    // Initialize
    num_moved = 0;
    
    if (grid.Dimension() != DIM3) {
      error.AddMessage("Programming error. Dimension not equal to 3.");
      error.AddMessage("  grid.Dimension(): ", grid.Dimension());
      throw error;
    }
    
    for (NTYPEC i0 = 0; i0 < num_active_cubes; i0++) {
      const TABLE_INDEX table_index0 =
        active_cube_list[i0].table_index;
      
      if (isodual_table.NumAmbiguousFacets(table_index0) != 1)
        { continue; }

      if (isodual_table.NumIsoVertices(table_index0) != 1)
        { continue; }

      if (!active_cube_list[i0].is_adjacent_cube_set) {
        // Ambiguous facet is a boundary facet.
        continue;
      }

      const ISOV_TYPE isov0 = active_cube_list[i0].first_isov;
      const int ifacet_ambiguous =
        isodual_table.IndexOfSomeAmbiguousFacet(table_index0);
      const int facet_orth_dir =
        IJK::cube_facet_orth_dir(DIM3, ifacet_ambiguous);
      const int d1 = (facet_orth_dir+1)%DIM3;
      const int d2 = (facet_orth_dir+2)%DIM3;
      CTYPE * isov0_coord = coord + isov0*DIM3;

      const NTYPEC i1 =
        active_cube_list[i0].list_loc_of_adjacent_cube;

      const TABLE_INDEX table_index1 =
        active_cube_list[i1].table_index;

      const CTYPE c0d1 = isov0_coord[d1];
      const CTYPE c0d2 = isov0_coord[d2];

      if (isodual_table.NumIsoVertices(table_index1) == 2) {
        const ISOV_TYPE isov1A = active_cube_list[i1].first_isov;
        const ISOV_TYPE isov1B = isov1A+1;

        move_isov_coord_between_3D
          (isov0, isov1A, isov1B, d1, d2, offset, coord);
      }
      else {
        // NumIsoVertices(table_index1) > 2.
        const CUBE_INDEX_TYPE icube0 =
          active_cube_list[i0].cube_index;
        grid.ComputeCubeCenterCoord(icube0, cube_center_coord);
        isov0_coord[d1] = cube_center_coord[d1];
        isov0_coord[d2] = cube_center_coord[d2];
      }

      if ((c0d1 != isov0_coord[d1]) || (c0d2 != isov0_coord[d2])) 
        { num_moved++; }
    }
  }


  /*!
   *  @overload
   *  Reposition all doubly connected isosurface vertices. (C++ vector.)
   *  - Version using C++ vectors for arrays active_cube_list[] and coord[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename GRID_CUBE_DATA_TYPE,
            typename OFFSET_TYPE, typename CTYPE, typename NTYPEM>
  void reposition_all_doubly_connected_3D_centerIII
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord, NTYPEM & num_moved)
  {
    reposition_all_doubly_connected_3D_centerIII
      (grid, isodual_table,
       IJK::vector2pointer(active_cube_list), active_cube_list.size(),
       offset, IJK::vector2pointerNC(coord), num_moved);
  }


  /*!
   *  Reposition all doubly connected isosurface vertices.
   *  - Align all vertices with cube center.
   *  - An isosurface vertex is doubly connected to cube c if it is adjacent
   *    to two (or more) isosurface vertices in cube c.
   *  - NOTE: Only for dimension 3.
   *  @pre Isosurface vertices in adjacent cube are already
   *    appropriately separated by at least 2*offset in each direction.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename GRID_CUBE_DATA_TYPE, typename NTYPEC,
            typename OFFSET_TYPE, typename CTYPE, typename NTYPEM>
  void reposition_all_doubly_connected_3D_center_all
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const GRID_CUBE_DATA_TYPE active_cube_list[],
   const NTYPEC num_active_cubes,
   const OFFSET_TYPE offset,
   CTYPE * coord, NTYPEM & num_moved)
  {
    typedef typename GRID_CUBE_DATA_TYPE::ISOV_TYPE ISOV_TYPE;
    typedef typename
      GRID_CUBE_DATA_TYPE::CUBE_INDEX_TYPE CUBE_INDEX_TYPE;
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;

    const int DIM3(3);
    CTYPE cube_center_coord[DIM3];
    IJK::PROCEDURE_ERROR
      error("reposition_all_doubly_connected_3D_center_all");

    // Initialize
    num_moved = 0;
    
    if (grid.Dimension() != DIM3) {
      error.AddMessage("Programming error. Dimension not equal to 3.");
      error.AddMessage("  grid.Dimension(): ", grid.Dimension());
      throw error;
    }
    
    for (NTYPEC i0 = 0; i0 < num_active_cubes; i0++) {
      const TABLE_INDEX table_index0 =
        active_cube_list[i0].table_index;
      
      if (isodual_table.NumAmbiguousFacets(table_index0) != 1)
        { continue; }

      if (isodual_table.NumIsoVertices(table_index0) != 1)
        { continue; }

      if (!active_cube_list[i0].is_adjacent_cube_set) {
        // Ambiguous facet is a boundary facet.
        continue;
      }

      const ISOV_TYPE isov0 = active_cube_list[i0].first_isov;
      const int ifacet_ambiguous =
        isodual_table.IndexOfSomeAmbiguousFacet(table_index0);
      const int facet_orth_dir =
        IJK::cube_facet_orth_dir(DIM3, ifacet_ambiguous);
      const int d1 = (facet_orth_dir+1)%DIM3;
      const int d2 = (facet_orth_dir+2)%DIM3;
      CTYPE * isov0_coord = coord + isov0*DIM3;

      const CTYPE c0d1 = isov0_coord[d1];
      const CTYPE c0d2 = isov0_coord[d2];

      // NumIsoVertices(table_index1) > 2.
      const CUBE_INDEX_TYPE icube0 =
        active_cube_list[i0].cube_index;
      grid.ComputeCubeCenterCoord(icube0, cube_center_coord);
      isov0_coord[d1] = cube_center_coord[d1];
      isov0_coord[d2] = cube_center_coord[d2];

      if ((c0d1 != isov0_coord[d1]) || (c0d2 != isov0_coord[d2])) 
        { num_moved++; }
    }
  }

  
  /*!
   *  @overload
   *  Reposition all doubly connected isosurface vertices. (C++ vector.)
   *  - Version using C++ vectors for arrays active_cube_list[] and coord[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename GRID_CUBE_DATA_TYPE,
            typename OFFSET_TYPE, typename CTYPE, typename NTYPEM>
  void reposition_all_doubly_connected_3D_center_all
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord, NTYPEM & num_moved)
  {
    reposition_all_doubly_connected_3D_center_all
      (grid, isodual_table,
       IJK::vector2pointer(active_cube_list), active_cube_list.size(),
       offset, IJK::vector2pointerNC(coord), num_moved);
  }

  
  /*!
   *  Reposition all doubly connected isosurface vertices.
   *  - For isosurface vertices doubly connected to cubes containing 3
   *    or more isosurface vertices, clamp to [clamp_min,clamp_max] of cube.
   *  - An isosurface vertex is doubly connected to cube c if it is adjacent
   *    to two (or more) isosurface vertices in cube c.
   *  - NOTE: Only for dimension 3.
   *  @pre Isosurface vertices in adjacent cube are already
   *    appropriately separated by at least 2*offset in each direction.
   *  @param clamp_min Clamp lower bound.
   *    @pre In range [0,clamp_max].
   *  @param clamp_max Clamp upper bound.
   *    @pre In range [clamp_min, 1].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename GRID_CUBE_DATA_TYPE, typename NTYPEC,
            typename OFFSET_TYPE, typename CLAMP_TYPE,
            typename CTYPE, typename NTYPEM>
  void reposition_all_doubly_connected_3D_clampIII
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const GRID_CUBE_DATA_TYPE active_cube_list[],
   const NTYPEC num_active_cubes,
   const CLAMP_TYPE & clamp_min, const CLAMP_TYPE & clamp_max,
   const OFFSET_TYPE offset,
   CTYPE * coord, NTYPEM & num_moved)
  {
    typedef typename GRID_CUBE_DATA_TYPE::ISOV_TYPE ISOV_TYPE;
    typedef typename
      GRID_CUBE_DATA_TYPE::CUBE_INDEX_TYPE CUBE_INDEX_TYPE;
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;

    const int DIM3(3);
    CTYPE cube_coord[DIM3];
    IJK::PROCEDURE_ERROR
      error("reposition_all_doubly_connected_3D_clampIII");

    // Initialize
    num_moved = 0;
    
    if (grid.Dimension() != DIM3) {
      error.AddMessage("Programming error. Dimension not equal to 3.");
      error.AddMessage("  grid.Dimension(): ", grid.Dimension());
      throw error;
    }
    
    for (NTYPEC i0 = 0; i0 < num_active_cubes; i0++) {
      const TABLE_INDEX table_index0 =
        active_cube_list[i0].table_index;
      
      if (isodual_table.NumAmbiguousFacets(table_index0) != 1)
        { continue; }

      if (isodual_table.NumIsoVertices(table_index0) != 1)
        { continue; }

      if (!active_cube_list[i0].is_adjacent_cube_set) {
        // Ambiguous facet is a boundary facet.
        continue;
      }

      const ISOV_TYPE isov0 = active_cube_list[i0].first_isov;
      const int ifacet_ambiguous =
        isodual_table.IndexOfSomeAmbiguousFacet(table_index0);
      const int facet_orth_dir =
        IJK::cube_facet_orth_dir(DIM3, ifacet_ambiguous);
      const int d1 = (facet_orth_dir+1)%DIM3;
      const int d2 = (facet_orth_dir+2)%DIM3;
      CTYPE * isov0_coord = coord + isov0*DIM3;

      const NTYPEC i1 =
        active_cube_list[i0].list_loc_of_adjacent_cube;

      const TABLE_INDEX table_index1 =
        active_cube_list[i1].table_index;

      const CTYPE c0d1 = isov0_coord[d1];
      const CTYPE c0d2 = isov0_coord[d2];

      if (isodual_table.NumIsoVertices(table_index1) == 2) {
        const ISOV_TYPE isov1A = active_cube_list[i1].first_isov;
        const ISOV_TYPE isov1B = isov1A+1;

        move_isov_coord_between_3D
          (isov0, isov1A, isov1B, d1, d2, offset, coord);
      }
      else {
        // NumIsoVertices(table_index1) > 2.
        // All isosurface vertices in icube1 have degree 3
        // and should be within distance 1/3 of cube facets incident
        // on the cube corner that is cut off by the isosurface patch.
        const CUBE_INDEX_TYPE icube0 =
          active_cube_list[i0].cube_index;
        grid.ComputeCoord(icube0, cube_coord);
        isov0_coord[d1] =
          clamp_between_offset
          (c0d1, cube_coord[d1] + clamp_min,
           cube_coord[d1] + clamp_max, offset);
        isov0_coord[d2] =
          clamp_between_offset
          (c0d2, cube_coord[d2] + clamp_min,
           cube_coord[d2] + clamp_max, offset);
      }

      if ((c0d1 != isov0_coord[d1]) || (c0d2 != isov0_coord[d2])) 
        { num_moved++; }
    }
  }


  /*!
   *  Reposition all doubly connected isosurface vertices.
   *  - For isosurface vertices doubly connected to cubes containing 3
   *    or more isosurface vertices, clamp to [1/3,2/3] of cube.
   *  - An isosurface vertex is doubly connected to cube c if it is adjacent
   *    to two (or more) isosurface vertices in cube c.
   *  - NOTE: Only for dimension 3.
   *  @pre Isosurface vertices in adjacent cube are already
   *    appropriately separated by at least 2*offset in each direction.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename GRID_CUBE_DATA_TYPE, typename NTYPEC,
            typename OFFSET_TYPE, typename CTYPE, typename NTYPEM>
  void reposition_all_doubly_connected_3D_clampIII_one_third
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const GRID_CUBE_DATA_TYPE active_cube_list[],
   const NTYPEC num_active_cubes,
   const OFFSET_TYPE offset,
   CTYPE * coord, NTYPEM & num_moved)
  {
    const CTYPE ONE_THIRD(1.0/3.0);
    const CTYPE TWO_THIRDS(2.0/3.0);
    
    reposition_all_doubly_connected_3D_clampIII
      (grid, isodual_table, active_cube_list, num_active_cubes,
       ONE_THIRD, TWO_THIRDS, offset, coord, num_moved);
  }


  /*!
   *  @overload
   *  Reposition all doubly connected isosurface vertices. (C++ vector.)
   *  - Version using C++ vectors for arrays active_cube_list[] and coord[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename GRID_CUBE_DATA_TYPE,
            typename OFFSET_TYPE, typename CTYPE, typename NTYPEM>
  void reposition_all_doubly_connected_3D_clampIII_one_third
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord, NTYPEM & num_moved)
  {
    reposition_all_doubly_connected_3D_clampIII_one_third
      (grid, isodual_table,
       IJK::vector2pointer(active_cube_list), active_cube_list.size(),
       offset, IJK::vector2pointerNC(coord), num_moved);
  }

  
  /*!
   *  Reposition all doubly connected isosurface vertices.
   *  @pre Isosurface vertices in adjacent cube are already
   *    appropriately separated by at least 2*offset in each direction.
   *  @param doubly_connected_isov_reposition_method
   *    Method for repositioning doubly connected isosurface vertices.
   *    @pre Method is NOT DC_POS_CLAMP_ADJACENT since that method
   *      requires DUAL_TABLE_VERTEX_INFO.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename GRID_CUBE_DATA_TYPE,
            typename OFFSET_TYPE, typename CTYPE, typename NTYPEM>
  void reposition_all_doubly_connected_isov_cube3D
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   const OFFSET_TYPE offset,
   const DOUBLY_CONNECTED_ISOV_REPOSITION_METHOD
   doubly_connected_isov_reposition_method,
   std::vector<CTYPE> & coord, NTYPEM & num_moved)
  {
    IJK::PROCEDURE_ERROR
      error("reposition_all_doubly_connected_isov_cube3D");
    
    switch(doubly_connected_isov_reposition_method) {

    case DC_POS_CLAMPIII_ONE_THIRD:
      reposition_all_doubly_connected_3D_clampIII_one_third
        (grid, isodual_table, active_cube_list, offset, coord, num_moved);
      break;

    case DC_POS_CENTERIII:
      reposition_all_doubly_connected_3D_centerIII
        (grid, isodual_table, active_cube_list, offset, coord, num_moved);
      break;

    case DC_POS_CENTER_ALL:
      reposition_all_doubly_connected_3D_center_all
        (grid, isodual_table, active_cube_list, offset, coord, num_moved);
      break;      
      
    case DC_POS_CLAMP_ALL:
      error.AddMessage
        ("Programming error. Method DC_POS_CLAMP_ALL not available from this routine.");
      error.AddMessage
        ("Method DC_POS_CLAMP_ALL requires DUAL_TABLE_VERTEX_INFO.");
      throw error;
      break;
      
    default:
      error.AddMessage
        ("Programming error. Illegal doubly_connected_isov_reposition_method.");
      throw error;
    }

  }


  //@}
  
   
  // ******************************************************************
  //! @name Separate isosurface vertices by cube centers or planes
  // ******************************************************************

  //@{


  /*!
   *  @brief Determine the relative position of two isosurface vertices in a cube.
   *    - Isosurface vertices are identified by index of patch containing vertex.
   *  @param[out] isov0_coord_lt_isov1_coord
   *    If true, then isov0_coord[dir] should be less than isov1_coord[dir].
   *    - isov0_coord_lt_isov1_coord and isov0_coord_gt_isov1_coord cannot
   *      both be true.
   *    - isov0_coord_lt_isov1_coord and isov0_coord_gt_isov1_coord could
   *      both be false, if there is no restriction on their relative positions.
   *  @param[out] isov0_coord_gt_isov1_coord
   *    If true, then isov_coord[dir] should be greater than isov1_coord[dir].
   */
  template <typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename TABLE_INDEX_TYPE, typename IPATCH_INDEX_TYPE>
  void determine_relative_position_of_isov_in_cube
  (const int dimension,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const TABLE_INDEX_TYPE table_index,
   const IPATCH_INDEX_TYPE ipatch0,
   const IPATCH_INDEX_TYPE ipatch1,
   const int dir,
   bool & flag_isov0_coord_lt_isov1_coord,
   bool & flag_isov0_coord_gt_isov1_coord)
  {
    // Initialize
    flag_isov0_coord_lt_isov1_coord = false;
    flag_isov0_coord_gt_isov1_coord = false;

    const bool flag_intersects_lower_facet0 =
      isodual_table_vinfo.VertexInfo(table_index, ipatch0).IntersectsFacet(dir);
    const bool flag_intersects_upper_facet0 =
      isodual_table_vinfo.VertexInfo(table_index, ipatch0).IntersectsFacet(dir+dimension);
    const bool flag_intersects_lower_facet1 =
      isodual_table_vinfo.VertexInfo(table_index, ipatch1).IntersectsFacet(dir);
    const bool flag_intersects_upper_facet1 =
      isodual_table_vinfo.VertexInfo(table_index, ipatch1).IntersectsFacet(dir+dimension);

    if (flag_intersects_lower_facet0 &&
        flag_intersects_upper_facet1) {
      if (!flag_intersects_upper_facet0 ||
          !flag_intersects_lower_facet1) {
        flag_isov0_coord_lt_isov1_coord = true;
      }
    }
    else if (flag_intersects_upper_facet0 &&
             flag_intersects_lower_facet1) {
      if (!flag_intersects_lower_facet0 ||
          !flag_intersects_upper_facet1) {
        flag_isov0_coord_gt_isov1_coord = true;
      }
    }
  }


  /*!
   *  @brief Separate two coordinates by a single axis parallel plane.
     *  @param plane_coord_d Coordinate of plane in direction orthogonal 
   *    to plane.
   *  @param offset Offset from separating plane.
   *  @param[out] coord0_d Coordinate that must be left/below plane.
   *  @param[out] coord1_d Coordinate that must be right/above plane.
   *  @param[out] flag_changed0 True if coord0_d changes.
   *  @param[out] flag_changed1 True if coord1_d changes.
   */
  template <typename COORDP_TYPE, typename OFFSET_TYPE,
            typename COORDC_TYPE>
  inline void separate_coord_by_axis_parallel_plane
  (const COORDP_TYPE plane_coord_d,
   const OFFSET_TYPE offset,
   COORDC_TYPE & coord0_d,
   COORDC_TYPE & coord1_d,
   bool & flag_changed0,
   bool & flag_changed1)
  {
    if (coord0_d > plane_coord_d - offset) {
      coord0_d = plane_coord_d - offset;
      flag_changed0 = true;
    }
    else
      { flag_changed0 = false; }
    
    if (coord1_d < plane_coord_d + offset) {
      coord1_d = plane_coord_d + offset;
      flag_changed1 = true;
    }
    else
      { flag_changed1 = false; }
  }
  
 
  /*!
   *  @brief Separate two isosurface vertices by axis parallel planes.
   *  - Reposition two isosurface vertices in same cube.
   *  - Let midc[d] = (isov0_coord[d]+isov1_coord[d])/2.
   *  - Coordinate plane_coord[d] of separating plane orthogonal to d is:
   *    -# midc[d], if midc[d] is in range cube_coord[d]+[2*offset,1-2*offset],
   *    -# cube_cood[d]+2*offset if midc[d] < cube_coord[d]+2*offset,   
   *    -# cube_coord[d]+1-2*offset if midc[d] > cube_coord[d]-2*offset.
   *  - Reposition based on the isosurface edges incident 
   *    on each isosurface vertex.
   *  - If isov0 has an incident edge dual to facet f,
   *    but has no incident edge dual to the facet f' parallel to f,
   *    and isov1 has an incident edge dual to f or f' or both,
   *    then reposition isov0 closer to f and isov1 closer to f'.
   *  - If isov1 has an incident edge dual to facet f,
   *    but has no incident edge dual to the facet f' parallel to f,
   *    and isov0 has an incident edge dual to f or f' or both,
   *    then reposition isov0 closer to f' and isov1 closer to f.
   *  - Reposition isov0 so that:
   *      -# |isov0_coord[d] - isov1_coord[d]| >= 2*offset, and
   *      -# |isov0_coord[d] - plane_coord[d]| >= offset, and
   *      -# |isov1_coord[d] - plane_coord[d]| >= offset.
   *  @param offset Offset. Must be in range [0,0.25].
   *  @param flag_moved0[out] True if position of isov0 changes.
   *  @param flag_moved1[out] True if position of isov1 changes.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename ISOV_INDEX_TYPE, typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE,
            typename CTYPE>
  void separate_dual_isov_by_axis_parallel_planes
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const ISOV_INDEX_TYPE isov0,
   const ISOV_INDEX_TYPE isov1,
   const DUAL_ISOV_TYPE isov_list[],
   const OFFSET_TYPE offset,
   CTYPE * coord,
   bool & flag_moved0,
   bool & flag_moved1)
  {
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_ISOV_TYPE::CUBE_INDEX_TYPE CUBE_INDEX;
    typedef typename DUAL_ISOV_TYPE::PATCH_INDEX_TYPE PATCH_INDEX;
    
    const int dimension = grid.Dimension();
    const CUBE_INDEX icube0 = isov_list[isov0].cube_index;
    const TABLE_INDEX table_index = isov_list[isov0].table_index;
    const PATCH_INDEX ipatch0 = isov_list[isov0].patch_index;
    const PATCH_INDEX ipatch1 = isov_list[isov1].patch_index;    
    IJK::ARRAY<CTYPE> cube_coord(dimension);

    grid.ComputeCoord(icube0, cube_coord.Ptr());

    CTYPE * isov0_coord = coord + dimension*isov0;
    CTYPE * isov1_coord = coord + dimension*isov1;

    flag_moved0 = false;
    flag_moved1 = false;
    for (int d = 0; d < dimension; d++) {

      // At most one of these two flags should be true.
      bool flag_isov0_coord_lt_isov1_coord = false;
      bool flag_isov0_coord_gt_isov1_coord = false;

      determine_relative_position_of_isov_in_cube
        (dimension, isodual_table, isodual_table_vinfo,
         table_index, ipatch0, ipatch1, d,
         flag_isov0_coord_lt_isov1_coord,
         flag_isov0_coord_gt_isov1_coord);

      CTYPE plane_coord_d =
        (isov0_coord[d] + isov1_coord[d])/2.0;
      const CTYPE cube_coord_d = cube_coord[d];

      // Ensure that plane is not too close to cube facets.
      if (plane_coord_d < cube_coord_d + 2*offset) 
        { plane_coord_d = cube_coord_d + 2*offset; }
      if (plane_coord_d > cube_coord_d + 1 - 2*offset)
        { plane_coord_d = cube_coord_d + 1 - 2*offset; }

      bool flag_moved0_d(false), flag_moved1_d(false);
      if (flag_isov0_coord_lt_isov1_coord) {
        separate_coord_by_axis_parallel_plane
          (plane_coord_d, offset, isov0_coord[d], isov1_coord[d],
           flag_moved0_d, flag_moved1_d);
      }
      else if (flag_isov0_coord_gt_isov1_coord) {
        separate_coord_by_axis_parallel_plane
          (plane_coord_d, offset, isov1_coord[d], isov0_coord[d],
           flag_moved0_d, flag_moved1_d);
      }

      flag_moved0 = (flag_moved0 || flag_moved0_d);
      flag_moved1 = (flag_moved1 || flag_moved1_d);
    }
  }


  /*!
   *  @brief Separate all pairs of isosurface vertices by axis parallel planes.
   *  - Reposition isosurface vertices in cube containing
   *    multiple isosurface vertices.
   *  - Reposition based on the isosurface edges incident 
   *    on each isosurface vertex.
   *  - If isov0 has an incident edge dual to facet f,
   *    but has no incident edge dual to the facet f' parallel to f,
   *    and isov1 has an incident edge dual to f or f' or both,
   *    then reposition isov0 and isov1 so that isov0 is closer to f.
   *  - If isov1 has an incident edge dual to facet f,
   *    but has no incident edge dual to the facet f' parallel to f,
   *    and isov0 has an incident edge dual to f or f' or both,
   *    then reposition isov0 and isov1 so that isov1 is closer to f.
   *  - Reposition isov0 so that:
   *      |isov0_coord[d] - isov1_coord[d]| >= 2*offset.
   *  @param offset Offset. Must be in range [0,0.25].
   *  @param num_moved[out] Number of isosurface vertex pairs
   *    where at least one isosurface vertex is moved.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_CUBE_DATA_TYPE, typename NTYPEC,
            typename ISOV_INDEX_TYPE, typename NTYPEP,
            typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE,
            typename CTYPE, typename NTYPEM>
  void separate_all_dual_isov_by_axis_parallel_planes
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const GRID_CUBE_DATA_TYPE active_cube_list[],
   const NTYPEC num_active_cubes,
   const ISOV_INDEX_TYPE isopoly_vert[],
   const NTYPEP num_isopoly,
   const DUAL_ISOV_TYPE isov_list[],
   const OFFSET_TYPE offset,
   CTYPE * coord,
   NTYPEM & num_moved)
  {
    bool flag_moved0, flag_moved1;

    for (NTYPEC i = 0; i < num_active_cubes; i++) {
      for (int j0 = 0; j0+1 < active_cube_list[i].num_isov; j0++) {
        for (int j1 = j0+1; j1 < active_cube_list[i].num_isov; j1++) {
          const ISOV_INDEX_TYPE isov0 = active_cube_list[i].first_isov+j0;
          const ISOV_INDEX_TYPE isov1 = active_cube_list[i].first_isov+j1;
          separate_dual_isov_by_axis_parallel_planes
            (grid, isodual_table, isodual_table_vinfo,
             isov0, isov1, isov_list, offset, coord,
             flag_moved0, flag_moved1);
          if (flag_moved0) { num_moved++; }
          if (flag_moved1) { num_moved++; }
        }
      }
    }

  }

  
  /*!
   *  @brief Separate all pairs of isosurface vertices by axis parallel planes.
   *  - Version with C++ STL vector format for array isov_list[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_CUBE_DATA_TYPE,
            typename ISOV_INDEX_TYPE,
            typename DUAL_ISOV_TYPE, typename OFFSET_TYPE,
            typename CTYPE, typename NTYPEM>
  void separate_all_dual_isov_by_axis_parallel_planes
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const OFFSET_TYPE offset,
   CTYPE * coord, NTYPEM & num_moved)
  {
    typedef typename std::vector<ISOV_INDEX_TYPE>::size_type SIZE_TYPE;

    const SIZE_TYPE num_vert_per_isopoly = grid.NumCubeFacetVertices();
    const SIZE_TYPE num_isopoly = isopoly_vert.size()/num_vert_per_isopoly;
    const SIZE_TYPE num_active_cubes = active_cube_list.size();

    separate_all_dual_isov_by_axis_parallel_planes
    (grid, isodual_table, isodual_table_vinfo,
     IJK::vector2pointer(active_cube_list), num_active_cubes,
     IJK::vector2pointer(isopoly_vert), num_isopoly,
     IJK::vector2pointer(isov_list),
     offset, coord, num_moved);
  }


  /*!
   *  @brief Separate all pairs of isosurface vertices by axis parallel planes.
   *  - Version with C++ STL vector format for array isov_list[].
   *  - Version with C++ STL vector format for array coord[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_CUBE_DATA_TYPE,
            typename ISOV_INDEX_TYPE,
            typename DUAL_ISOV_TYPE, typename OFFSET_TYPE,
            typename CTYPE, typename NTYPEM>
  void separate_all_dual_isov_by_axis_parallel_planes
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord, NTYPEM & num_moved)
  {
    separate_all_dual_isov_by_axis_parallel_planes
    (grid, isodual_table, isodual_table_vinfo, active_cube_list,
     isopoly_vert, isov_list, offset,
     IJK::vector2pointerNC(coord), num_moved);
  }


  /*!
   *  @brief Separate isosurface vertex by cube center from other isosurface vertices.
   *  - Reposition isosurface vertices in cube containing
   *    multiple isosurface vertices.
   *  - Reposition based on the isosurface edges incident 
   *    on each isosurface vertex.
   *  - If isov0 has an incident edge dual to facet f,
   *    but has no incident edge dual to the facet f' parallel to f,
   *    and isov1 has an incident edge dual to f or f' or both,
   *    then reposition isov0 and isov1 so that isov0 is closer to f.
   *  - If isov1 has an incident edge dual to facet f,
   *    but has no incident edge dual to the facet f' parallel to f,
   *    and isov0 has an incident edge dual to f or f' or both,
   *    then reposition isov0 and isov1 so that isov1 is closer to f.
   *  - Reposition isov0 so that isov_coord[d] <= cube_coord[d]+0.5-offset
   *    or isov_coord[d] >= cube_coord[d]+0.5+offset, as appropriate.
   *  @param active_cube Active cube containing isosurface vertex isov0.
   *    @pre active_cube.first_isov <= isov0 < active_cube.first_isov+active_cube.num_isov.
   *  @param offset Offset. Must be in range [0,0.25].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename ISOV_INDEX_TYPE, typename PATCH_INDEX_TYPE,
            typename GRID_CUBE_DATA_TYPE, typename OFFSET_TYPE,
            typename CTYPE>
  void separate_dual_isov_by_cube_center
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const ISOV_INDEX_TYPE isov0,
   const PATCH_INDEX_TYPE ipatch0,
   const GRID_CUBE_DATA_TYPE & active_cube,
   const OFFSET_TYPE offset,
   CTYPE * coord)
  {
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename GRID_CUBE_DATA_TYPE::CUBE_INDEX_TYPE CUBE_INDEX;

    const int dimension = grid.Dimension();
    const CUBE_INDEX icube0 = active_cube.CubeIndex();
    const TABLE_INDEX it0 = active_cube.table_index;
    const PATCH_INDEX_TYPE num_patches =
      isodual_table.NumIsoVertices(it0);
    IJK::ARRAY<CTYPE> cube_coord(dimension);

    if (num_patches < 2) {
      // No repositioning necessary.
      return;
    }

    grid.ComputeCoord(icube0, cube_coord.Ptr());
    CTYPE * isov0_coord = coord + dimension*isov0;

    for (PATCH_INDEX_TYPE ipatch1 = 0; ipatch1 < num_patches; ipatch1++) {

      if (ipatch0 == ipatch1) { continue; }

      for (int d = 0; d < dimension; d++) {

        // At most one of these two flags should be true.
        bool flag_isov0_coord_lt_isov1_coord = false;
        bool flag_isov0_coord_gt_isov1_coord = false;

        determine_relative_position_of_isov_in_cube
          (dimension, isodual_table, isodual_table_vinfo,
           it0, ipatch0, ipatch1, d,
           flag_isov0_coord_lt_isov1_coord,
           flag_isov0_coord_gt_isov1_coord);

        if (flag_isov0_coord_lt_isov1_coord) {
          const CTYPE cmax = cube_coord[d] + 0.5 - offset;
          if (isov0_coord[d] > cmax)
            { isov0_coord[d] = cmax; }
        }
        else if (flag_isov0_coord_gt_isov1_coord) {
            const CTYPE cmin = cube_coord[d] + 0.5 + offset;
            
            if (isov0_coord[d] < cmin)
              { isov0_coord[d] = cmin; }
        }
      }

    }
  }
  
    
  /*!
   *  @brief Separate all isosurface vertices by cube centers from other isosurface vertices.
   *  - Reposition isosurface vertices in cubes containing multiple isosurface vertices.
   *  - If isosurface vertices isov0 and isov1 are in the same cube
   *    and each has an incident edge dual to the same cube facet f,
   *    then reposition isov0, if necessary.
   *  - Reposition isov0 so that isov_coord[d] <= cube_coord[d]+0.5-offset
   *    or isov_coord[d] >= cube_coord[d]+0.5+offset, as appropriate.
   *  @param offset Offset. Must be in range [0,0.25].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_CUBE_DATA_TYPE, typename NTYPEC,
            typename ISOV_INDEX_TYPE, typename NTYPEP,
            typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE,
            typename CTYPE>
  void separate_all_dual_isov_by_cube_centers
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const GRID_CUBE_DATA_TYPE active_cube_list[],
   const NTYPEC num_active_cubes,
   const ISOV_INDEX_TYPE isopoly_vert[],
   const NTYPEP num_isopoly,
   const DUAL_ISOV_TYPE isov_list[],
   const OFFSET_TYPE offset,
   CTYPE * coord)
  {
    for (NTYPEC i = 0; i < num_active_cubes; i++) {
      ISOV_INDEX_TYPE isov = active_cube_list[i].first_isov;
      for (int j = 0; j < active_cube_list[i].num_isov; j++) {
        separate_dual_isov_by_cube_center
          (grid, isodual_table, isodual_table_vinfo,
           isov, j, active_cube_list[i], offset, coord);
        isov++;
      }
    }
  }


  /*!
   *  @brief Separate all isosurface vertices by cube centers from other isosurface vertices.
   *  @brief Reposition all isosurface vertices in cubes containing multiple isosurface vertices.
   *  - Version with C++ STL vector format for array isov_list[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_CUBE_DATA_TYPE,
            typename ISOV_INDEX_TYPE,
            typename DUAL_ISOV_TYPE, typename OFFSET_TYPE,
            typename CTYPE>
  void separate_all_dual_isov_by_cube_centers
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const OFFSET_TYPE offset,
   CTYPE * coord)
  {
    typedef typename std::vector<ISOV_INDEX_TYPE>::size_type SIZE_TYPE;

    const SIZE_TYPE num_vert_per_isopoly = grid.NumCubeFacetVertices();
    const SIZE_TYPE num_isopoly = isopoly_vert.size()/num_vert_per_isopoly;
    const SIZE_TYPE num_active_cubes = active_cube_list.size();

    separate_all_dual_isov_by_cube_centers
    (grid, isodual_table, isodual_table_vinfo,
     IJK::vector2pointer(active_cube_list), num_active_cubes,
     IJK::vector2pointer(isopoly_vert), num_isopoly,
     IJK::vector2pointer(isov_list),
     offset, coord);
  }


  /*!
   *  @brief Separate all isosurface vertices by cube centers from other isosurface vertices.
   *  - Version with C++ STL vector format for array isov_list[].
   *  - Version with C++ STL vector format for array coord[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_CUBE_DATA_TYPE,
            typename ISOV_INDEX_TYPE,
            typename DUAL_ISOV_TYPE, typename OFFSET_TYPE,
            typename CTYPE>
  void separate_all_dual_isov_by_cube_centers
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord)
  {
    separate_all_dual_isov_by_cube_centers
    (grid, isodual_table, isodual_table_vinfo, active_cube_list,
     isopoly_vert, isov_list, offset, IJK::vector2pointerNC(coord));
  }

  //@}


  // ******************************************************************
  //! @name Separate isosurface edges crossing ambiguous facets
  // ******************************************************************

  //@{

  /*!
   *  @brief Separate isosurface vertex isovB0 
   *    from isosurface edge (isovA0,isovA1).
   *  - Moves isovB0, if not already separated.
   *  - Does not move isovA0 or isovA1.
   *  - Will not move isovB0 closer than offset to cube facet.
   *  @pre Isosurface vertices isovA0 and isovB0 are in the same cube.
   *  @param shared_facet_orth_dir Direction orthogonal to shared facet.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename DUAL_ISOV_TYPE, typename ISOV_INDEX_TYPE,
            typename CUBE_CTYPE, typename DIR_TYPE,
            typename OFFSET_TYPE, typename CTYPE>
  bool separate_isov_from_iso_edge
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const DUAL_ISOV_TYPE isov_list[], 
   const ISOV_INDEX_TYPE isovB0,
   const ISOV_INDEX_TYPE isovA0,
   const ISOV_INDEX_TYPE isovA1,
   const CUBE_CTYPE cube_coord[],
   const DIR_TYPE shared_facet_orth_dir,
   const OFFSET_TYPE offset,
   CTYPE * coord)
  {
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_ISOV_TYPE::PATCH_INDEX_TYPE PATCH_INDEX;
    
    const int dimension = grid.Dimension();
    const TABLE_INDEX table_index = isov_list[isovA0].table_index;
    const PATCH_INDEX ipatchA = isov_list[isovA0].patch_index;
    const PATCH_INDEX ipatchB = isov_list[isovB0].patch_index;
    CTYPE * isovB0_coord = coord + dimension*isovB0;
    CTYPE * isovA0_coord = coord + dimension*isovA0;
    CTYPE * isovA1_coord = coord + dimension*isovA1;
    bool flag_isovA_coord_lt_isovB_coord;
    bool flag_isovA_coord_gt_isovB_coord;

    bool flag_moved = false;

    for (int dir = 0; dir < dimension; dir++) {

      if (dir == shared_facet_orth_dir) {
        // Skip coordinate corresponding to shared_facet_orth_dir.
        continue;
      }

      determine_relative_position_of_isov_in_cube
        (dimension, isodual_table, isodual_table_vinfo,
         table_index, ipatchA, ipatchB, dir,
         flag_isovA_coord_lt_isovB_coord,
         flag_isovA_coord_gt_isovB_coord);

      if (flag_isovA_coord_lt_isovB_coord) {
        const CTYPE cmax =
          std::max(isovA0_coord[dir], isovA1_coord[dir]);
        const CTYPE plane_coord_d =
          std::min(cube_coord[dir]+1-2*offset, cmax);
        
        if (isovB0_coord[dir] < plane_coord_d + offset) {
          isovB0_coord[dir] = plane_coord_d + offset;
          flag_moved = true;
        }
      }
      else if (flag_isovA_coord_gt_isovB_coord) {
        const CTYPE cmin =
          std::min(isovA0_coord[dir], isovA1_coord[dir]);
        const CTYPE plane_coord_d =
          std::max(cube_coord[dir]+2*offset, cmin);
        
        if (isovB0_coord[dir] > plane_coord_d - offset) {
          isovB0_coord[dir] = plane_coord_d - offset;
          flag_moved = true;
        }
      }
    }

    return flag_moved;
  }

  
  /*!
   *  @brief Separate isosurface vertices in cube containing isovA0
   *    from isosurface edge (isovA0,isovA1).
   *    - Separate only vertices incident on edges crossing ishared_facet.
   *  @param ishared_facet Index of shared facet.
   *    - In range [0,..., 2*dimension-1].
   *  @param[out] num_moved Number of isosurface vertices moved.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_CUBE_DATA_TYPE,
            typename DUAL_ISOV_TYPE, typename ISOV_INDEX_TYPE,
            typename CUBE_CTYPE, typename ITYPE,
            typename OFFSET_TYPE, typename CTYPE, typename NTYPE>
  void separate_isov_in_cube_from_iso_edge
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   const DUAL_ISOV_TYPE isov_list[],
   const ISOV_INDEX_TYPE isovA0,
   const ISOV_INDEX_TYPE isovA1,
   const CUBE_CTYPE cube_coord[],
   const ITYPE ishared_facet,
   const OFFSET_TYPE offset,
   CTYPE * coord, NTYPE & num_moved)
  {
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename std::vector<GRID_CUBE_DATA_TYPE>::size_type SIZE_TYPE;
    typedef typename GRID_CUBE_DATA_TYPE::NUMBER_TYPE NUMBER_TYPE;

    const int dimension = grid.Dimension();
    const TABLE_INDEX table_index = isov_list[isovA0].table_index;
    const SIZE_TYPE cube0_list_index = isov_list[isovA0].CubeListIndex();
    const ISOV_INDEX_TYPE first_isov =
      active_cube_list[cube0_list_index].first_isov;
    const NUMBER_TYPE num_isov =
      active_cube_list[cube0_list_index].num_isov;
    const int shared_facet_orth_dir =
      IJK::cube_facet_orth_dir(dimension, ishared_facet);
    
    // Initialize
    num_moved = 0;
    
    if (num_isov <= 1) { return; }

    for (NUMBER_TYPE i = 0; i < num_isov; i++) {
      const ISOV_INDEX_TYPE isovB0 = first_isov + i;

      const bool flag_intersects_facet =
        isodual_table_vinfo.VertexInfo(table_index, i).IntersectsFacet(ishared_facet);

      if ((isovB0 != isovA0) && flag_intersects_facet) {
        const bool flag_moved =
          separate_isov_from_iso_edge
          (grid, isodual_table, isodual_table_vinfo,
           isov_list, isovB0, isovA0, isovA1, cube_coord,
           shared_facet_orth_dir, offset, coord);
        if (flag_moved) { num_moved++; }
      }      
    }
    
    return;
  }

  
  /*!
   *  @brief Separate isosurface edges crossing ambiguous facets.
   *  - For each isosurface edge (isovA0,isovA1) do
   *    -# ifacet = Grid facet crossed by (isovA0,isovA).
   *    -# orth_dir = Direction orthgonal to ifacet.
   *    -# c0 = cube containing isovA0.
   *    -# for each vertex isovB0 != isovA0 in cube c0 do
   *      -# for each direction d != orth_dir do
   *        -# Strictly separate isovB0 from isovA0 and isovA1
   *          by plane H orthogonal to d.
   *        -# isovA0, isovA1 and isovB1 should be distance
   *          at least offset from plane H.
   *  @param[out] num_moved Number of isosurface edges moved.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_CUBE_DATA_TYPE,
            typename ISOV_INDEX_TYPE, typename DUAL_ISOV_TYPE,
            typename NTYPEI,
            typename OFFSET_TYPE, typename CTYPE, typename NTYPE>
  void separate_all_iso_edges_crossing_ambiguous_facets
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   const ISOV_INDEX_TYPE isopoly_vert[],
   const NTYPEI num_isopoly,
   const DUAL_ISOV_TYPE isov_list[],
   const OFFSET_TYPE offset,
   CTYPE * coord, NTYPE & num_moved)
  {
    typedef typename DUAL_ISOV_TYPE::CUBE_INDEX_TYPE CUBE_INDEX_TYPE;

    const int dimension = grid.Dimension();
    const IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);
    const NTYPE num_vert_per_isopoly = cube.NumFacetVertices();
    const NTYPE num_facet_edges = cube.NumFacetEdges();
    std::vector<CTYPE> cube_coord(dimension);

    // Note: Unfortunately, this checks/processes each edge twice,
    //   but, typically, few cubes have multiple isosurface vertices,
    //   so it's not worth creating a separate list of edges.
    for (NTYPEI ipoly = 0; ipoly < num_isopoly; ipoly++) {

      // first_vert: First vertex in isosurface polytope ipoly.
      const ISOV_INDEX_TYPE * first_isov =
        isopoly_vert + ipoly*num_vert_per_isopoly;

      for (NTYPE j = 0; j < num_facet_edges; j++) {
        const NTYPE je = cube.FacetEdge(dimension-1, j);
        const NTYPE jend0 = cube.EdgeEndpoint(je, 0);
        const NTYPE jend1 = cube.EdgeEndpoint(je, 1);

        const ISOV_INDEX_TYPE isov0 = first_isov[jend0];
        const ISOV_INDEX_TYPE isov1 = first_isov[jend1];
        const CUBE_INDEX_TYPE icube0 = isov_list[isov0].cube_index;
        const CUBE_INDEX_TYPE icube1 = isov_list[isov1].cube_index;

        // Get facet of icube0 shared with icube1.
        const int ishared_facet0 =
          grid.SharedFacet(icube0, icube1);

        // Separate isosurface edges from isosurface vertices
        //   in cube containing isov0.
        grid.ComputeCoord(icube0, cube_coord.data());
        NTYPE num_moved0;
        separate_isov_in_cube_from_iso_edge
          (grid, isodual_table, isodual_table_vinfo,
           active_cube_list, isov_list, isov0, isov1,
           cube_coord.data(), ishared_facet0, offset,
           coord, num_moved0);
        num_moved += num_moved0;

        const int ishared_facet1 =
          cube.OppositeFacet(ishared_facet0);
                                                      
        // Separate isosurface edges from isosurface vertices
        //   in cube containing isov1.
        grid.ComputeCoord(icube1, cube_coord.data());
        NTYPE num_moved1;
          separate_isov_in_cube_from_iso_edge
          (grid, isodual_table, isodual_table_vinfo,
           active_cube_list, isov_list, isov1, isov0,
           cube_coord.data(), ishared_facet1, offset,
           coord, num_moved1);
      }
    }

  }


  /*!
   *  @overload
   *  @brief Separate isosurface edges crossing ambiguous facets. 
   *    (C++ STL vector isopoly_vert[], isov_list[].)
   *  - C++ STL vector format for arrays isopoly_vert[] and isov_list[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_CUBE_DATA_TYPE,
            typename ISOV_INDEX_TYPE, typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE, typename CTYPE, typename NTYPE>
  void separate_all_iso_edges_crossing_ambiguous_facets
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const OFFSET_TYPE offset,
   CTYPE * coord, NTYPE & num_moved)
  {
    typedef typename std::vector<ISOV_INDEX_TYPE>::size_type SIZE_TYPE;

    const SIZE_TYPE num_vert_per_isopoly = grid.NumCubeFacetVertices();
    const SIZE_TYPE num_isopoly = isopoly_vert.size()/num_vert_per_isopoly;

    separate_all_iso_edges_crossing_ambiguous_facets
      (grid, isodual_table, isodual_table_vinfo, active_cube_list,
       IJK::vector2pointer(isopoly_vert), num_isopoly,
       IJK::vector2pointer(isov_list), offset, coord, num_moved);
  }

  
  /*!
   *  @overload
   *  @brief Separate isosurface edges crossing ambiguous facets. 
   *    (C++ STL vector isopoly_vert[], isov_list[], coord[].)
   *  - C++ STL vector format for arrays isopoly_vert[] and isov_list[].
   *  - C++ STL vector format for array coord[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_CUBE_DATA_TYPE,
            typename ISOV_INDEX_TYPE, typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE, typename CTYPE, typename NTYPE>
  void separate_all_iso_edges_crossing_ambiguous_facets
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord, NTYPE & num_moved)
  {
    separate_all_iso_edges_crossing_ambiguous_facets
      (grid, isodual_table, isodual_table_vinfo, active_cube_list,
       isopoly_vert, isov_list, offset,
       IJK::vector2pointer(coord), num_moved);
  }
  
  //@}
  

  // ******************************************************************
  //! @name Compute location of isosurface vertex on grid edge.
  // ******************************************************************

  //@{

  /*!
   *  @brief Compute position of isosurface vertex on grid edge.
   *  - Compute position using linear interpolation on scalar values
   *    at grid edge endpoints.
   *  - Store coefficient indicating grid edge position in isopoly_info[].
   *  @tparam PTYPE  Precision type to be used in calculations.
   *    - Should be float or double.
   *  @param max_small_difference Maximum small difference.
   *    - If abs(s1-s0) <= max_small_diffence, return 0.5.
   *  @pre max_small_difference >= 0.
   */
  template <typename PTYPE, typename GRID_TYPE, typename STYPE,
	    typename ISOPOLY_INFO_TYPE, typename NTYPE>
  void compute_grid_edge_isov_position_interpolate_scalar
  (const GRID_TYPE & scalar_grid,
   const STYPE isovalue,
   const ISOPOLY_INFO_TYPE isopoly_info[],
   const NTYPE num_poly,
   const PTYPE max_small_difference)
  {
    for (NTYPE ipoly = 0; ipoly < num_poly; ipoly++) {
      isopoly_info[ipoly].grid_edge_isov_position_coef =
        IJK::compute_linear_interpolation_coef_scalar_grid_edge
        (scalar_grid, isopoly_info[ipoly].DualGridEdge(),
         isovalue, max_small_difference);
    }
  }


  //@}


#ifdef INCLUDE_RANDOM

  // *****************************************************************
  //! @name Position iso vertices at random locations.
  // *****************************************************************

  //@{

  /*!
   *  @brief Position dual isosurface vertex at random location in cube.
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   *  @param offset Position isosurface vertices at distance offset
   *      from cube facets.
   *  @param random_engine Pseudorandom number generator.
   *  @param random_pos Random number generator with a given distribution.
   *  @pre 0 <= offset <= 0.5.
   *  @param cube_coord[] Coordinate of lower left vertex of cube.
   */
  template <typename CTYPE_OFFSET,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE0, typename CTYPE1>
  void position_isov_random01
  (const int dimension,
   const CTYPE0 cube_coord[],
   const CTYPE_OFFSET offset,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   CTYPE1 * isov_coord)
  {
    // Put isosurface vertex in random position.
    for (int d = 0; d < dimension; d++) {
      isov_coord[d] =
        cube_coord[d] + random_pos.Random01(offset, random_engine);
    }
  }


  /*!
   *  @brief Position dual isosurface vertices at random location.
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   *  @param offset Position isosurface vertices at distance offset
   *      from cube facets.
   *  @pre 0 <= offset <= 0.5.
   *  @param random_pos Distribution is determined by class random_pos.
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename CTYPE_OFFSET,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
	    typename CTYPE>
  void position_all_isov_random
  (const GRID_TYPE & grid,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const CTYPE_OFFSET offset,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   CTYPE * coord)
  {
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = grid.Dimension();
    IJK::ARRAY<CTYPE> cube_coord(dimension);

    // @pre 0 <= offset <= 0.5, but clamp offset, just in case.
    const CTYPE_OFFSET offsetB =
      IJK::clamp_coord_to_range(offset, 0, 0.5);

    for (SIZE_TYPE isov = 0; isov < isov_list.size(); isov++) {
      const auto icube =
	IJK::get_isovert_grid_cube_index(isov_list[isov]);
      CTYPE * isov_coord = coord + dimension*isov;

      grid.ComputeCoord(icube, cube_coord.Ptr());

      // Put isosurface vertex in random position.
      position_isov_random01
        (dimension, cube_coord.PtrConst(), offsetB, random_engine, random_pos,
         isov_coord);
    }
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices at random location.
   *  - C++ STL vector format for array coord[].
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   */
  template <typename GRID_TYPE, typename VTYPE,
            typename CTYPE_OFFSET, typename RANDOM_ENGINE,
	    typename RANDOM_POS_TYPE, typename CTYPE>
  void position_all_isov_random
  (const GRID_TYPE & grid,
   const std::vector<VTYPE> & vlist,
   const CTYPE_OFFSET offset,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   std::vector<CTYPE> & coord)
  {
    const int dimension = grid.Dimension();

    coord.resize(vlist.size()*dimension);
    position_all_isov_random
      (grid, vlist, offset, random_engine, random_pos,
       IJK::vector2pointerNC(coord));
  }


  /*!
   *  @brief Position dual isosurface vertices at random location.
   *    (Initalize random seed.)
   *  - Initialize random engine with seed.
   *  - C++ STL vector format for array coord[].
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   */
  template <typename GRID_TYPE, typename VTYPE, typename CTYPE_OFFSET,
            typename RANDOM_SEED, typename RANDOM_POS_TYPE,
	    typename CTYPE>
  void position_all_isov_random_init
  (const GRID_TYPE & grid,
   const std::vector<VTYPE> & vlist,
   const CTYPE_OFFSET offset,
   const RANDOM_SEED random_seed,
   RANDOM_POS_TYPE & random_pos,
   std::vector<CTYPE> & coord)
  {
    std::minstd_rand random_engine(random_seed);

    position_all_isov_random
      (grid, vlist, offset, random_engine, random_pos, coord);
  }


  /*!
    *  @brief Position dual isosurface vertices at random location.
    *  - Use uniform distribution.
    *  - C++ STL vector format for array coord[].
    *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
    */
   template <typename GRID_TYPE, typename VTYPE,
             typename CTYPE_OFFSET, typename RANDOM_SEED,
	     typename CTYPE>
   void position_all_isov_random_uniform
   (const GRID_TYPE & grid,
    const std::vector<VTYPE> & vlist,
    const CTYPE_OFFSET offset,
    const RANDOM_SEED random_seed,
    std::vector<CTYPE> & coord)
   {
     IJK::GENERATE_UNIFORM_RANDOM random_pos;

     position_all_isov_random_init
       (grid, vlist, offset, random_seed, random_pos, coord);
   }


  /*!
   *  @brief Position dual isosurface vertices at random location.
   *  - Use U-quadratic distribution.
   *  - C++ STL vector format for array coord[].
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   */
  template <typename GRID_TYPE, typename VTYPE,
            typename CTYPE_OFFSET, typename RANDOM_SEED,
	    typename CTYPE>
  void position_all_isov_random_U_quadratic
  (const GRID_TYPE & grid,
   const std::vector<VTYPE> & vlist,
   const CTYPE_OFFSET offset,
   const RANDOM_SEED random_seed,
   std::vector<CTYPE> & coord)
  {
    IJK::GENERATE_U_QUADRATIC_RANDOM random_pos;

    position_all_isov_random_init
      (grid, vlist, offset, random_seed, random_pos, coord);
  }


  /*!
   *  @brief Position dual isosurface vertex at random location in cube.
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   *  @param cube_coord[] Coordinate of lower left vertex of cube
   *    containing isosurface vertex.
   *  @param bwidth Boundary width. Controls boundary sampling.
   *    - Generate random number in range [-boundary_width,boundary_width]
   *      and then clamp to [0,1].
   *  @param random_engine Pseudorandom number generator.
   *  @param random_pos Random number generator with a given distribution.
   */
  template <typename BWIDTH,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE0, typename CTYPE1>
  void position_isov_random01_sample_boundary
  (const int dimension,
   const CTYPE0 cube_coord[],
   const BWIDTH boundary_width,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   CTYPE1 * isov_coord)
  {
     // Put isosurface vertex in random position.
     for (int d = 0; d < dimension; d++) {
       const CTYPE1 rpos01 = random_pos.RandomAndClamp
         (0, 1, boundary_width, random_engine);

       isov_coord[d] = cube_coord[d] + rpos01;
     }
   }


  /*!
   *  @brief Position dual isosurface vertex at random location in cube.
   *  - Includes offset. Boundary vertices are sampled on the offset.
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   *  @param cube_coord[] Coordinate of lower left vertex of cube
   *    containing isosurface vertex.
   *  @param offset Offset.
   *  @param bwidth Boundary width. Controls boundary sampling.
   *    - Generate random number in range
   *      [offset-boundary_width,(1-offset)+boundary_width]and
   *      then clamp to [offset,1-offset].
   *  @param random_engine Pseudorandom number generator.
   *  @param random_pos Random number generator with a given distribution.
    */
  template <typename CTYPE_OFFSET, typename BWIDTH,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE0, typename CTYPE1>
  void position_isov_random_offset_sample_boundary
  (const int dimension,
   const CTYPE0 cube_coord[],
   const CTYPE_OFFSET offset,
   const BWIDTH boundary_width,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   CTYPE1 * isov_coord)
   {
     // Put isosurface vertex in random position.
     for (int d = 0; d < dimension; d++) {
       const CTYPE1 rpos = random_pos.RandomAndClamp
         (offset, 1.0-offset, boundary_width, random_engine);

       isov_coord[d] = cube_coord[d] + rpos;
     }
   }


  /*!
   *  @brief Position dual isosurface vertices at random location.
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   *  @param bwidth Boundary width. Controls boundary sampling.
   *     - Generate random number in range [-boundary_width,boundary_width]
   *       and then clamp to [0,1].
   *  @param random_pos Distribution is determined by class random_pos.
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE, typename BWIDTH,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_all_isov_random_sample_boundary
  (const GRID_TYPE & grid,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const BWIDTH boundary_width,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   CTYPE * coord)
  {
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = grid.Dimension();
    std::vector<CTYPE> cube_coord(dimension);
    CTYPE * const cube_coord_ptr = IJK::vector2pointerNC(cube_coord);

    for (SIZE_TYPE isov = 0; isov < isov_list.size(); isov++) {
      CTYPE * const isov_coord = coord + isov*dimension;
      const auto icube =
        IJK::get_isovert_grid_cube_index(isov_list[isov]);
      grid.ComputeCoord(icube, cube_coord_ptr);

      position_isov_random01_sample_boundary
        (dimension, cube_coord_ptr, boundary_width,
         random_engine, random_pos, isov_coord);
    }
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices at random location.
   *  - C++ STL vector format for array coord[].
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   */
  template <typename GRID_TYPE, typename VTYPE,
            typename BWIDTH, typename RANDOM_ENGINE,
            typename RANDOM_POS_TYPE, typename CTYPE>
  void position_all_isov_random_sample_boundary
  (const GRID_TYPE & grid,
   const std::vector<VTYPE> & vlist,
   const BWIDTH boundary_width,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   std::vector<CTYPE> & coord)
  {
    const int dimension = grid.Dimension();

    coord.resize(vlist.size()*dimension);
    position_all_isov_random_sample_boundary
    (grid, vlist, boundary_width, random_engine, random_pos,
     IJK::vector2pointerNC(coord));
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices at random location.
   *  - Initialize random engine with seed.
   *  - C++ STL vector format for array coord[].
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   */
  template <typename GRID_TYPE, typename VTYPE, typename BWIDTH,
            typename RANDOM_SEED, typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_all_isov_random_sample_boundary_init
  (const GRID_TYPE & grid,
   const std::vector<VTYPE> & vlist,
   const BWIDTH boundary_width,
   const RANDOM_SEED random_seed,
   RANDOM_POS_TYPE & random_pos,
   std::vector<CTYPE> & coord)
  {
    std::minstd_rand random_engine(random_seed);

    position_all_isov_random_sample_boundary
    (grid, vlist, boundary_width, random_engine, random_pos, coord);
  }

  //@}
  

  // *****************************************************************
  //! @name Position at random locations, flag centroid positioning.
  // *****************************************************************

  // NEED TO PASS position_offset AND iso_grid_edgeI_offset.
  /*!
   *  @brief Position dual isosurface vertices at random location.
   *      (Use centroid positioning.)
   *    - Use centroid positioning in all ambiguous cubes
   *      or cubes with ambiguous facets.
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   *  @param offset Position isosurface vertices at distance offset
   *      from cube facets.
   *  @pre 0 <= offset <= 0.5.
   *  @param random_pos Distribution is determined by class random_pos.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename DUAL_ISOV_TYPE,
            typename CTYPE_OFFSET,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_all_isov_multi_random_use_centroid_pos
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const CTYPE_OFFSET offset,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   CTYPE * coord)
  {
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();
    IJK::ARRAY<CTYPE> cube_coord(dimension);
    IJK::ARRAY<CTYPE> temp_coord(dimension);

    // @pre 0 <= offset <= 0.5, but clamp offset, just in case.
    const CTYPE_OFFSET offsetB =
      IJK::clamp_coord_to_range(offset, 0, 0.5);

    for (SIZE_TYPE isov = 0; isov < isov_list.size(); isov++) {

      const TABLE_INDEX table_index = isov_list[isov].table_index;
      CTYPE * isov_coord = coord + dimension*isov;
      
      if (isodual_table.IsAmbiguous(table_index) ||
          (isodual_table.NumAmbiguousFacets(table_index) > 0)) {
        IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);

        position_isov_in_cube_centroid_offset_multi_checkI
          (scalar_grid, isodual_table, isovalue, isov_list[isov],
           offsetB, cube, isov_coord, temp_coord.Ptr());
      }
      else {
        const auto icube =
          IJK::get_isovert_grid_cube_index(isov_list[isov]);

        scalar_grid.ComputeCoord(icube, cube_coord.Ptr());

        // Put isosurface vertex in random position.
        position_isov_random01
          (dimension, cube_coord.PtrConst(), offsetB,
           random_engine, random_pos, isov_coord);
      }
    }
  }


  // NEED TO PASS position_offset AND iso_grid_edgeI_offset.
  /*!
   *  @overload
   *  @brief Position dual isosurface vertices at random location.
   *    Use centroid positioning. (C++ vector.)
   *  - C++ STL vector format for array coord[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename DUAL_ISOV_TYPE,
            typename CTYPE_OFFSET,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_all_isov_multi_random_use_centroid_pos
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const CTYPE_OFFSET offset,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   std::vector<CTYPE> & coord)
  {
    position_all_isov_multi_random_use_centroid_pos
      (scalar_grid, isodual_table, isovalue, isov_list, offset,
       random_engine, random_pos, IJK::vector2pointerNC(coord));
  }

  
  // NEED TO PASS position_offset AND iso_grid_edgeI_offset.
  /*!
   *  @brief Position dual isosurface vertices at random location.
   *      Use centroid positioning. Sample boundary.
   *    - Sample boundary of offset region.
   *    - Use centroid positioning in all ambiguous cubes
   *      or cubes with ambiguous facets.
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   *  @param offset Position isosurface vertices at distance offset
   *      from cube facets.
   *  @param boundary_width Boundary width.
   *  - Generate coordinates in range [rmin-boundary_width,rmax+boundary_width]
   *    and then clamp to [rmin,rmax].   
   *  @pre 0 <= offset <= 0.5.
   *  @param random_pos Distribution is determined by class random_pos.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename DUAL_ISOV_TYPE,
            typename CTYPE_OFFSET, typename BWIDTH,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_all_isov_multi_random_sampleB_use_centroid_pos
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const CTYPE_OFFSET offset,
   const BWIDTH boundary_width,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   CTYPE * coord)
  {
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();
    IJK::ARRAY<CTYPE> cube_coord(dimension);
    IJK::ARRAY<CTYPE> temp_coord(dimension);

    // @pre 0 <= offset <= 0.5, but clamp offset, just in case.
    const CTYPE_OFFSET offsetB =
      IJK::clamp_coord_to_range(offset, 0, 0.5);

    for (SIZE_TYPE isov = 0; isov < isov_list.size(); isov++) {

      const TABLE_INDEX table_index = isov_list[isov].table_index;
      CTYPE * isov_coord = coord + dimension*isov;
      
      if (isodual_table.IsAmbiguous(table_index) ||
          (isodual_table.NumAmbiguousFacets(table_index) > 0)) {
        IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);

        position_isov_in_cube_centroid_offset_multi_checkI
          (scalar_grid, isodual_table, isovalue, isov_list[isov],
           offsetB, cube, isov_coord, temp_coord.Ptr());
      }
      else {
        const auto icube =
          IJK::get_isovert_grid_cube_index(isov_list[isov]);

        scalar_grid.ComputeCoord(icube, cube_coord.Ptr());

        // Put isosurface vertex in random position.
        position_isov_random_offset_sample_boundary
          (dimension, cube_coord.PtrConst(), offsetB, boundary_width,
           random_engine, random_pos, isov_coord);
      }
    }
  }


  // NEED TO PASS position_offset AND iso_grid_edgeI_offset.
  /*!
   *  @overload
   *  @brief Position dual isosurface vertices at random location.
   *    Use centroid positioning. (C++ vector.)
   *  - C++ STL vector format for array coord[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename DUAL_ISOV_TYPE,
            typename CTYPE_OFFSET, typename BTYPE,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_all_isov_multi_random_sampleB_use_centroid_pos
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const CTYPE_OFFSET offset,
   const BTYPE boundary_width,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   std::vector<CTYPE> & coord)
  {
    position_all_isov_multi_random_sampleB_use_centroid_pos
      (scalar_grid, isodual_table, isovalue, isov_list,
       offset, boundary_width, random_engine, random_pos,
       IJK::vector2pointerNC(coord));
  }


  // NEED TO PASS position_offset AND iso_grid_edgeI_offset.
  /*!
   *  @overload
   *  @brief Position dual isosurface vertices at random location.
   *      (Use centroid positioning.)
   *    - Use centroid positioning in all ambiguous cubes
   *      or cubes with ambiguous facets.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename DUAL_ISOV_TYPE,
            typename CTYPE_OFFSET, typename BWIDTH,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_all_isov_multi_random_use_centroid_pos
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const CTYPE_OFFSET offset,
   const BWIDTH boundary_width,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   std::vector<CTYPE> & coord)
  {
    if (boundary_width == 0) {
      position_all_isov_multi_random_use_centroid_pos
        (scalar_grid, isodual_table, isovalue, isov_list, offset,
         random_engine, random_pos, coord);
    }
    else {
      position_all_isov_multi_random_sampleB_use_centroid_pos
        (scalar_grid, isodual_table, isovalue, isov_list,
         offset, boundary_width, random_engine, random_pos, coord);
    }
  }
  
  // *****************************************************************
  //! @name Position at random locations, seperate isov by cube centers.
  // *****************************************************************

  //@{
  
  /*!
   *  @brief Position dual isosurface vertices at random location
   *  - More than one vertex can be in a cube.
   *  - If cube contains multiple isosurface then vertices are positioned
   *    in separate halves/quarters/eighths of cube.
   *  - In 3D, this routine calls
   *    reposition_all_doubly_connected_dual_isov_cube3D().
   *  - Use position_all_isov_random_uniform()
   *    or position_all_isov_random_sample_boundary_init()
   *    to generate random positions without automatically calling
   *    reposition_all_doubly_connected_dual_isov_cube3D().
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   *  @param isodual_table_vinfo Isodual table vertex information.
   *    @pre compute_dual_cube_isotable_vertex_connectivity()
   *       has been called on isodual_table_vinfo.
   *  @param offset Position isosurface vertices at distance offset
   *      from cube facets.
   *    @pre 0 <= offset <= 0.25.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_CUBE_DATA_TYPE, typename NTYPEC,
            typename DUAL_ISOV_TYPE, typename CTYPE_OFFSET,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_all_isov_multi_random_sep_by_cc
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const GRID_CUBE_DATA_TYPE active_cube_list[],
   const NTYPEC num_active_cubes,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const CTYPE_OFFSET offset,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   CTYPE * coord)
  {
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;

    const int DIM3(3);
    const int dimension = grid.Dimension();
    IJK::ARRAY<CTYPE> cube_coord(dimension);
    IJK::PROCEDURE_ERROR error
      ("position_all_isov_multi_random_sep_by_cc");

    if (offset < 0.0 || offset > 0.25) {
      error.AddMessage("Programming error. Offset not in range [0,0.25].");
      throw error;
    }

    for (SIZE_TYPE isov = 0; isov < isov_list.size(); isov++) {
      const VTYPE icube = isov_list[isov].cube_index;
      const NTYPE ipatch = isov_list[isov].patch_index;
      const TABLE_INDEX it = isov_list[isov].table_index;

      CTYPE * isov_coord = coord + dimension*isov;

      grid.ComputeCoord(icube, cube_coord.Ptr());

      if (isodual_table.NumIsoVertices(it) == 1) {

        // Put isosurface vertex in random position.
        for (int d = 0; d < dimension; d++) {
          isov_coord[d] =
            cube_coord[d] + random_pos.Random01(offset, random_engine);
        }

      }
      else {

        NTYPE num_orth_dir_with_two_intersected_facets = 0;
        for (int d = 0; d < dimension; d++) {
          if (isodual_table_vinfo.VertexInfo(it,ipatch).IntersectsBothFacetsOrthogonalTo
              (d, dimension))
            { num_orth_dir_with_two_intersected_facets++; }
        }

        for (int d = 0; d < dimension; d++) {

          // At most one of these two flags should be true.
          bool flag_lower_pos = false;
          bool flag_upper_pos = false;

          // Note: This may not guarantee that all self intersections
          //   are avoided in 4D.
          const bool flag_intersects_lower_facet =
            isodual_table_vinfo.VertexInfo(it, ipatch).IntersectsFacet(d);
          const bool flag_intersects_upper_facet =
            isodual_table_vinfo.VertexInfo(it, ipatch).IntersectsFacet(d+dimension);

          if (flag_intersects_lower_facet != flag_intersects_upper_facet) {
            flag_lower_pos = flag_intersects_lower_facet;
            flag_upper_pos = flag_intersects_upper_facet;
          }
          else {
            const NTYPE ipatch1 =
              (ipatch + 1)%(isodual_table.NumIsoVertices(it));

            // Use ipatch1 to determine position of vertex in ipatch0.
            // Note: This does not correctly handle all cases in 4D.
            // In particular, there may be 3 iso vertices,
            //   and the vertex should be positioned between two others.
            const bool flag_intersects_lower_facet =
              isodual_table_vinfo.VertexInfo(it, ipatch1).IntersectsFacet(d);
            const bool flag_intersects_upper_facet =
              isodual_table_vinfo.VertexInfo(it, ipatch1).IntersectsFacet(d+dimension);

            if (flag_intersects_lower_facet != flag_intersects_upper_facet) {
              flag_lower_pos = !flag_intersects_lower_facet;
              flag_upper_pos = !flag_intersects_upper_facet;
            }            
          }
          
          if (flag_lower_pos) {
            isov_coord[d] =
              cube_coord[d] +
              random_pos.Random(0.0, 0.5, offset, random_engine);
          }
          else if (flag_upper_pos) {
            isov_coord[d] =
              cube_coord[d] +
              random_pos.Random(0.5, 1.0, offset, random_engine);
          }
          else {
            isov_coord[d] =
              cube_coord[d] + random_pos.Random01(offset, random_engine);
          }

        }
      }
    }

    if (dimension == DIM3) {
      NTYPE num_moved;
      reposition_all_doubly_connected_3D_centerIII
        (grid, isodual_table, active_cube_list, num_active_cubes,
         offset, coord, num_moved);
    }
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices at random location. (C++ vector.)
   *  - Version using C++ STL vector for arrays active_cube_list[] and coord[].
   *  - More than one vertex can be in a cube.
   *  - If cube contains multiple isosurface then vertices are positioned
   *    in separate halves/quarters/eighths of cube.
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
 	    typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_CUBE_DATA_TYPE, typename DUAL_ISOV_TYPE,
	    typename CTYPE_OFFSET,
	    typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
	    typename CTYPE>
  void position_all_isov_multi_random_sep_by_cc
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const CTYPE_OFFSET offset,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   std::vector<CTYPE> & coord)
   {
     typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

     const DTYPE dimension = grid.Dimension();

     coord.resize(isov_list.size()*dimension);
     position_all_isov_multi_random_sep_by_cc
       (grid, isodual_table, isodual_table_vinfo,
        IJK::vector2pointer(active_cube_list), active_cube_list.size(),
        isov_list, offset, random_engine, random_pos,
        IJK::vector2pointerNC(coord));
   }


  /*!
    *  @brief Position dual isosurface vertices at random location.
    *    (C++ vector and create random_engine.)
    *  - Version using C++ STL vector for arrays active_cube_list[] 
    *      and coord[].
    *  - Version that creates random_engine initialized with random_seed.
    *  - More than one vertex can be in a cube.
    *  - If cube contains multiple isosurface then vertices are positioned
    *    in separate halves/quarters/eighths of cube.
    *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
    */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
	    typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_CUBE_DATA_TYPE, typename DUAL_ISOV_TYPE,
 	    typename SEED_TYPE, typename CTYPE_OFFSET,
 	    typename RANDOM_POS_TYPE, typename CTYPE>
   void position_all_isov_multi_random_sep_by_cc_init
   (const GRID_TYPE & grid,
    const ISODUAL_TABLE_TYPE & isodual_table,
    const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
    const std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
    const std::vector<DUAL_ISOV_TYPE> & isov_list,
    const CTYPE_OFFSET offset,
    const SEED_TYPE random_seed,
    RANDOM_POS_TYPE & random_pos,
    std::vector<CTYPE> & coord)
  {
    std::minstd_rand random_engine(random_seed);

    position_all_isov_multi_random_sep_by_cc
    (grid, isodual_table, isodual_table_vinfo,
     active_cube_list, isov_list, offset,
     random_engine, random_pos, coord);
  }

  //@}
  

  // *****************************************************************
  //! @name Position at random locations, separate isov by degree.
  // *****************************************************************

  //@{
  
  /*!
   *  @brief Generate a random coordinate in range.
   *  - If flag_lower_pos is true, generate random coordinate
   *    in range [offset, lowerB-offset].
   *  - Otherwise, generate random coordinate
   *    in range [upperA+offset, 1.0].
   *  @param offset Range offset.
   *    @pre 2*offset <= 2*lowerB.
   *    @pre 2*offset <= 1-upperA.
   */
  template <typename CTYPE, typename CTYPER, typename OFFSET_TYPE,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE>
  CTYPE random_coord_in_range
  (const bool flag_lower_pos,
   const CTYPER lowerB,
   const CTYPER upperA,
   const OFFSET_TYPE offset,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos)
  {
    if (flag_lower_pos)
      { return random_pos.Random(0.0, lowerB, offset, random_engine); }
    else
      { return random_pos.Random(upperA, 1.0, offset, random_engine); }
  }


  /*!
   *  @brief Position degree 3 isosurface vertex in random location in a cube.
   *  - Position degree 3 vertex in (1/3)x(1/3)x(1/3) subcube.
   *  - NOTE: Only for dimension 3.
   *  @param isov Index of isosurface vertex.
   *    @pre Isosurface vertex isov has degree 3.
   *  @param ipatch Patch index of isosurface vertex isov.
   *  @isov_coord[] Coordinates of isosurface vertex isov.
   */
  template <typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_COORD_TYPE, typename TABLE_INDEX_TYPE,
            typename CTYPE_OFFSET,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_degree3_isov_random_3D
  (const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const GRID_COORD_TYPE cube_coord[],
   const TABLE_INDEX_TYPE table_index,
   const int ipatch,
   const CTYPE_OFFSET offset,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   CTYPE isov_coord[])
  {
    const int DIM3(3);
    const CTYPE ONE_THIRD(1.0/3.0);
    const CTYPE TWO_THIRDS(2.0/3.0);
    
    for (int d = 0; d < DIM3; d++) {
      // Since isov degree is 3, some edge incident on isov
      //   intersects lower facet or upper facet, but not both.
      const bool flag_intersects_lower_facet =
        isodual_table_vinfo.VertexInfo(table_index, ipatch).IntersectsFacet(d);

      const CTYPE c = random_coord_in_range<CTYPE>
        (flag_intersects_lower_facet, ONE_THIRD, TWO_THIRDS, offset,
         random_engine, random_pos);
      
      isov_coord[d] = cube_coord[d] + c;
    }
  }


  /*!
   *  @brief Position degree 4 isosurface vertex in random location
   *    in a cube.
   *  - Degree 4 vertices are either positioned on an axis-oriented
   *    square incident on the cube center or on an axis-oriented
   *    line segment, depending on the facets crossed by their
   *    incident edges.
   *  - NOTE: Only for dimension 3.
   *  @param isov Index of isosurface vertex.
   *    @pre Isosurface vertex isov has degree 4.
   *  @param ipatch Patch index of isosurface vertex isov.
   *  @isov_coord[] Coordinates of isosurface vertex isov.
   */
  template <typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_COORD_TYPE,
            typename ISOV_INDEX_TYPE, typename TABLE_INDEX_TYPE,
            typename CTYPE_OFFSET,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_degree4_isov_random_3D
  (const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const GRID_COORD_TYPE cube_coord[],
   const ISOV_INDEX_TYPE isov,
   const TABLE_INDEX_TYPE table_index,
   const int ipatch,
   const CTYPE_OFFSET offset,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   CTYPE isov_coord[])
  {
    const int DIM3(3);
    const CTYPE ONE_HALF(1.0/2.0);

    for (int d = 0; d < DIM3; d++) {
      const bool flag_intersects_lower_facet =
        isodual_table_vinfo.VertexInfo(table_index, ipatch).IntersectsFacet(d);
      const bool flag_intersects_upper_facet =
        isodual_table_vinfo.VertexInfo(table_index, ipatch).IntersectsFacet(d+DIM3);

      if (flag_intersects_lower_facet && flag_intersects_upper_facet) {
        // Isosurface vertex is halfway between lower and upper facets.
        isov_coord[d] = cube_coord[d] + 0.5;
      }
      else if (flag_intersects_lower_facet || flag_intersects_upper_facet) {
        // Incident edges intersect one facet but not the other.
        const CTYPE c = random_coord_in_range<CTYPE>
          (flag_intersects_lower_facet, ONE_HALF, ONE_HALF, offset,
           random_engine, random_pos);
      
        isov_coord[d] = cube_coord[d] + c;
      }
      else {
        // No incident edge intersects facet orthogonal to d.
        // Put isosurface vertex in random position.
        isov_coord[d] =
          cube_coord[d] + random_pos.Random01(offset, random_engine);
      }
    }
  }


  /*!
   *  @brief Position degree 5 or 6 isosurface vertex in random location
   *    in a cube.
   *  - If cube has two vertices, position degree 5 or 6 vertex
   *    in (2/3)x(2/3)x(2/3) subcube that does not contain degree 3 vertex.
   *  - Otherwise, position vertex anywhere in cube.
   *  - NOTE: Only for dimension 3.
   *  @pre Cube has at most 2 vertices.
   *  @param isov Index of isosurface vertex.
   *    @pre Isosurface vertex isov has degree 4.
   *  @param ipatch Patch index of isosurface vertex isov.
   *    @pre ipatch = 0 or 1.
   *  @isov_coord[] Coordinates of isosurface vertex isov.
   */
  template <typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_COORD_TYPE, typename GRID_CUBE_DATA_TYPE,
            typename CTYPE_OFFSET,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_degree5_or_degree6_isov_random_3D
  (const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const GRID_COORD_TYPE cube_coord[],
   const GRID_CUBE_DATA_TYPE & active_cube,
   const int ipatch,
   const CTYPE_OFFSET offset,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   CTYPE isov_coord[])
  {
    typedef typename GRID_CUBE_DATA_TYPE::TABLE_INDEX_TYPE TABLE_INDEX_TYPE;
    
    const int DIM3(3);
    const CTYPE ONE_THIRD(1.0/3.0);
    const CTYPE TWO_THIRDS(2.0/3.0);
    
    if (active_cube.num_isov == 2) {
      const TABLE_INDEX_TYPE table_index = active_cube.table_index;
      const int ipatch_other = 1-ipatch;
      
      for (int d = 0; d < DIM3; d++) {
        const bool flag_intersects_lower_facet =
          isodual_table_vinfo.VertexInfo(table_index, ipatch_other).IntersectsFacet(d);

        // Move isov to region NOT containing ipatch_other.
        // Region has size (2/3)x(2/3)x(2/3).
        const CTYPE c = random_coord_in_range<CTYPE>
          (!flag_intersects_lower_facet, TWO_THIRDS, ONE_THIRD, offset,
           random_engine, random_pos);
      
        isov_coord[d] = cube_coord[d] + c;
      }
    }
    else if (active_cube.num_isov == 1) {
      // Put isosurface vertex in random position.
      position_isov_random01
        (DIM3, cube_coord, offset, random_engine, random_pos,
         isov_coord);      
    }
    else {
      IJK::PROCEDURE_ERROR
        error("position_degree5_or_degree6_isov_random_3D");
      
      if (active_cube.num_isov == 0) {
        error.AddMessage
          ("Programming error. Active cube contains no isosurface vertices.");
        throw error;
      }
      else {
        error.AddMessage
          ("Programming error. Active cube contains more than 2 isosurface vertices.");
        throw error;        
      }
    }

  }


  /*!
   *  @brief Position two dual isosurface vertices in random location in a cube.
   *  - Cube contains exactly two isosurface vertices.
   *  - Position based on degrees.
   *  - Degree 3 vertices are placed in (1/3)x(1/3)x(1/3) subcube.
   *  - Degree 4 vertices are placed on an axis-oriented square
   *    incident on the cube center.
   *  - Degree 5 and 6 vertices are placed in (2/3)x(2/3)x(2/3) subcubes.
   *  - NOTE: Only for dimension 3.
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   *  @param isodual_table_vinfo Isodual table vertex information.
   *    @pre compute_dual_cube_isotable_vertex_connectivity()
   *       has been called on isodual_table_vinfo.
   *  @param active_cube Information about grid cube containing
   *     the 2 isosurface vertices.
   *    @pre active_cube.num_isov == 2.
   *  @param offset Position isosurface vertices at distance offset
   *      from cube facets.
   *    @pre 0 <= offset <= 1/6.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_CUBE_DATA_TYPE, typename CTYPE_OFFSET,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_isovII_random_sep_by_degree_3D
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const GRID_CUBE_DATA_TYPE & active_cube,
   const CTYPE_OFFSET offset,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   CTYPE * coord)  
  {
    typedef typename GRID_CUBE_DATA_TYPE::CUBE_INDEX_TYPE CUBE_INDEX_TYPE;
    typedef typename GRID_CUBE_DATA_TYPE::ISOV_TYPE ISOV_TYPE;
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;

    const int DIM3(3);
    const TABLE_INDEX table_index = active_cube.table_index;
    const CUBE_INDEX_TYPE icube = active_cube.cube_index;
    CTYPE cube_coord[DIM3];

    grid.ComputeCoord(icube, cube_coord);
    
    for (int ipatch = 0; ipatch < active_cube.num_isov; ipatch++) {
      const int degree =
        isodual_table_vinfo.VertexInfo(table_index,ipatch).degree;
      const ISOV_TYPE isov = active_cube.first_isov + ipatch;
      CTYPE * isov_coord = coord + DIM3*isov;

      // *** DEBUG ***
      /*
      using namespace std;
      cerr << "Reposition degree " << degree << " isov " << isov;
      grid.PrintIndexAndCoord(cerr, " in cube ", icube, ".\n");
      const int table_index = active_cube.table_index;
      cerr << "  table_index: " << table_index;
      cerr << "  ambiguous facet: "
           << isodual_table.IndexOfSomeAmbiguousFacet(table_index)
           << endl;
      */
      
      if (degree == 3) {
        position_degree3_isov_random_3D
          (isodual_table_vinfo, cube_coord, table_index, ipatch,
           offset, random_engine, random_pos, isov_coord);        
      }
      else if (degree == 4) {
        position_degree4_isov_random_3D
          (isodual_table_vinfo, cube_coord, isov, table_index,
           ipatch, offset, random_engine, random_pos, isov_coord);        
      }
      else {
        // degree = 5 or 6 and 1-ipatch has degree 3.
        position_degree5_or_degree6_isov_random_3D
          (isodual_table_vinfo, cube_coord, active_cube, ipatch,
           offset, random_engine, random_pos, isov_coord);        
      }
    }
  }




  /*!
   *  @brief Position dual isosurface vertices at random location.
   *  - More than one vertex can be in a cube.
   *  - If cube contains multiple isosurface then vertices are positioned
   *    separated based on degree.
   *  - Degree 3 vertices are placed in (1/3)x(1/3)x(1/3) subcube.
   *  - Degree 4 vertices are placed on an axis-oriented square
   *    incident on the cube center.
   *  - Degree 5 and 6 vertices are placed in (2/3)x(2/3)x(2/3) subcubes.
   *  - NOTE: Only for dimension 3.
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   *  @param isodual_table_vinfo Isodual table vertex information.
   *    @pre compute_dual_cube_isotable_vertex_connectivity()
   *       has been called on isodual_table_vinfo.
   *  @param offset Position isosurface vertices at distance offset
   *      from cube facets.
   *    @pre 0 <= offset <= 1/6.
   *  @param boundary_width Boundary width.
   *  - Generate coordinates in range [rmin-boundary_width,rmax+boundary_width]
   *    and then clamp to [rmin,rmax].   
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_CUBE_DATA_TYPE, typename NTYPEC,
            typename CTYPE_OFFSET,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_all_isov_multi_random_sep_by_degree_3D
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const GRID_CUBE_DATA_TYPE active_cube_list[],
   const NTYPEC num_active_cubes,
   const CTYPE_OFFSET offset,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   CTYPE * coord)
  {
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename GRID_CUBE_DATA_TYPE::ISOV_TYPE ISOV_TYPE;
    
    const int DIM3(3);
    const CTYPE ONE_SIXTH(1.0/6.0);
    CTYPE cube_coord[DIM3];

    IJK::PROCEDURE_ERROR error
      ("position_all_isov_multi_random_sep_by_degree_3D");

    if (grid.Dimension() != DIM3) {
      error.AddMessage("Programming error. Dimension not equal to 3.");
      error.AddMessage("  grid.Dimension(): ", grid.Dimension());
      throw error;
    }
    
    if (offset < 0.0 || offset > ONE_SIXTH) {
      error.AddMessage("Programming error. Offset not in range [0,1/6].");
      throw error;
    }

    for (NTYPEC i = 0; i < num_active_cubes; i++) {
      const TABLE_INDEX table_index =
        active_cube_list[i].table_index;
      const int num_isov_in_cube =
        isodual_table.NumIsoVertices(table_index);
      const VTYPE icube = active_cube_list[i].cube_index;
      
      grid.ComputeCoord(icube, cube_coord);
      
      if (num_isov_in_cube == 1) {
        const ISOV_TYPE isov =
          active_cube_list[i].first_isov;
        CTYPE * isov_coord = coord + DIM3*isov;
          
        // Put isosurface vertex in random position.
        position_isov_random01
          (DIM3, cube_coord, offset, random_engine, random_pos,
           isov_coord);
      }
      else if (num_isov_in_cube == 2) {
        position_isovII_random_sep_by_degree_3D
          (grid, isodual_table, isodual_table_vinfo, active_cube_list[i],
           offset, random_engine, random_pos, coord);
      }
      else {
        for (int ipatch = 0; ipatch < num_isov_in_cube; ipatch++) {
          const ISOV_TYPE isov =
            active_cube_list[i].first_isov + ipatch;
          CTYPE * isov_coord = coord + DIM3*isov;

          position_degree3_isov_random_3D
            (isodual_table_vinfo, cube_coord, table_index, ipatch,
             offset, random_engine, random_pos, isov_coord);
        }
      }
    }
  }

  //@}


  // *****************************************************************
  //! @name Position at random locations, sample boundary, separate isov by degree.
  // *****************************************************************

  //@{

  /*!
   *  @brief Generate a random coordinate in range.
   *  - Sample boundary.
   *  - If flag_lower_pos is true, generate random coordinate
   *    in range [offset, lowerB-offset].
   *  - Otherwise, generate random coordinate
   *    in range [upperA+offset, 1.0].
   *  @param offset Range offset.
   *    @pre 2*offset <= 2*lowerB.
   *    @pre 2*offset <= 1-upperA.
   */
  template <typename CTYPE, typename CTYPER,
            typename OFFSET_TYPE, typename BWIDTH,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE>
  CTYPE random_coord_in_range_sampleB
  (const bool flag_lower_pos,
   const CTYPER lowerB,
   const CTYPER upperA,
   const OFFSET_TYPE offset,
   const BWIDTH boundary_width,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos)
  {
    if (flag_lower_pos) {
      return random_pos.RandomAndClamp
        (offset, lowerB-offset, boundary_width, random_engine);
    }
    else {
      return random_pos.RandomAndClamp
        (upperA+offset, 1.0-offset, boundary_width, random_engine);
    }
  }


  /*!
   *  @brief Position degree 3 isosurface vertex in random location in a cube.
   *  - Position degree 3 vertex in (1/3)x(1/3)x(1/3) subcube.
   *  - Sample boundary.
   *  - NOTE: Only for dimension 3.
   *  @param isov Index of isosurface vertex.
   *    @pre Isosurface vertex isov has degree 3.
   *  @param ipatch Patch index of isosurface vertex isov.
   *  @isov_coord[] Coordinates of isosurface vertex isov.
   */
  template <typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_COORD_TYPE, typename TABLE_INDEX_TYPE,
            typename CTYPE_OFFSET, typename BWIDTH,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_degree3_isov_random_sampleB_3D
  (const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const GRID_COORD_TYPE cube_coord[],
   const TABLE_INDEX_TYPE table_index,
   const int ipatch,
   const CTYPE_OFFSET offset,
   const BWIDTH boundary_width,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   CTYPE isov_coord[])
  {
    const int DIM3(3);
    const CTYPE ONE_THIRD(1.0/3.0);
    const CTYPE TWO_THIRDS(2.0/3.0);
    
    for (int d = 0; d < DIM3; d++) {
      // Since isov degree is 3, some edge incident on isov
      //   intersects lower facet or upper facet, but not both.
      const bool flag_intersects_lower_facet =
        isodual_table_vinfo.VertexInfo(table_index, ipatch).IntersectsFacet(d);

      const CTYPE c = random_coord_in_range_sampleB<CTYPE>
        (flag_intersects_lower_facet, ONE_THIRD, TWO_THIRDS,
         offset, boundary_width, random_engine, random_pos);
      
      isov_coord[d] = cube_coord[d] + c;
    }
  }


  /*!
   *  @brief Position degree 4 isosurface vertex in random location
   *    in a cube.
   *  - Degree 4 vertices are either positioned on an axis-oriented
   *    square incident on the cube center or on an axis-oriented
   *    line segment, depending on the facets crossed by their
   *    incident edges.
   *  - NOTE: Only for dimension 3.
   *  @param isov Index of isosurface vertex.
   *    @pre Isosurface vertex isov has degree 4.
   *  @param ipatch Patch index of isosurface vertex isov.
   *  @isov_coord[] Coordinates of isosurface vertex isov.
   */
  template <typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_COORD_TYPE,
            typename ISOV_INDEX_TYPE, typename TABLE_INDEX_TYPE,
            typename CTYPE_OFFSET, typename BWIDTH,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_degree4_isov_random_separateB_3D
  (const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const GRID_COORD_TYPE cube_coord[],
   const ISOV_INDEX_TYPE isov,
   const TABLE_INDEX_TYPE table_index,
   const int ipatch,
   const CTYPE_OFFSET offset,
   const BWIDTH boundary_width,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   CTYPE isov_coord[])
  {
    const int DIM3(3);
    const CTYPE ONE_HALF(1.0/2.0);

    for (int d = 0; d < DIM3; d++) {
      const bool flag_intersects_lower_facet =
        isodual_table_vinfo.VertexInfo(table_index, ipatch).IntersectsFacet(d);
      const bool flag_intersects_upper_facet =
        isodual_table_vinfo.VertexInfo(table_index, ipatch).IntersectsFacet(d+DIM3);

      if (flag_intersects_lower_facet && flag_intersects_upper_facet) {
        // Isosurface vertex is halfway between lower and upper facets.
        isov_coord[d] = cube_coord[d] + 0.5;
      }
      else if (flag_intersects_lower_facet || flag_intersects_upper_facet) {
        // Incident edges intersect one facet but not the other.
        const CTYPE c = random_coord_in_range_sampleB<CTYPE>
          (flag_intersects_lower_facet, ONE_HALF, ONE_HALF,
           offset, boundary_width, random_engine, random_pos);
      
        isov_coord[d] = cube_coord[d] + c;
      }
      else {
        const CTYPE rpos = random_pos.RandomAndClamp
         (offset, 1.0-offset, boundary_width, random_engine);
        
        // No incident edge intersects facet orthogonal to d.
        // Put isosurface vertex in random position.
        isov_coord[d] = cube_coord[d] + rpos;
      }
    }
  }


  /*!
   *  @brief Position degree 5 or 6 isosurface vertex in random location
   *    in a cube.
   *  - If cube has two vertices, position degree 5 or 6 vertex
   *    in (2/3)x(2/3)x(2/3) subcube that does not contain degree 3 vertex.
   *  - Otherwise, position vertex anywhere in cube.
   *  - NOTE: Only for dimension 3.
   *  @pre Cube has at most 2 vertices.
   *  @param isov Index of isosurface vertex.
   *    @pre Isosurface vertex isov has degree 4.
   *  @param ipatch Patch index of isosurface vertex isov.
   *    @pre ipatch = 0 or 1.
   *  @isov_coord[] Coordinates of isosurface vertex isov.
   */
  template <typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_COORD_TYPE, typename GRID_CUBE_DATA_TYPE,
            typename CTYPE_OFFSET, typename BWIDTH,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_degree5_or_degree6_isov_random_separateB_3D
  (const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const GRID_COORD_TYPE cube_coord[],
   const GRID_CUBE_DATA_TYPE & active_cube,
   const int ipatch,
   const CTYPE_OFFSET offset,
   const BWIDTH boundary_width,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   CTYPE isov_coord[])
  {
    typedef typename GRID_CUBE_DATA_TYPE::TABLE_INDEX_TYPE TABLE_INDEX_TYPE;
    
    const int DIM3(3);
    const CTYPE ONE_THIRD(1.0/3.0);
    const CTYPE TWO_THIRDS(2.0/3.0);

    if (active_cube.num_isov == 2) {
      const TABLE_INDEX_TYPE table_index = active_cube.table_index;
      const int ipatch_other = 1-ipatch;
      
      for (int d = 0; d < DIM3; d++) {
        const bool flag_intersects_lower_facet =
          isodual_table_vinfo.VertexInfo(table_index, ipatch_other).IntersectsFacet(d);

        // Move isov to region NOT containing ipatch_other.
        // Region has size (2/3)x(2/3)x(2/3).
        const CTYPE c = random_coord_in_range_sampleB<CTYPE>
          (!flag_intersects_lower_facet, TWO_THIRDS, ONE_THIRD,
           offset, boundary_width, random_engine, random_pos);
      
        isov_coord[d] = cube_coord[d] + c;
      }
    }
    else if (active_cube.num_isov == 1) {
      // Put isosurface vertex in random position.
      position_isov_random_offset_sample_boundary
        (DIM3, cube_coord, offset, boundary_width,
         random_engine, random_pos, isov_coord);
    }
    else {
      IJK::PROCEDURE_ERROR
        error("position_degree5_or_degree6_isov_random_3D");
      
      if (active_cube.num_isov == 0) {
        error.AddMessage
          ("Programming error. Active cube contains no isosurface vertices.");
        throw error;
      }
      else {
        error.AddMessage
          ("Programming error. Active cube contains more than 2 isosurface vertices.");
        throw error;        
      }
    }
  }

    
  /*!
   *  @brief Position two dual isosurface vertices in random location in a cube.
   *  - Cube contains exactly two isosurface vertices.
   *  - Position based on degrees.
   *  - Sample region boundaries.
   *  - Degree 3 vertices are placed in (1/3)x(1/3)x(1/3) subcube.
   *  - Degree 4 vertices are placed on an axis-oriented square
   *    incident on the cube center.
   *  - Degree 5 and 6 vertices are placed in (2/3)x(2/3)x(2/3) subcubes.
   *  - NOTE: Only for dimension 3.
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   *  @param isodual_table_vinfo Isodual table vertex information.
   *    @pre compute_dual_cube_isotable_vertex_connectivity()
   *       has been called on isodual_table_vinfo.
   *  @param active_cube Information about grid cube containing
   *     the 2 isosurface vertices.
   *    @pre active_cube.num_isov == 2.
   *  @param offset Position isosurface vertices at distance offset
   *      from cube facets.
   *    @pre 0 <= offset <= 1/6.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_CUBE_DATA_TYPE,
            typename CTYPE_OFFSET, typename BWIDTH,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_isovII_random_sampleB_sep_by_degree_3D
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const GRID_CUBE_DATA_TYPE & active_cube,
   const CTYPE_OFFSET offset,
   const BWIDTH boundary_width,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   CTYPE * coord)  
  {
    typedef typename GRID_CUBE_DATA_TYPE::CUBE_INDEX_TYPE CUBE_INDEX_TYPE;
    typedef typename GRID_CUBE_DATA_TYPE::ISOV_TYPE ISOV_TYPE;
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;

    const int DIM3(3);
    const TABLE_INDEX table_index = active_cube.table_index;
    const CUBE_INDEX_TYPE icube = active_cube.cube_index;
    CTYPE cube_coord[DIM3];

    grid.ComputeCoord(icube, cube_coord);
    
    for (int ipatch = 0; ipatch < active_cube.num_isov; ipatch++) {
      const int degree =
        isodual_table_vinfo.VertexInfo(table_index,ipatch).degree;
      const ISOV_TYPE isov = active_cube.first_isov + ipatch;
      CTYPE * isov_coord = coord + DIM3*isov;

      // *** DEBUG ***
      /*
      using namespace std;
      cerr << "Reposition degree " << degree << " isov " << isov;
      grid.PrintIndexAndCoord(cerr, " in cube ", icube, ".\n");
      const int table_index = active_cube.table_index;
      cerr << "  table_index: " << table_index;
      cerr << "  ambiguous facet: "
           << isodual_table.IndexOfSomeAmbiguousFacet(table_index)
           << endl;
      */
      
      if (degree == 3) {
        position_degree3_isov_random_sampleB_3D
          (isodual_table_vinfo, cube_coord, table_index, ipatch,
           offset, boundary_width, random_engine, random_pos, isov_coord);
      }
      else if (degree == 4) {
        position_degree4_isov_random_separateB_3D
          (isodual_table_vinfo, cube_coord, isov, table_index, ipatch,
           offset, boundary_width, random_engine, random_pos, isov_coord);
      }
      else {
        // degree = 5 or 6 and other isov has degree 3.
        position_degree5_or_degree6_isov_random_separateB_3D
          (isodual_table_vinfo, cube_coord, active_cube, ipatch,
           offset, boundary_width, random_engine, random_pos, isov_coord);
      }
    }
  }

  
  /*!
   *  @brief Position dual isosurface vertices at random location.
   *  - Sample  boundary.
   *  - More than one vertex can be in a cube.
   *  - If cube contains multiple isosurface then vertices are positioned
   *    separated based on degree.
   *  - Degree 3 vertices are placed in (1/3)x(1/3)x(1/3) subcube.
   *  - Degree 4 vertices are placed on an axis-oriented square
   *    incident on the cube center.
   *  - Degree 5 and 6 vertices are placed in (2/3)x(2/3)x(2/3) subcubes.
   *  - NOTE: Only for dimension 3.
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   *  @param isodual_table_vinfo Isodual table vertex information.
   *    @pre compute_dual_cube_isotable_vertex_connectivity()
   *       has been called on isodual_table_vinfo.
   *  @param offset Position isosurface vertices at distance offset
   *      from cube facets.
   *    @pre 0 <= offset <= 1/6.
   *  @param boundary_width Boundary width.
   *  - Generate coordinates in range [rmin-boundary_width,rmax+boundary_width]
   *    and then clamp to [rmin,rmax].   
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_CUBE_DATA_TYPE, typename NTYPEC,
            typename CTYPE_OFFSET, typename CTYPE_BWIDTH,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_all_isov_multi_random_sampleB_sep_by_degree_3D
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const GRID_CUBE_DATA_TYPE active_cube_list[],
   const NTYPEC num_active_cubes,
   const CTYPE_OFFSET offset,
   const CTYPE_BWIDTH boundary_width,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   CTYPE * coord)
  {
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename GRID_CUBE_DATA_TYPE::ISOV_TYPE ISOV_TYPE;
    
    const int DIM3(3);
    const CTYPE ONE_SIXTH(1.0/6.0);
    CTYPE cube_coord[DIM3];

    IJK::PROCEDURE_ERROR error
      ("position_all_isov_multi_random_sampleB_sep_by_degree_3D");

    if (grid.Dimension() != DIM3) {
      error.AddMessage("Programming error. Dimension not equal to 3.");
      error.AddMessage("  grid.Dimension(): ", grid.Dimension());
      throw error;
    }
    
    if (offset < 0.0 || offset > ONE_SIXTH) {
      error.AddMessage("Programming error. Offset not in range [0,1/6].");
      throw error;
    }

    for (NTYPEC i = 0; i < num_active_cubes; i++) {
      const TABLE_INDEX table_index =
        active_cube_list[i].table_index;
      const int num_isov_in_cube =
        isodual_table.NumIsoVertices(table_index);
      const VTYPE icube = active_cube_list[i].cube_index;
      
      grid.ComputeCoord(icube, cube_coord);
      
      if (num_isov_in_cube == 1) {
        const ISOV_TYPE isov =
          active_cube_list[i].first_isov;
        CTYPE * isov_coord = coord + DIM3*isov;

        // Put isosurface vertex in random position.
        // Apply offset.
        position_isov_random_offset_sample_boundary
          (DIM3, cube_coord, offset, boundary_width,
           random_engine, random_pos, isov_coord);
      }
      else if (num_isov_in_cube == 2) {
        position_isovII_random_sampleB_sep_by_degree_3D
          (grid, isodual_table, isodual_table_vinfo, active_cube_list[i],
           offset, boundary_width, random_engine, random_pos, coord);
      }
      else {
        for (int ipatch = 0; ipatch < num_isov_in_cube; ipatch++) {
          const ISOV_TYPE isov =
            active_cube_list[i].first_isov + ipatch;
          CTYPE * isov_coord = coord + DIM3*isov;

          position_degree3_isov_random_sampleB_3D
            (isodual_table_vinfo, cube_coord, table_index, ipatch,
             offset, boundary_width, random_engine, random_pos, isov_coord);
        }
      }
    }
  }
  
  /*!
   *  @overload
   *  @brief Position dual isosurface vertices at random location. (C++ vector.)
   *  - Version using C++ STL vector for arrays active_cube_list[] and coord[].
   *  - More than one vertex can be in a cube.
   *  - If cube contains multiple isosurface then vertices are positioned
   *    separated based on degree.
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   *  @param boundary_width Boundary width.
   *  - Generate coordinates in range [rmin-boundary_width,rmax+boundary_width]
   *    and then clamp to [rmin,rmax].   
   *  @param[out] coord[] Array of isosurface vertex coordinates.
   *    @pre coord[] has size at least (number of isosurface vertices)*3.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
 	    typename ISODUAL_TABLE_VINFO_TYPE,
            typename GRID_CUBE_DATA_TYPE, typename CTYPE_OFFSET,
            typename CTYPE_BWIDTH,
	    typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
	    typename CTYPE>
  void position_all_isov_multi_random_sep_by_degree_3D
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   const CTYPE_OFFSET offset,
   const CTYPE_BWIDTH boundary_width,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   std::vector<CTYPE> & coord)
  {
    if (boundary_width == 0) {
      position_all_isov_multi_random_sep_by_degree_3D
        (grid, isodual_table, isodual_table_vinfo,
         IJK::vector2pointer(active_cube_list), active_cube_list.size(),
         offset, random_engine, random_pos,
         IJK::vector2pointerNC(coord));
    }
    else {
      // Sample boundary.
      position_all_isov_multi_random_sampleB_sep_by_degree_3D
        (grid, isodual_table, isodual_table_vinfo,
         IJK::vector2pointer(active_cube_list), active_cube_list.size(),
         offset, boundary_width, random_engine, random_pos,
         IJK::vector2pointerNC(coord));
    }
  }

  //@}


  // *****************************************************************
  //! @name Position at random locations, allow multiple isov per cube.
  // *****************************************************************

  //@{
  
  /*!
   *  @brief Position dual isosurface vertices at random location.
   *  - More than one vertex can be in a cube.
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename STYPE,
            typename GRID_CUBE_DATA_TYPE, typename DUAL_ISOV_TYPE,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_all_isov_multi_random
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const STYPE isovalue,
   const std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const RANDOM_POS_PARAM & random_pos_param,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   std::vector<CTYPE> & coord)
  {
    const int DIM2(2);
    const int DIM3(3);
    IJK::PROCEDURE_ERROR error("position_all_isov_multi_random");

    const int dimension = scalar_grid.Dimension();
    coord.resize(isov_list.size()*dimension);

    COORD_TYPE position_offset = 0.0;
    if (random_pos_param.flag_gen_isov_apply_offset)
      { position_offset = random_pos_param.position_offset; }

    // Initialize boundary_width to 0.0.
    CTYPE boundary_width = 0.0;
    if (random_pos_param.flag_gen_isov_on_region_boundary) {
      boundary_width = random_pos_param.boundary_width;
    }
      
    if (random_pos_param.random_pos_isov_separation_method ==
        RANDOM_POS_SEPARATE_USING_CENTROID_POS) {

      /* OBSOLETE?
      const bool flag_use_centroid_pos_on_multi_isov =
        random_pos_param.flag_use_centroid_pos_on_multi_isov;
      const bool flag_use_centroid_pos_on_doubly_connected_isov =
        random_pos_param.flag_use_centroid_pos_on_doubly_connected_isov;
      */

      position_all_isov_multi_random_use_centroid_pos
        (scalar_grid, isodual_table, isovalue,
         isov_list, position_offset, boundary_width,
         random_engine, random_pos,
         coord);
    }
    else if (random_pos_param.random_pos_isov_separation_method ==
             RANDOM_POS_SEPARATE_BY_CUBE_CENTER) {
      position_all_isov_multi_random_sep_by_cc
        (scalar_grid, isodual_table, isodual_table_vinfo,
         active_cube_list, isov_list,
         position_offset, random_engine, random_pos, coord);
    }
    else if (random_pos_param.random_pos_isov_separation_method ==
             RANDOM_POS_SEPARATE_BASED_ON_DEGREE) {

      if (dimension == DIM2) {
        // In 2D, separation based on degree is same
        //   as separation by cube center.
        position_all_isov_multi_random_sep_by_cc
          (scalar_grid, isodual_table, isodual_table_vinfo,
           active_cube_list, isov_list,
           position_offset, random_engine, random_pos, coord);
      }
      else if (dimension == DIM3) {
        position_all_isov_multi_random_sep_by_degree_3D
          (scalar_grid, isodual_table, isodual_table_vinfo,
           active_cube_list, position_offset, boundary_width,
           random_engine, random_pos, coord);
      }
      else {
        error.AddMessage("Programming error. Dimension greater than 3.");
        error.AddMessage
          ("Separation by degree is not implemented for grid dimension > 3.");
        throw error;
      }
    }
    else {
      if (random_pos_param.flag_gen_isov_on_region_boundary) {
        position_all_isov_random_sample_boundary
          (scalar_grid, isov_list, random_pos_param.boundary_width,
           random_engine, random_pos, coord);
      }
      else {
        position_all_isov_random
          (scalar_grid, isov_list, position_offset,
           random_engine, random_pos, coord);
      }
    }
  }
  

  /*!
   *  @overload
   *  @brief Position dual isosurface vertices at random location.
   *  - More than one vertex can be in a cube.
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename STYPE,
            typename GRID_CUBE_DATA_TYPE, typename DUAL_ISOV_TYPE,
            typename CTYPE>
  void position_all_isov_multi_random
  (const GRID_TYPE & scalar_grid,   
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const STYPE isovalue,
   const std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const RANDOM_POS_PARAM & random_pos_param,
   std::vector<CTYPE> & coord)
  {
    const RANDOM_DISTRIBUTION random_pos_distribution =
      random_pos_param.distribution;
    const RANDOM_SEED_TYPE random_seed = random_pos_param.random_seed;
    IJK::PROCEDURE_ERROR
      error("position_all_isov_multi_random");

    std::minstd_rand random_engine(random_seed);

    if (random_pos_distribution == UNIFORM_DISTRIBUTION) {
      IJK::GENERATE_UNIFORM_RANDOM random_pos;
      
      position_all_isov_multi_random
        (scalar_grid, isodual_table, isodual_table_vinfo,
         isovalue, active_cube_list, isov_list, random_pos_param,
         random_engine, random_pos, coord);
    }
    else if (random_pos_distribution == U_QUADRATIC_DISTRIBUTION) {
      IJK::GENERATE_U_QUADRATIC_RANDOM random_pos;

      position_all_isov_multi_random
        (scalar_grid, isodual_table, isodual_table_vinfo,
         isovalue, active_cube_list, isov_list, random_pos_param,
         random_engine, random_pos, coord);
    }
    else {
      error.AddMessage
        ("Programming error. Random position distribution not set.");
      throw error;
    }

  }


  //@}

#endif

}

#endif
