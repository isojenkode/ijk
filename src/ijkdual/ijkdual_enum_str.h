/*!
 *  @file ijkdual_enum_str.h
 *  @brief Strings corresponding to ijkdual enumerated types.
 *  - Version 0.6.0
 */


/*
  IJK: Isosurface Jeneration Kode
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

#ifndef _IJKDUAL_ENUM_STR_H_
#define _IJKDUAL_ENUM_STR_H_

#include <vector>

#include "ijkenum.tpp"

#include "ijkdual_types.h"
#include "ijktri_interior_vertex.tpp"

namespace IJKDUAL {

  // **************************************************
  // ENUM VALUE-STRING PAIRS
  // **************************************************

  /// @brief Enum-string pairs for enum type VERTEX_POSITION_METHOD.
  static std::vector<IJK::ENUM_STR<VERTEX_POSITION_METHOD> >
  vertex_position_method_strings =
    { {CUBE_CENTER, "CubeCenter"},
      {CENTROID_EDGE_ISO, "CentroidEdge"},
      {DIAGONAL_INTERPOLATION, "DiagonalInterpolation"},
      {IVOL_LIFTED02, "IntervalVolumeLifted"},
      {QDUAL_INTERPOLATION, "QualityDualInterpolation"},
      {RANDOM_ISOV_POS, "Random"},
      {UNDEFINED_VERTEX_POSITION_METHOD, "Undefined"} };

  /// @brief Enum-string pairs for enum type
  ///   ISO_GRID_EDGE_INTERSECTION_OFFSET_METHOD.
  static std::vector<IJK::ENUM_STR<ISO_GRID_EDGE_INTERSECTION_OFFSET_METHOD> >
  iso_grid_edge_intersection_offset_method_strings =
    { {NO_OFFSET, "NoOffset"},
      {OFFSET_ONLY_MULTI_IN_CUBE, "OffsetOnlyMultiInCube"},
      {OFFSET_ALL, "OffsetAll"},
      {UNDEFINED_ISO_GRID_EDGE_INTERSECTION_OFFSET_METHOD, "Undefined"} };
  
  /// @brief Enum-string pairs for enum type QUAD_TRI_METHOD.
  static std::vector<IJK::ENUM_STR<QUAD_TRI_METHOD> >
  quad_tri_method_strings =
    { {UNDEFINED_TRI, "Undefined"},
      {UNIFORM_TRI, "Uniform"},
      {SPLIT_MAX_ANGLE, "SplitMaxAngle"},
      {MAX_MIN_ANGLE, "MaxMinAngle"},
      {TRI4_ALL_QUADS, "Tri4AllQuads"},
      {TRI4_MAX_MIN_ANGLE, "MaxMinAngle"},
      {TRI_DIAGONALS_OUTSIDE_ENVELOPE, "TriDiagonalsOutsideEnvelope"} };

  /// @brief Enum-string pairs for enum type ENVELOPE_QUAD_TRI_METHOD.
  static std::vector<IJK::ENUM_STR<ENVELOPE_QUAD_TRI_METHOD> >
  envelope_quad_tri_method_strings =
    { {ENVELOPE_UNDEFINED_TRI, "Undefined"},
      {ENVELOPE_ONLY_TRI2, "OnlyTwoTriangles"},
      {ENVELOPE_ONLY_TRI4, "OnlyFourTriangles"},
      {ENVELOPE_PREFER_TRI2, "PreferTwoTriangles"},
      {ENVELOPE_MAX_MIN_ANGLE, "MaxMinAngle"} };

  /// @brief Enum-string pairs for enum type ISOV_SEPARATION_METHOD.
  static std::vector<IJK::ENUM_STR<ISOV_SEPARATION_METHOD> >
  isov_separation_method_strings =
    { {SEPARATE_BY_CUBE_CENTER, "ByCubeCenter"},
      {SEPARATE_BY_PLANES, "ByAxisParallelPlanes"},
      {UNDEFINED_ISOV_SEPARATION_METHOD, "Undefined" } };

  /// @brief Enum-string pairs for enum type
  ///   DOUBLY_CONNECTED_ISOV_REPOSITION_METHOD.
  static std::vector<IJK::ENUM_STR<DOUBLY_CONNECTED_ISOV_REPOSITION_METHOD> >
  doubly_connected_isov_reposition_method_strings =
    { {DC_POS_CLAMPIII_ONE_THIRD, "ClampIIIOneThird"},
      {DC_POS_CLAMP_ALL, "ClampAll"},
      {DC_POS_CENTERIII, "CenterIII"},
      {DC_POS_CENTER_ALL, "CenterAll"},
      {UNDEFINED_DC_ISOV_REPOSITION_METHOD, "Undefined"} };
  
  /// @brief Enum-string pairs for enum type POLY_EDGE_INTERSECTION_METHOD.
  /// - IJK::POLY_EDGE_INTERSECTION_METHOD is defined in file
  ///   ijktri_interior_vertex.tpp.
  static std::vector<IJK::ENUM_STR<IJK::POLY_EDGE_INTERSECTION_METHOD> >
  quad_edge_intersection_method_strings =
    { {IJK::POLY_EDGE_MULTILINEAR_INTERPOLATION, "BilinearSurface"},
      {IJK::INTERPOLATE_EDGE_ENDPOINT_SCALARS, "InterpolateScalar"},
      {IJK::POLY_EDGE_AVERAGE_PROJECTION, "AverageProjection"},
      {IJK::POLY_EDGE_WEIGHTED_AVERAGE_PROJECTION,
       "WeightedAverageProjection"},
      {IJK::UNDEFINED_POLY_EDGE_INTERSECTION_METHOD, "Undefined"} };

  /// @brief Enum-string pairs for enum type INTERIOR_VERTEX_POSITION_METHOD.
  /// - IJK::INTERIOR_VERTEX_POSITION_METHOD is defined
  ///   in file ijktri_interior_vertex.tpp.
  static std::vector<IJK::ENUM_STR<IJK::INTERIOR_VERTEX_POSITION_METHOD> >
  interior_vertex_position_method_strings =
    { {IJK::INTERIOR_VERTEX_ON_GRID_EDGE, "OnGridEdge"},
      {IJK::INTERIOR_VERTEX_IN_ENVELOPE, "InEnvelope"},
      {IJK::INTERIOR_VERTEX_AT_CENTROID, "AtCentroid"},
      {IJK::UNDEFINED_INTERIOR_VERTEX_POSITION_METHOD, "Undefined" } };

  /// @brief Enum-string pairs for enum type RANDOM_DISTRIBUTION.
  static std::vector<IJK::ENUM_STR<RANDOM_DISTRIBUTION> >
  random_distribution_strings =
    { {UNIFORM_DISTRIBUTION, "Uniform"},
      {U_QUADRATIC_DISTRIBUTION, "U-Quadratic"},
      {UNDEFINED_DISTRIBUTION, "Undefined"} };

  /// @brief Enum-string pairs for enum type RANDOM_POS_ISOV_SEPARATION
  static std::vector<IJK::ENUM_STR<RANDOM_POS_ISOV_SEPARATION_METHOD> >
  random_pos_isov_separation_method_strings =
    { {RANDOM_POS_NO_SEPARATION, "NoSeparation" },
      {RANDOM_POS_SEPARATE_USING_CENTROID_POS, "CentroidPos" },
      {RANDOM_POS_SEPARATE_BY_CUBE_CENTER, "ByCubeCenter" },
      {RANDOM_POS_SEPARATE_BASED_ON_DEGREE, "BasedOnDegree" },
      {UNDEFINED_RANDOM_POS_ISOV_SEPARATION_METHOD, "Undefined"} };


  // **************************************************
  // CLASS IJKDUAL_ENUM_STRINGS
  // **************************************************

  /// @brief Class for converting between enum types and strings.
  /// - Do not confuse with COMMAND_LINE_ENUM_STRINGS.
  class IJKDUAL_ENUM_STRINGS {

  protected:

    /// @brief Strings corresponding to enum VERTEX_POSITION_METHOD.
    /// - Set vertex_position_method_list
    ///   from vertex_position_method_strings in constructor.
    IJK::ENUM_LIST<VERTEX_POSITION_METHOD> vertex_position_method_list;

    /// @brief Strings corresponding to enum
    ///   ISO_GRID_EDGE_INTERSECTION_OFFSET_METHOD.
    /// - Set iso_grid_edge_intersection_offset_method_list
    ///   from iso_grid_edge_intersection_offset_method_strings
    ///   in constructor.
    IJK::ENUM_LIST<ISO_GRID_EDGE_INTERSECTION_OFFSET_METHOD>
    iso_grid_edge_intersection_offset_method_list;
    
    /// @brief Strings corresponding to enum QUAD_TRI_METHOD.
    /// - Set quad_tri_method_list from quad_tri_method_strings in constructor.
    IJK::ENUM_LIST<QUAD_TRI_METHOD> quad_tri_method_list;

    /// @brief Strings corresponding to enum ENVELOPE_QUAD_TRI_METHOD.
    /// - Set envelope_quad_tri_method_list
    ///   from envelope_quad_tri_method_strings in constructor.
    IJK::ENUM_LIST<ENVELOPE_QUAD_TRI_METHOD> envelope_quad_tri_method_list;
    
    /// @brief Strings corresponding to enum ISOV_SEPARATION_METHOD.
    /// - Set isov_separation_method_list
    ///   from isov_separation_method_strings in constructor.
    IJK::ENUM_LIST<ISOV_SEPARATION_METHOD> isov_separation_method_list;

    /// @brief Strings corresponding to enum
    ///   DOUBLY_CONNECTED_ISOV_REPOSITION_METHOD.
    /// - Set doubly_connected_isov_reposition_method_list
    ///   from doubly_connected_isov_reposition_method_strings
    ///   in constructor.
    IJK::ENUM_LIST<DOUBLY_CONNECTED_ISOV_REPOSITION_METHOD>
    doubly_connected_isov_reposition_method_list;

    /// @brief Strings describing quad edge intersection method.
    /// - Set quad_edge_intersection_method_list
    ///   from quad_edge_intersection_strings in constructor.
    IJK::ENUM_LIST<IJK::POLY_EDGE_INTERSECTION_METHOD>
    quad_edge_intersection_method_list;
    
    /// @brief Strings corresponding to enum INTERIOR_VERTEX_POSITION_METHOD.
    /// - Set interior_vertex_position_method_list
    ///   from interior_vertex_position_method_strings in constructor.
    IJK::ENUM_LIST<IJK::INTERIOR_VERTEX_POSITION_METHOD>
    interior_vertex_position_method_list;

    /// @brief Strings corresponding to enum RANDOM_DISTRIBUTION.
    /// - Set random_distribution_list
    ///   from random_distribution_strings in constructor.
    IJK::ENUM_LIST<RANDOM_DISTRIBUTION> random_distribution_list;

    /// @brief Strings corresponding to enum
    ///   RANDOM_POS_ISOV_SEPARATION_METHOD.
    /// - Set random_pos_isov_separation_method_list
    ///   from random_pos_isov_separation_method_strings in constructor.
    IJK::ENUM_LIST<RANDOM_POS_ISOV_SEPARATION_METHOD>
    random_pos_isov_separation_method_list;

    
  public:
    
    /// @brief Constructor
    IJKDUAL_ENUM_STRINGS():
      vertex_position_method_list(UNDEFINED_VERTEX_POSITION_METHOD,
                                  vertex_position_method_strings),
      iso_grid_edge_intersection_offset_method_list
      (UNDEFINED_ISO_GRID_EDGE_INTERSECTION_OFFSET_METHOD,
       iso_grid_edge_intersection_offset_method_strings),
      quad_tri_method_list(UNDEFINED_TRI, quad_tri_method_strings),
      envelope_quad_tri_method_list(ENVELOPE_UNDEFINED_TRI,
                                    envelope_quad_tri_method_strings),
      isov_separation_method_list(UNDEFINED_ISOV_SEPARATION_METHOD,
                                  isov_separation_method_strings),
      doubly_connected_isov_reposition_method_list
      (UNDEFINED_DC_ISOV_REPOSITION_METHOD,
       doubly_connected_isov_reposition_method_strings),
      quad_edge_intersection_method_list
      (IJK::UNDEFINED_POLY_EDGE_INTERSECTION_METHOD,
       quad_edge_intersection_method_strings),
      interior_vertex_position_method_list
      (IJK::UNDEFINED_INTERIOR_VERTEX_POSITION_METHOD,
       interior_vertex_position_method_strings),
      random_distribution_list(UNDEFINED_DISTRIBUTION,
                               random_distribution_strings),
      random_pos_isov_separation_method_list
      (UNDEFINED_RANDOM_POS_ISOV_SEPARATION_METHOD,
       random_pos_isov_separation_method_strings)
    {};

    /// @brief Return string of vertex position method.
    std::string VertexPositionMethodString
    (const VERTEX_POSITION_METHOD _vpos_method) const
    { return vertex_position_method_list.String(_vpos_method); }

    /// @brief Return iso_grid_edgeI_offset_method_list.
    std::string IsoGridEdgeIntersectionOffsetMethodString
    (const ISO_GRID_EDGE_INTERSECTION_OFFSET_METHOD _offset_method) const
    { return iso_grid_edge_intersection_offset_method_list.String
        (_offset_method); }
    
    /// @brief Return string of quadrilateral triangulation method.
    std::string QuadTriMethodString
    (const QUAD_TRI_METHOD _quad_tri_method) const
    { return quad_tri_method_list.String(_quad_tri_method); }

    /// @brief Return string of envelope quadrilateral triangulation method.
    std::string EnvelopeQuadTriMethodString
    (const ENVELOPE_QUAD_TRI_METHOD _envelope_quad_tri_method) const
    {
      return envelope_quad_tri_method_list.String
        (_envelope_quad_tri_method);
    }

    /// @brief Return string of quad edge intersection method.
    std::string QuadEdgeIntersectionMethodString
    (const IJK::POLY_EDGE_INTERSECTION_METHOD _qei_method) const
    { return quad_edge_intersection_method_list.String(_qei_method); }
    
    /// @brief Return string of isosurface vertices separation method.
    std::string IsovSeparationMethodString
    (const ISOV_SEPARATION_METHOD _isov_separation_method) const
    { return isov_separation_method_list.String(_isov_separation_method); }

    /// @brief Return string of doubly connected isov reposition method.
    std::string DoublyConnectedIsovRepositionMethodString
    (const DOUBLY_CONNECTED_ISOV_REPOSITION_METHOD _reposition_method) const
    {
      return doubly_connected_isov_reposition_method_list.String
        (_reposition_method);
    }
    
    /// @brief Return string of interior vertex position method.
    std::string InteriorVertexPositionMethodString
    (const IJK::INTERIOR_VERTEX_POSITION_METHOD _position_method) const
    { return
        interior_vertex_position_method_list.String(_position_method); }

    /// @brief Return string of random distribution.
    std::string RandomDistributionString
    (const RANDOM_DISTRIBUTION _random_distribution) const
    { return random_distribution_list.String(_random_distribution); }
    
    /// @brief Return string of random position isosurface vertex
    ///   separation method.
    std::string RandomPosIsovSeparationMethodString
    (const RANDOM_POS_ISOV_SEPARATION_METHOD _separation_method) const
    { return random_pos_isov_separation_method_list.String
        (_separation_method); }
  };
  
}

#endif
