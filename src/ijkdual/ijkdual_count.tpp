/*!
 *  @file ijkdual_count.tpp
 *  @brief Templates for counting ambiguous cubes and facets.
 *  - Version 0.6.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2025-2025 Rephael Wenger

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

#ifndef _IJKDUAL_COUNT_
#define _IJKDUAL_COUNT_

/// @brief Count number of cubes with ambiguous facets.
/// - Note: All such cubes are active.
template <typename ISODUAL_TABLE_TYPE,
          typename GRID_CUBE_DATA_TYPE>
unsigned int count_num_cubes_with_ambig_facets
(const ISODUAL_TABLE_TYPE & isodual_table,
 const std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list)
{
  typedef typename std::vector<GRID_CUBE_DATA_TYPE>::size_type
    SIZE_TYPE;
  typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX_TYPE;

  SIZE_TYPE num_cubes_with_ambig_facets = 0;
  
  for (SIZE_TYPE i = 0; i < active_cube_list.size(); i++) {
    const TABLE_INDEX_TYPE table_index =
      active_cube_list[i].table_index;

    if (isodual_table.NumAmbiguousFacets(table_index) > 0)
      { ++num_cubes_with_ambig_facets; }
  }

  return num_cubes_with_ambig_facets;
}


/// @brief Count number of ambiguous cubes and
///    number of ambiguous grid facets.
template <typename GRID_TYPE,
          typename ISODUAL_TABLE_TYPE,
          typename GRID_CUBE_DATA_TYPE,
          typename NTYPE1, typename NTYPE2>
void count_num_ambig
(const GRID_TYPE & grid,
 const ISODUAL_TABLE_TYPE & isodual_table,
 const std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
 NTYPE1 & num_cubes_with_ambig_facets,
 NTYPE2 & num_ambig_grid_facets)
{
  typedef typename std::vector<GRID_CUBE_DATA_TYPE>::size_type
    SIZE_TYPE;
  typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
  typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX_TYPE;

  using IJKDUAL::BOUNDARY_BITS_TYPE;
  
  const int dimension = grid.Dimension();
  NTYPE2 num_ambig_interior_grid_facets = 0;
  NTYPE2 num_ambig_boundary_grid_facets = 0;

  // Initialize
  num_cubes_with_ambig_facets = 0;
  num_ambig_grid_facets = 0;
  
  for (SIZE_TYPE i = 0; i < active_cube_list.size(); i++) {
    const TABLE_INDEX_TYPE table_index =
      active_cube_list[i].table_index;

    if (isodual_table.NumAmbiguousFacets(table_index) > 0) {

      ++num_cubes_with_ambig_facets;
      
      BOUNDARY_BITS_TYPE boundary_bits;
      const VERTEX_INDEX_TYPE icube = active_cube_list[i].cube_index;
      grid.ComputeBoundaryCubeBits(icube, boundary_bits);

      if (boundary_bits == 0) {
        // Interior cube.
        num_ambig_interior_grid_facets +=
          isodual_table.NumAmbiguousFacets(table_index);
      }
      else {
        // Cube on the boundary.
        BOUNDARY_BITS_TYPE mask = BOUNDARY_BITS_TYPE(1);
        for (int jf = 0; jf < grid.NumCubeFacets(); jf++) {

          // Unfortunately, facets are indexed diffently
          //   in grid boundary_bits than in ijkdualtable.
          const int orth_dir = jf/2;
          const int iside = jf%2;
          const int isotable_facet_index =
            iside*dimension + orth_dir;

          if (isodual_table.IsFacetAmbiguous
              (table_index, isotable_facet_index)) {
          
            if ((boundary_bits | mask) == 1) {
              // facet jf is an interior facet.
              ++num_ambig_interior_grid_facets;
            }
            else {
              // facet jf is a boundary facet.
              ++num_ambig_boundary_grid_facets;
            }            
          }
        }
      }
    }
  }

  // Previous loop double counts interior grid facets.
  num_ambig_interior_grid_facets = num_ambig_interior_grid_facets/2;
  
  num_ambig_grid_facets =
    num_ambig_interior_grid_facets + num_ambig_boundary_grid_facets;
}

#endif
