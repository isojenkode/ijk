/*!
 *  @file ijkdualIO.tpp
 *  @brief IO templates for ijkdual.
 * - Version 0.6.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2017-2025 Rephael Wenger

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

#ifndef _IJKDUALIO_TPP_
#define _IJKDUALIO_TPP_

#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>

#include "ijk.tpp"
#include "ijkisopoly.tpp"
#include "ijkprint.tpp"

#include "ijkdual_types.h"


namespace IJKDUAL {

  // ******************************************************************
  //! @name OUTPUT ISOSURFACE
  // ******************************************************************

  //@{

  /*!
   *  @brief Output dual isosurface.
   *  - Use \a output_info.mesh_type to determine type of isosurface.
   */
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE, 
            typename DUAL_ISOSURFACE_TYPE, typename DUALISO_INFO_TYPE,
            typename IO_TIME_TYPE>
  void output_dual_isosurface
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const DUAL_ISOSURFACE_TYPE & dual_isosurface,
   const DUALISO_INFO_TYPE & dualiso_info, IO_TIME_TYPE & io_time)
  {
    const int THREE(3);
    const int dimension = output_info.Dimension();
    
    if ((dimension == THREE) &&
        (output_info.mesh_type == SIMPLICIAL_COMPLEX)) {
      output_dual_tri_isosurface
        (output_info, dualiso_data, dual_isosurface.vertex_coord, 
         dual_isosurface.simplex_vert, dualiso_info, io_time);
    }
    else if ((dimension == THREE) &&
             (output_info.mesh_type == MIXED_MESH)) {
      output_dual_quad_tri_isosurface
        (output_info, dualiso_data, dual_isosurface.vertex_coord, 
         dual_isosurface.isopoly_vert, dual_isosurface.simplex_vert, 
         dualiso_info, io_time);
    }
    else if ((dimension == THREE) && output_info.flag_dual_collapse) {
      output_dual_quad_tri_isosurface
        (output_info, dualiso_data, dual_isosurface.vertex_coord, 
         dual_isosurface.isopoly_vert, dual_isosurface.simplex_vert,
         dualiso_info, io_time);
    }
    else {
      output_dual_isosurface
        (output_info, dualiso_data, dual_isosurface.vertex_coord, 
         dual_isosurface.isopoly_vert, dualiso_info, io_time);
    }
  }


  /*!
   *  @brief Output dual isosurface and report isosurface information.
   *  - Output file(s) and format(s) depend on flags in \a output_info.
   *  - May write to more than one file.
   *  - Output (report) isosurface information to stdout.
   */
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE, 
            typename DUALISO_INFO_TYPE, typename IO_TIME_TYPE>
  void output_dual_isosurface
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const std::vector<COORD_TYPE> & vertex_coord, 
   const std::vector<VERTEX_INDEX> & slist,
   const DUALISO_INFO_TYPE & dualiso_info, 
   IO_TIME_TYPE & io_time)
  {
    if (!output_info.flag_use_stdout && !output_info.flag_silent) {
      report_iso_info(output_info, dualiso_data, 
                      vertex_coord, slist, dualiso_info);
    }

    if (!output_info.flag_nowrite) 
      { write_dual_mesh(output_info, vertex_coord, slist, io_time); }
  }


  /*!
   *  @brief Output isosurface of  cubes (line segments, quads, hexahedra, ...)
   *  - Isosurface contains only cubes.
   */
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE,
            typename DUAL_ISOSURFACE_TYPE, typename DUALISO_INFO_TYPE,
            typename IO_TIME_TYPE>
  void output_dual_cube_complex_isosurface
  (const OUTPUT_INFO_TYPE & output_info,
   const DUALISO_DATA_TYPE & dualiso_data,
   const DUAL_ISOSURFACE_TYPE & dual_isosurface,
   const DUALISO_INFO_TYPE & dualiso_info, IO_TIME_TYPE & io_time)
  {
    output_dual_isosurface
      (output_info, dualiso_data, dual_isosurface.vertex_coord,
       dual_isosurface.isopoly_vert, dualiso_info, io_time);
  }


  /// Output isosurface of triangles.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE, 
            typename DUALISO_INFO_TYPE, typename IO_TIME_TYPE>
  void output_dual_tri_isosurface
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & tri_vert,
   const DUALISO_INFO_TYPE & dualiso_info, IO_TIME_TYPE & io_time)
  {
    if (!output_info.flag_use_stdout && !output_info.flag_silent) {
      report_iso_info(output_info, dualiso_data,
                      vertex_coord, tri_vert, dualiso_info);
    }

    if (!output_info.flag_nowrite) {
      write_dual_mesh(output_info, vertex_coord, tri_vert, io_time);
    }
  }


  /// Output isosurface of quadrilaterals and triangles.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE, 
            typename DUALISO_INFO_TYPE, typename IO_TIME_TYPE>
  void output_dual_quad_tri_isosurface
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const std::vector<COORD_TYPE> & vertex_coord,
   const std::vector<VERTEX_INDEX> & quad_vert,
   const std::vector<VERTEX_INDEX> & tri_vert,
   const DUALISO_INFO_TYPE & dualiso_info, IO_TIME_TYPE & io_time)
  {
    if (!output_info.flag_use_stdout && !output_info.flag_silent) {
      report_quad_tri_iso_info
        (output_info, dualiso_data, vertex_coord, quad_vert, tri_vert, 
         dualiso_info);
    }

    if (!output_info.flag_nowrite) {
      write_dual_quad_tri_mesh
        (output_info, vertex_coord, quad_vert, tri_vert, io_time);
    }
  }


  /// Output isosurface with colored facets.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE, 
            typename DUALISO_INFO_TYPE, typename COLOR_TYPE,
            typename IO_TIME_TYPE>
  void output_dual_isosurface_color
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const std::vector<COORD_TYPE> & vertex_coord, 
   const std::vector<VERTEX_INDEX> & slist,
   const COLOR_TYPE * front_color, 
   const COLOR_TYPE * back_color,
   const DUALISO_INFO_TYPE & dualiso_info, IO_TIME_TYPE & io_time)
  {
    if (!output_info.flag_use_stdout && !output_info.flag_silent) {
      report_iso_info
        (output_info, dualiso_data, vertex_coord, slist, dualiso_info);
    }
  
    if (!output_info.flag_nowrite) {
      write_dual_mesh_color
        (output_info, vertex_coord, slist, front_color, back_color, io_time);
    }
  }


  /// Output isosurface with colored facets.
  /// - Version with input DUAL_ISOSURFACE_TYPE.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE, 
            typename DUAL_ISOSURFACE_TYPE, typename DUALISO_INFO_TYPE,
            typename COLOR_TYPE, typename IO_TIME_TYPE>
  void output_dual_isosurface_color
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const DUAL_ISOSURFACE_TYPE & dual_isosurface,
   const COLOR_TYPE * front_color, const COLOR_TYPE * back_color,
   const DUALISO_INFO_TYPE & dualiso_info, IO_TIME_TYPE & io_time)
  {
    output_dual_isosurface_color
      (output_info, dualiso_data, 
       dual_isosurface.vertex_coord, dual_isosurface.isopoly_vert,
       front_color, back_color, dualiso_info, io_time);
  }

  //@}


  // ******************************************************************
  //! @name REPORT SCALAR FIELD OR ISOSURFACE INFORMATION
  // ******************************************************************

  //@{

  /// Report information about cubes with multiple vertices.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_INFO_TYPE>  
  void report_isov_info
  (const OUTPUT_INFO_TYPE & output_info, const DUALISO_INFO_TYPE & dualiso_info)
  {
    const int DIM3(3);
    const char * indent4 = "    ";

    using namespace std;

    if (output_info.Dimension() == DIM3) {
      if (output_info.FlagRepositionDoublyConnectedIsosurfaceVertices()) {
        cout << indent4
             << "# repositioned doubly connected isosurface vertices: "
             << dualiso_info.isov.num_repositioned_doubly_connected_isov
             << endl;
      }
    }
    
    
    if (output_info.Dimension() >= DIM3) {
      if (output_info.FlagSeparateIsovNearSharedGridCubeFacets() ||
          output_info.FlagMoveAllIsovAwayFromCubeBoundaries()) {
        cout << indent4
             << "# isosurface vertices moved away from cube facets: "
             << dualiso_info.isov.num_isov_moved_away_from_grid_cube_facets
             << endl;
      }

      if (output_info.FlagSeparateIsovNearSharedGridCubeRidges() ||
          output_info.FlagMoveAllIsovAwayFromCubeRidges()) {
        cout << indent4 << "# isosurface vertices moved away from cube";
        if (output_info.Dimension() == DIM3) { cout << " edges: "; }
        else { cout << " ridges: "; }
        cout << dualiso_info.isov.num_isov_moved_away_from_grid_cube_ridges
             << endl;
      }

      if (output_info.FlagSeparateIsovNearSharedGridVertices()) {
        cout << indent4
             << "# isosurface vertices moved away from grid vertices: ";
        cout << dualiso_info.isov.num_isov_moved_away_from_grid_vertices
             << endl;
      }
      
      if (output_info.FlagSeparateIsoPolyNearSharedGridVertices()) {
        cout << indent4 << "# isosurface vertices moved to separate iso";
        if (output_info.Dimension() == DIM3) { cout << " quads: "; }
        else { cout << " poly: "; }
        cout << dualiso_info.isov.num_isov_moved_to_separate_iso_poly
             << endl;
      }      
    }

  }

  
  /// Report information about cubes with multiple vertices.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_INFO_TYPE>  
  void report_multi_isov_info
  (const OUTPUT_INFO_TYPE & output_info, const DUALISO_INFO_TYPE & dualiso_info)
  {
    const char * indent4 = "    ";
    const char * indent6 = "      ";

    using namespace std;

    if (output_info.allow_multiple_iso_vertices) {
      cout << indent4 << "# active (non-empty) cubes: "
           << dualiso_info.scalar.num_non_empty_cubes << endl;

      if (output_info.flag_report_info) {
        if (dualiso_info.scalar.num_cubes_with_ambig_facets.IsSet()) {
          cout << indent6 << "# active cubes with ambiguous facets: "
               << dualiso_info.scalar.num_cubes_with_ambig_facets.Value()
               << endl;
        }

        if (dualiso_info.scalar.num_ambig_grid_facets.IsSet()) {
          cout << indent4 << "# ambiguous grid facets: "
               << dualiso_info.scalar.num_ambig_grid_facets.Value()
               << endl;          
        }
      }
      
      cout << indent4 << "# cubes with single isosurface vertex: "
           << dualiso_info.multi_isov.num_cubes_single_isov << endl;
      cout << indent4 << "# cubes with multiple isosurface vertices: "
           << dualiso_info.multi_isov.num_cubes_multi_isov << endl;

      if (output_info.flag_split_non_manifold) {
        cout << indent4 << "# cubes changed to 2 iso vertices to avoid non-manifold edges: "
             << dualiso_info.multi_isov.num_non_manifold_split << endl;
        cout << indent4 
             << "# ambiguous ridge cubes changed to avoid non-manifold edges: "
             << dualiso_info.multi_isov.num_ambig_ridge_cubes_changed << endl;
        cout << indent4 
             << "# non-ambiguous ridge cubes changed to avoid non-manifold edges: "
             << dualiso_info.multi_isov.num_non_ambig_ridge_cubes_changed << endl;
      }

      if (output_info.flag_select_split) {
        cout << indent4 << "# cubes changed in selecting isosurface patch splits: "
             << dualiso_info.multi_isov.num_1_2_changed << endl;
      }

      if (output_info.flag_connect_ambiguous) {
        cout << indent4 << "# ambiguous cubes changed to connect isosurface: "
             << dualiso_info.multi_isov.num_connect_changed << endl;
      }
    }

  }


  /// Report triangulation information.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_INFO_TYPE>  
  void report_triangulation_info
  (const OUTPUT_INFO_TYPE & output_info, const DUALISO_INFO_TYPE & dualiso_info)
  {
    const char * indent4 = "    ";
    const char * indent6 = "      ";
    const char * indent8 = "        ";

    using namespace std;

    if (output_info.flag_check_envelope) {
      cout << indent4 << "# quads with some diagonal outside envelope: "
           << dualiso_info.triangulation.num_iso_cubes_with_diag_outside_envelope
           << endl;
    }
    
    if (output_info.flag_trimesh ||
        (output_info.flag_check_envelope &&
         (output_info.flag_allow_tri2_envelope ||
          output_info.flag_allow_tri4_envelope))) {

      cout << indent4 << "# triangulated quads: "
           << dualiso_info.triangulation.num_iso_cubes_tri_total
           << endl;
      
      if (output_info.quad_tri_method != TRI4_ALL_QUADS ||
          (output_info.flag_check_envelope &&
           output_info.flag_allow_tri2_envelope)) {
        
        cout << indent6 << "# quads split into 2 triangles: "
           << dualiso_info.triangulation.num_iso_cubes_tri_no_add
           << endl;

        if (output_info.flag_check_envelope &&
            output_info.flag_allow_tri2_envelope) {

          cout << indent8 << "# quads split into 2 triangles bcuz diagonal outside envelope: "
               << dualiso_info.triangulation.num_tri_no_add_with_diag_outside_envelope
               << endl;
        }
      }

      if (output_info.flag_tri4_quad ||
          (output_info.flag_check_envelope &&
           output_info.flag_allow_tri4_envelope)) {
        
        cout << indent6 << "# quads split into 4 triangles: "
           << dualiso_info.triangulation.num_iso_cubes_tri_add_interior1
           << endl;

        if (output_info.flag_check_envelope &&
            output_info.flag_allow_tri4_envelope) {
          cout << indent8 << "# quads split into 4 triangles bcuz diagonal outside envelope: "
               << dualiso_info.triangulation.num_tri_add_interior1_with_diag_outside_envelope
               << endl;
        }
      }
    }
  }

  /// Report isosurface information.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE,
            typename DUALISO_INFO_TYPE>
  void report_iso_info
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const std::vector<COORD_TYPE> & vertex_coord, 
   const std::vector<VERTEX_INDEX> & plist, 
   const DUALISO_INFO_TYPE & dualiso_info)
  {
    const int TWO(2);
    const int THREE(3);
    const int FOUR(4);
    const int dimension = output_info.dimension;
    const int numv_per_isopoly = output_info.num_vertices_per_isopoly;
    const char * polytopes_name = "polytopes";
    
    using namespace std;

    VERTEX_INDEX numv = (vertex_coord.size())/dimension;
    VERTEX_INDEX num_poly = (plist.size())/numv_per_isopoly;
    IJK::PROCEDURE_ERROR error("report_iso_info");

    if (output_info.output_isovalue.size() < 1) {
      error.AddMessage
        ("Programming error. Output isovalue not set in output_info.");
      throw error;
    }
    
    if (output_info.flag_interval_volume) {
      cout << "  Interval volume [" 
           << output_info.isovalue[0] << ":"
           << output_info.isovalue[1] << "].  "
           << numv << " ivol vertices.  "
           << num_poly << " ivol polytopes." << endl;
    }
    else {

      if (output_info.Dimension() == TWO) {
        polytopes_name = "line segments";
      }
      else if (output_info.Dimension() == THREE) {

        // Default in 3D.
        polytopes_name = "polygons";

        if (output_info.mesh_type == SIMPLICIAL_COMPLEX) {
          polytopes_name = "triangles";
        }
        else if (output_info.mesh_type == CUBE_COMPLEX) {
          polytopes_name = "quadrilaterals";
        };
      }
      else if (output_info.Dimension() == FOUR) {
        polytopes_name = "hexahedra";

        // Triangulation not implemented.
      }

      cout << "  Isovalue " << output_info.output_isovalue[0] << ".  "
           << numv << " isosurface vertices.  "
           << num_poly << " isosurface " << polytopes_name
           << "." << endl;
    }

    if (!output_info.flag_use_stdout && !output_info.flag_silent &&
        !output_info.flag_terse) {             
      report_multi_isov_info(output_info, dualiso_info);
      report_isov_info(output_info, dualiso_info);
      report_triangulation_info(output_info, dualiso_info);
    }
  }


  /// Report information about isosurface quadrilaterals and triangles.
  template <typename OUTPUT_INFO_TYPE, typename DUALISO_DATA_TYPE,
            typename DUALISO_INFO_TYPE>
  void report_quad_tri_iso_info
  (const OUTPUT_INFO_TYPE & output_info, 
   const DUALISO_DATA_TYPE & dualiso_data,
   const std::vector<COORD_TYPE> & vertex_coord, 
   const std::vector<VERTEX_INDEX> & quad_list, 
   const std::vector<VERTEX_INDEX> & tri_list, 
   const DUALISO_INFO_TYPE & dualiso_info)
  {
    const int dimension = output_info.dimension;
    const int NUM_VERT_PER_QUAD(4);
    const int NUM_VERT_PER_TRI(3);
    const char * indent4 = "    ";

    using namespace std;

    const VERTEX_INDEX numv = (vertex_coord.size())/dimension;
    const VERTEX_INDEX num_quad = (quad_list.size())/NUM_VERT_PER_QUAD;
    const VERTEX_INDEX num_tri = (tri_list.size())/NUM_VERT_PER_TRI;

    cout << "  Isovalue " << output_info.isovalue[0] << ".  " 
         << numv << " isosurface vertices.  "
         << num_quad+num_tri << " isosurface polygons." << endl;
    cout << indent4 << num_quad << " quadrilaterals.  "
         << num_tri << " triangles." << endl;

    report_multi_isov_info(output_info, dualiso_info);
    report_triangulation_info(output_info, dualiso_info);
  }

  //@}


  // ******************************************************************
  //! @name WRITE ISOSURFACE VERTEX INFORMATION TO FILE
  // ******************************************************************

  //@{

  namespace {

    template <typename OUTPUT_INFO_TYPE, typename NTYPE>
    void output_isov_info_header
    (std::ostream & out, const OUTPUT_INFO_TYPE & output_info,
     const NTYPE num_isov)
    {
      using std::endl;
      
      if (output_info.flag_output_off) {
        out << "Isosurface vertex information for file: "
            << output_info.output_off_filename << endl;
      }
      else if (output_info.flag_output_ply) {
        out << "Isosurface vertex information for file: "
            << output_info.output_ply_filename << endl;      
      }
      else if (output_info.flag_output_fig) {
        out << "Isosurface vertex information for file: "
            << output_info.output_fig_filename << endl;      
      }
      else {
        out << "Isosurface vertex information for file: "
            << output_info.output_basename << ".*" << endl;
      }

      out << "Number of isosurface vertices: " << num_isov << endl;
      out << endl;
    }


    inline bool open_isov_info_file
    (const char * output_filename, std::ofstream & output_file)
    {
      using namespace std;
      
      output_file.open(output_filename, ios::out);

      if (output_file.is_open()) {
        return true;
      }
      else {
        cerr << "Warning: Unable to write isosurface vertex information to file."
             << endl;
        cerr << "  Unable to open file: " << output_filename << endl;
        cerr << "  Skipping write of isosurface vertex information to file."
             << endl;

        return false;
      }
    }
    
  };

  
  /*!
   *  @brief Output isosurface vertex info to out.
   */
  template <typename OUTPUT_INFO_TYPE, typename GRID_TYPE,
            typename DUAL_ISOV_TYPE, typename CTYPE>
  void output_isov_info
  (std::ostream & out,
   const OUTPUT_INFO_TYPE & output_info,
   const GRID_TYPE & grid,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const CTYPE vertex_coord[])
  {
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE CUBE_INDEX_TYPE;

    const int dimension = grid.Dimension();
    std::vector<CTYPE> cube_coord(dimension);
    CTYPE * cube_coord_ptr = IJK::vector2pointerNC(cube_coord);

    using std::endl;

    output_isov_info_header
      (out, output_info, isov_list.size());
      
    for (SIZE_TYPE isov = 0; isov < isov_list.size(); isov++) {
      const CTYPE * isov_coord = vertex_coord + isov*dimension;
      const CUBE_INDEX_TYPE icube =
        IJK::get_isovert_grid_cube_index(isov_list[isov]);
      const bool flag_cube_on_boundary = grid.IsCubeOnGridBoundary(icube);
        
      grid.ComputeCoord(icube, cube_coord_ptr);
      
      out << "Isosurface vertex: " << isov << endl;
      IJK::print_list
        (out, "  Coordinates: ", isov_coord, dimension, "\n");
      grid.PrintIndexAndCoord(out, "  In cube: ", icube, "");
      if (flag_cube_on_boundary)
        { out << ". Boundary cube."; }
      out << endl;
    }
  }


  /*!
   *  @brief Write to file isosurface vertex info,
   *    including lookup table info.
   */
  template <typename OUTPUT_INFO_TYPE, typename GRID_TYPE,
            typename DUAL_ISOV_TYPE, typename CTYPE>
  void write_isov_info
  (const OUTPUT_INFO_TYPE & output_info,
   const GRID_TYPE & grid,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const CTYPE vertex_coord[])
  {
    const char * output_filename =
      output_info.write_isov_filename.c_str();
    std::ofstream output_file;

    using namespace std;

    if (!open_isov_info_file(output_filename, output_file)) {
      // Open failed. Skip write.
      return;
    }

    output_isov_info
      (output_file, output_info, grid, isov_list, vertex_coord);
    
    output_file.close();

    if (!output_info.flag_silent)
      { cout << "Wrote output to file: " << output_filename << endl; }
  }


  /*!
   *  @overload
   *  @brief Write to file isosurface vertex info,
   *    including lookup table info. (C++ vector vertex_coord[].)
   *  - Version using C++ STL vector vertex_coord[].
   */
  template <typename OUTPUT_INFO_TYPE, typename GRID_TYPE,
            typename DUAL_ISOV_TYPE, typename CTYPE>
  void write_isov_info
  (const OUTPUT_INFO_TYPE & output_info,
   const GRID_TYPE & grid,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const std::vector<CTYPE> & vertex_coord)
  {
    const CTYPE * vertex_coord_ptr =
      IJK::vector2pointer(vertex_coord);

    write_isov_info(output_info, grid, isov_list, vertex_coord_ptr);
  }

  //@}

  
  // ******************************************************************
  //! @name WRITE ISOV INFO TO FILE, ALLOW MULTIPLE ISOV PER CUBE
  // ******************************************************************

  //@{
  
  /*!
   *  @brief Output isosurface vertex info to out, including dual isosurface
   *    lookup table information.
   */
  template <typename OUTPUT_INFO_TYPE, typename GRID_TYPE,
            typename ISODUAL_TABLE_TYPE,
            typename DUAL_ISOV_TYPE, typename CTYPE>
  void output_multi_isov_info
  (std::ostream & out,
   const OUTPUT_INFO_TYPE & output_info,
   const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const CTYPE vertex_coord[])
  {
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;
    typedef typename DUAL_ISOV_TYPE::CUBE_INDEX_TYPE CUBE_INDEX_TYPE;
    typedef typename DUAL_ISOV_TYPE::TABLE_INDEX_TYPE TABLE_INDEX_TYPE;

    const int dimension = grid.Dimension();
    const int DIM3(3);
    std::vector<CTYPE> cube_coord(dimension);
    CTYPE * cube_coord_ptr = IJK::vector2pointerNC(cube_coord);

    using std::endl;

    output_isov_info_header
      (out, output_info, isov_list.size());
      
    for (SIZE_TYPE isov = 0; isov < isov_list.size(); isov++) {
      const CTYPE * isov_coord = vertex_coord + isov*dimension;
      const CUBE_INDEX_TYPE icube =
        IJK::get_isovert_grid_cube_index(isov_list[isov]);
      const bool flag_cube_on_boundary = grid.IsCubeOnGridBoundary(icube);
      const TABLE_INDEX_TYPE table_index = isov_list[isov].table_index;
      const int ipatch = isov_list[isov].patch_index;
      const int num_isov = isodual_table.NumIsoVertices(table_index);
      const int num_incident_isopoly =
        isodual_table.NumIncidentIsoPoly(table_index, ipatch);
        
      grid.ComputeCoord(icube, cube_coord_ptr);
      
      out << "Isosurface vertex: " << isov << endl;
      IJK::print_list
        (out, "  Coordinates: ", isov_coord, dimension, "\n");
      grid.PrintIndexAndCoord(out, "  In cube: ", icube, "");
      if (flag_cube_on_boundary)
        { out << ". Boundary cube."; }
      out << endl;

      out << "  Table index: " << table_index
          << "  Patch index: " << ipatch
          << "  Degree: " << num_incident_isopoly
          << "  Num isov in cube: " << num_isov
          << endl;

      if (dimension == DIM3) {
        if ((isodual_table.NumAmbiguousFacets(table_index) == 1) &
            (num_isov == 1)) {
          out << "  Is doubly connected: True" << endl;
        }
      }
    }
  }


  /*!
   *  @brief Write to file isosurface vertex info,
   *    including lookup table info.
   */
  template <typename OUTPUT_INFO_TYPE, typename GRID_TYPE,
            typename ISODUAL_TABLE_TYPE, typename DUAL_ISOV_TYPE,
            typename CTYPE>
  void write_multi_isov_info
  (const OUTPUT_INFO_TYPE & output_info,
   const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const CTYPE vertex_coord[])
  {
    const char * output_filename =
      output_info.write_isov_filename.c_str();
    std::ofstream output_file;

    using namespace std;

    if (!open_isov_info_file(output_filename, output_file)) {
      // Open failed. Skip write.
      return;
    }

    output_multi_isov_info
      (output_file, output_info, grid, isodual_table, isov_list,
       vertex_coord);
    
    output_file.close();

    if (!output_info.flag_silent)
      { cout << "Wrote output to file: " << output_filename << endl; }
  }


  /*!
   *  @overload
   *  @brief Write to file isosurface vertex info,
   *    including lookup table info. (C++ vector vertex_coord[].)
   *  - Version using C++ STL vector vertex_coord[].
   */
  template <typename OUTPUT_INFO_TYPE, typename GRID_TYPE,
            typename ISODUAL_TABLE_TYPE, typename DUAL_ISOV_TYPE,
            typename CTYPE>
  void write_multi_isov_info
  (const OUTPUT_INFO_TYPE & output_info,
   const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const std::vector<CTYPE> & vertex_coord)
  {
    const CTYPE * vertex_coord_ptr =
      IJK::vector2pointer(vertex_coord);

    write_multi_isov_info
      (output_info, grid, isodual_table, isov_list, vertex_coord_ptr);
  }

  //@}
  
}

#endif
