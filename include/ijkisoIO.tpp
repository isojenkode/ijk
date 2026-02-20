/*!
 *  @file ijkisoIO.tpp
 *  @brief Templates for isosurface IO routines.
 *  - Class IO_PARAM.
 *  - Output formats: Geomview .off, OpenInventor .iv (3D), Fig .fig (2D)
 *    Stanford .ply, and Visualization Toolkit, .vtk.
 *  - Requires compilation with c++17 or later version.
 *  - Version 0.6.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2023-2025 Rephael Wenger

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

#ifndef _IJKIO_PARAM_
#define _IJKIO_PARAM_

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "ijk.tpp"
#include "ijkgrid_nrrd.tpp"
#include "ijkIO.tpp"
#include "ijkstring.tpp"


namespace IJK {

  // *****************************************************************
  // Read Nrrd file into SCALAR_GRID data structure
  // *****************************************************************

  /*!
   *  @brief Read Nrrd file into scalar_grid.
   *  @tparam SCALAR_GRID_TYPE Scalar grid type.
   *    - Must include SetSpacing() member function.
   */
  template <typename SCALAR_GRID_TYPE, typename DATA_TYPE,
            typename AXIS_SIZE_TYPE>
  void read_nrrd_file
  (const char * input_filename, SCALAR_GRID_TYPE & scalar_grid,
   IJK::NRRD_DATA<DATA_TYPE,AXIS_SIZE_TYPE> & nrrd_header)
  {
    typedef typename SCALAR_GRID_TYPE::SPACING_TYPE SPACING_TYPE;

    IJK::GRID_NRRD_IN<DATA_TYPE,AXIS_SIZE_TYPE> nrrd_in;
    IJK::PROCEDURE_ERROR error("read_nrrd_file");

    nrrd_in.ReadScalarGrid(input_filename, scalar_grid, nrrd_header, error);
    if (nrrd_in.ReadFailed()) { throw error; }

    std::vector<SPACING_TYPE> grid_spacing;
    nrrd_header.GetSpacing(grid_spacing);
    scalar_grid.SetSpacing(grid_spacing);
  }

  
  /*!
   *  @brief Read Nrrd file into scalar_grid.
   *  - Version using std::string for input_filename.
   */
  template <typename SCALAR_GRID_TYPE, typename DATA_TYPE,
            typename AXIS_SIZE_TYPE>
  void read_nrrd_file
  (const std::string & input_filename, SCALAR_GRID_TYPE & scalar_grid,
   IJK::NRRD_DATA<DATA_TYPE,AXIS_SIZE_TYPE> & nrrd_header)
  {
    read_nrrd_file(input_filename.c_str(), scalar_grid, nrrd_header);
  }


  // *****************************************************************
  // Class ISO_IO_PARAM_BASE
  // *****************************************************************

  /*!
   *  @base IO Parameters
   *  @tparam _OTHER_DATA Any other data that is stored
   *     in the same class.
   */
  template <typename _SCALAR_TYPE, typename _COORD_TYPE,
            typename _OTHER_DATA>
  class ISO_IO_PARAM_BASE:
    public OUTPUT_PARAM_BASE<_COORD_TYPE,_OTHER_DATA> {

  public:
    typedef _COORD_TYPE COORD_TYPE;
    typedef _SCALAR_TYPE SCALAR_TYPE;

  protected:

    /// Initialize data structure.
    void Init();


  public:

    // Input parameters.
    std::vector<SCALAR_TYPE> isovalue;     ///< List of isovalues.
    std::vector<std::string> isovalue_string;
    std::vector<COORD_TYPE> grid_spacing;
    int subsample_resolution;
    int supersample_resolution;
    int region_length;

    /// @brief List of high resolution arguments.
    /// - Currently ignored with ijkdual.
    std::vector<std::string> high_resolution_string;

    // File/directory names.
    std::string input_filename;
    
    // *** DEPRECATED ***
    std::string report_isov_filename;

    // Flags.
    bool flag_report_time;
    bool flag_report_time_detailed;
    bool flag_report_info;
    bool flag_nowrite;
    bool flag_terse;
    bool flag_verbose;
    bool flag_no_warn;
    bool flag_subsample;
    bool flag_supersample;
    bool flag_report_all_isov;

    /// If true, output filename includes isovalue.
    bool flag_label_with_isovalue;

    /// If true, triangulate mesh.
    bool flag_trimesh;

    /*!
     *  @brief If true, color isosurface boundary vertices. 
     *  - Currently not used (always false).
     */
    bool flag_color_boundary_vert;

    // Flags for options set.

    /// Is quadrilateral edge intersection method set?
    bool is_qei_method_set;

    bool is_tri4_position_method_set;
    bool is_file_format_set;
    bool is_random_seed_set;
    bool is_isov_separation_method_set;

    /// @brief Construct output filenames.
    /// @param i Construct filenames for isovalue i.
    void ConstructOutputFilenames
      (const int i, const bool flag_interval_volume = false);


  public:
    ISO_IO_PARAM_BASE() { Init(); };
    ~ISO_IO_PARAM_BASE() { Init(); };
  };


  // *****************************************************************
  // Class ISO_OUTPUT_PARAM_BASE
  // *****************************************************************

  /*!
   *  @brief Output parameters.
   *  @tparam _ISO_IO_PARAM IO parameters class.
   *    - Usually derived from ISO_IO_PARAM_BASE.
   *    - Must include type definition of SCALAR_TYPE.
   */
  template <typename _SCALAR_TYPE, typename _COORD_TYPE,
            typename _OTHER_DATA>
  class ISO_OUTPUT_PARAM_BASE:
    public ISO_IO_PARAM_BASE<_SCALAR_TYPE,_COORD_TYPE,_OTHER_DATA> {

  public:
    typedef _COORD_TYPE COORD_TYPE;
    typedef _SCALAR_TYPE SCALAR_TYPE;
    
  protected:
    void Init();

  public:
    /// Number of vertices per output isosurface polytope.
    int num_vertices_per_isopoly;

    /// Isovalues in this output.
    std::vector<SCALAR_TYPE> output_isovalue;

    /// Constructor.
    ISO_OUTPUT_PARAM_BASE() { Init(); };

    /// Destructor.
    ~ISO_OUTPUT_PARAM_BASE() { Init(); };


    // Add comment routines.

    /// Add scalar data file comment to array comment[].
    void AddScalarDataFileComment(std::vector<std::string> & comment) const;
    
    /// Add isovalue comment to array comment[].
    void AddIsovalueComment(std::vector<std::string> & comment) const;

    /// Add high_resolution_region comment(s) to array comment[].
    void AddHighResolutionComments(std::vector<std::string> & comment) const;

    /// @brief Add subsample resolution comment to array comment[].
    /// - Only adds comment if flag_subsample is true.
    void AddSubsampleResolutionComment
      (std::vector<std::string> & comment) const;

    /// @brief Add supersample resolution comment to array comment[].
    /// - Only adds comment if flag_supersample is true.
    void AddSupersampleResolutionComment
      (std::vector<std::string> & comment) const;

    /// @brief Set from ISO_IO_PARAM_BASE.
    void Set
    (const IJK::ISO_IO_PARAM_BASE<_SCALAR_TYPE,_COORD_TYPE,_OTHER_DATA>
     & io_param)
    { ISO_IO_PARAM_BASE<_SCALAR_TYPE,_COORD_TYPE,_OTHER_DATA>::operator=
        (io_param); }

  };


  // *****************************************************************
  // Class ISO_IO_PARAM_BASE member functions
  // *****************************************************************

  // Initialize ISO_IO_PARAM_BASE.
  template <typename SCALAR_TYPE, typename COORD_TYPE,
            typename ISO_DATA_PARAM>
  void ISO_IO_PARAM_BASE<SCALAR_TYPE,COORD_TYPE,ISO_DATA_PARAM>::Init()
  {
    // Initialize input parameters.
    isovalue.clear();
    isovalue_string.clear();
    input_filename .clear();
    subsample_resolution = 2;
    supersample_resolution = 2;
    region_length = 1;

    // Initialize IO flags.
    flag_report_time = false;
    flag_report_time_detailed = false;
    flag_report_info = false;
    flag_nowrite = false;
    flag_terse = false;
    flag_verbose = false;
    flag_no_warn = false;
    flag_supersample = false;
    flag_subsample = false;
    flag_report_all_isov = false;
    flag_label_with_isovalue = false;
    flag_trimesh = false;

    // color isosurface boundary vertices
    flag_color_boundary_vert = false;

    is_qei_method_set = false;
    is_tri4_position_method_set = false;
    is_file_format_set = false;
    is_random_seed_set = false;
    is_isov_separation_method_set = false;
  }


  template <typename SCALAR_TYPE, typename COORD_TYPE,
            typename ISO_DATA_PARAM>
  void ISO_IO_PARAM_BASE<SCALAR_TYPE,COORD_TYPE,ISO_DATA_PARAM>::
  ConstructOutputFilenames(const int i, const bool flag_interval_volume)
  {
    typedef typename std::vector<std::string>::size_type SIZE_TYPE;

    using std::string;

    string basename, suffix;
    string ofilename;

    // create output filename
    if (this->output_basename == "") {

      const string fname =
        std::filesystem::path(input_filename).filename().string();

      // construct output filename
      IJK::split_string(fname, '.', basename, suffix);
      if (suffix == "nrrd" || suffix == "nhdr")
        { this->output_basename = basename; }
      else { this->output_basename = string(input_filename); }
    }

    this->output_basename =
      this->output_filename_prefix + this->output_basename +
      this->output_filename_suffix;

    if (flag_interval_volume) {
      if (flag_label_with_isovalue || isovalue_string.size() > 2) {
        if (SIZE_TYPE(i+1) < isovalue_string.size()) {
          this->output_basename +=
            string(".") + string("ivol=") + isovalue_string[i] +
            string("_") + isovalue_string[i+1];
        }
      }
    }
    else {
      if (flag_label_with_isovalue || isovalue_string.size() > 1) {
        this->output_basename += string(".") + string("isov=") + isovalue_string[i];
      }
    }

    this->output_off_filename = this->output_basename + ".off";
    this->output_ply_filename = this->output_basename + ".ply";
    this->output_fig_filename = this->output_basename + ".fig";
    this->output_iv_filename = this->output_basename + ".iv";
  }


  // *****************************************************************
  // Class ISO_OUTPUT_PARAM_BASE member functions
  // *****************************************************************

  // TO BE CONTINUED...
  template <typename _SCALAR_TYPE, typename _COORD_TYPE,
            typename _OTHER_DATA>
  void ISO_OUTPUT_PARAM_BASE<_SCALAR_TYPE,_COORD_TYPE,_OTHER_DATA>::
  Init()
  {
    const int NUM_VERT_PER_TRIANGLE(3);

    num_vertices_per_isopoly = NUM_VERT_PER_TRIANGLE;
  }


  template <typename _SCALAR_TYPE, typename _COORD_TYPE,
            typename _OTHER_DATA>
  void ISO_OUTPUT_PARAM_BASE<_SCALAR_TYPE,_COORD_TYPE,_OTHER_DATA>::  
  AddScalarDataFileComment(std::vector<std::string> & comment) const
  {
    if (this->input_filename != "") {
      const std::string fname =
        std::filesystem::path(this->input_filename).filename().string();
      const std::string data_file_string =
        "ScalarDataFile: " + fname;
      comment.push_back(data_file_string);
    }
  }
  
  
  template <typename _SCALAR_TYPE, typename _COORD_TYPE,
            typename _OTHER_DATA>
  void ISO_OUTPUT_PARAM_BASE<_SCALAR_TYPE,_COORD_TYPE,_OTHER_DATA>::
  AddIsovalueComment(std::vector<std::string> & comment) const
  {
    if (this->output_isovalue.size() == 1) {
      std::string isovalue_string;
      IJK::create_string
        ("Isovalue: ", this->output_isovalue[0], isovalue_string);
      comment.push_back(isovalue_string);
    }
    else if (this->output_isovalue.size() > 1) {
      std::string isovalue0_string, isovalue1_string;
      std::string isovalue_string;
      IJK::val2string(this->output_isovalue[0], isovalue0_string);
      IJK::val2string(this->output_isovalue[1], isovalue1_string);
      isovalue_string = "Isovalues: " + isovalue0_string + " " + isovalue1_string;
      comment.push_back(isovalue_string);
    }
  }


  // Add high_resolution_region comment(s) to array comment[].
  template <typename _SCALAR_TYPE, typename _COORD_TYPE,
            typename _OTHER_DATA>
  void ISO_OUTPUT_PARAM_BASE<_SCALAR_TYPE,_COORD_TYPE,_OTHER_DATA>::
  AddHighResolutionComments(std::vector<std::string> & comment) const
  {
    for (unsigned int i = 0; i < this->high_resolution_string.size(); i++) {
      std::string str_i;
      IJK::val2string(i, str_i);
      std::string highres_string_i =
        std::string("HighResRegion") + str_i + ": " +
        this->high_resolution_string[i];
      comment.push_back(highres_string_i);
    }
  }


  // Add subsample resolution to array comment[].
  template <typename _SCALAR_TYPE, typename _COORD_TYPE,
            typename _OTHER_DATA>
  void ISO_OUTPUT_PARAM_BASE<_SCALAR_TYPE,_COORD_TYPE,_OTHER_DATA>::
  AddSubsampleResolutionComment
    (std::vector<std::string> & comment) const
  {
    if (this->flag_subsample) {
      std::string str;
      IJK::create_string("SubsampleResolution: ",
                         this->subsample_resolution, str);
      comment.push_back(str);
    }
  }


  // Add supersample resolution to array comment[].
  template <typename _SCALAR_TYPE, typename _COORD_TYPE,
            typename _OTHER_DATA>
  void ISO_OUTPUT_PARAM_BASE<_SCALAR_TYPE,_COORD_TYPE,_OTHER_DATA>::
  AddSupersampleResolutionComment
    (std::vector<std::string> & comment) const
  {
    if (this->flag_supersample) {
      std::string str;
      IJK::create_string("SupersampleResolution: ",
                         this->supersample_resolution, str);
      comment.push_back(str);
    }
  }

}

#endif
  
