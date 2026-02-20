/*!
 *  @file ijkisoinfo.cpp
 *  @brief Get info about a level set in a scalar grid.
 *  - Input: isovalue and scalar grid.
 *  - Version 0.6.0
 */


/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2024-2025 Rephael Wenger

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


#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

#include "ijk.tpp"
#include "ijkcommand_line.tpp"
#include "ijkcoord.tpp"
#include "ijkisopoly.tpp"
#include "ijkprint.tpp"
#include "ijkscalar_grid.tpp"
#include "ijkstring.tpp"

#include "ijkgrid_macros.h"
#include "ijkNrrd.h"

#include "ijkisoinfo.h"

using namespace IJKISOINFO;
using namespace std;

// output routines
void report_info
(const INPUT_ARG & input_arg, const SGRID & scalar_grid,
 const int scalar_type);

// misc routines
void memory_exhaustion();
void help(), usage_error();
void parse_command_line
(int argc, char **argv, INPUT_ARG & input_arg);


// *******************************************************************
// MAIN
// *******************************************************************

int main(int argc, char **argv)
{
  INPUT_ARG input_arg;
  SGRID_NRRD_IN grid_nrrd_in;
  SGRID scalar_grid;
  IJK::ERROR error("ijklevelsetinfo");

  try {

    std::set_new_handler(memory_exhaustion);

    parse_command_line(argc, argv, input_arg);

    grid_nrrd_in.ReadScalarGrid
      (input_arg.nrrd_filename, scalar_grid, error);

    if (scalar_grid.Dimension() < 1) {
      cerr << "Illegal scalar grid dimension.  Dimension must be at least 1." << endl;
      exit(20);
    };

    const int scalar_type = grid_nrrd_in.DataPtrConst()->type;

    if (!input_arg.flag_subsample) {
      report_info(input_arg, scalar_grid, scalar_type);
    }
    else {
      // subsample grid
      SGRID subsampled_grid;
      subsampled_grid.Subsample
        (scalar_grid, input_arg.subsample_resolution);
      report_info(input_arg, subsampled_grid, scalar_type);
    }
  }
  catch(IJK::ERROR error) {
    if (error.NumMessages() == 0) {
      cerr << "Unknown error." << endl;
    }
    else { error.Print(cerr); }
    cerr << "Exiting." << endl;
    exit(20);
  }
  catch (...) {
    cerr << "Unknown error." << endl;
    exit(50);
  };

}


// *******************************************************************
// Count Routines
// *******************************************************************

/// Return true if grid vertex scalar value is >= isovalue.
inline bool is_positive
(const SCALAR_TYPE isovalue, const SGRID & scalar_grid,
 const VERTEX_INDEX iv)
{
  return ((scalar_grid.Scalar(iv) >= isovalue)? true: false);
}

  
template <bool flag_pos>
bool is_isolated_grid_vertex
(const SCALAR_TYPE isovalue, const SGRID & scalar_grid,
 const VERTEX_INDEX ivA)
{
  if (is_positive(isovalue, scalar_grid, ivA) != flag_pos) {
    // Vertex has label !flag_pos.
    return false;
  }
  
  for (int d = 0; d < scalar_grid.Dimension(); d++) {
    for (int side = 0; side < 2; side++) {
      const VERTEX_INDEX ivB =
        scalar_grid.AdjacentVertex(ivA, d, side);
      if (is_positive(isovalue, scalar_grid, ivB) == flag_pos) {
        // Vertex is not isolated.
        return false;
      }
    }
  }

  return true;
}


template <bool flag_pos>
int count_num_isolated_grid_vertices
(const SCALAR_TYPE isovalue, const SGRID & scalar_grid)
{
  int num_isolated = 0;
  
  IJK_FOR_EACH_INTERIOR_GRID_VERTEX(iv, scalar_grid, VERTEX_INDEX) {
    if (is_isolated_grid_vertex<flag_pos>
        (isovalue, scalar_grid, iv))
      { num_isolated++; }
  }

  return num_isolated;
}


int count_num_active_grid_cubes
(const SCALAR_TYPE isovalue, const SGRID & scalar_grid)
{
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  const NUM_TYPE num_cube_vertices = scalar_grid.NumCubeVertices();
  const VERTEX_INDEX * cube_vertex_increment =
    scalar_grid.CubeVertexIncrement();

  int num_active = 0;

  IJK_FOR_EACH_GRID_CUBE(icube, scalar_grid, VERTEX_INDEX) {
    if (IJK::intersects_poly
        (scalar, isovalue, icube, cube_vertex_increment, 
         num_cube_vertices))
    { num_active++; }
  }

  return num_active;
}

  
// *******************************************************************
// Report/Output Routines
// *******************************************************************

// Forward definitions.
void output_num_active_grid_cubes
(const SCALAR_TYPE isovalue, const SGRID & scalar_grid);
void output_num_isolated_grid_vertices
(const SCALAR_TYPE isovalue, const SGRID & scalar_grid);


void report_info
(const INPUT_ARG & input_arg, const SGRID & scalar_grid,
 const int scalar_type)
{
  for (int i = 0; i < input_arg.isovalue.size(); i++) {
    const SCALAR_TYPE isovalue = input_arg.isovalue[i];
    cout << "Isovalue: " << isovalue << endl;
    output_num_active_grid_cubes (isovalue, scalar_grid);
    output_num_isolated_grid_vertices(isovalue, scalar_grid);
    
    cout << endl;
  }
}


void output_num_active_grid_cubes
(const SCALAR_TYPE isovalue, const SGRID & scalar_grid)
{
  const NUM_TYPE num_active_grid_cubes = 
    count_num_active_grid_cubes(isovalue, scalar_grid);

  cout << "  # active grid cubes: " << num_active_grid_cubes << endl;
}

void output_num_isolated_grid_vertices
(const SCALAR_TYPE isovalue, const SGRID & scalar_grid)
{
  const NUM_TYPE num_pos_isolated =
    count_num_isolated_grid_vertices<true>(isovalue, scalar_grid);
  const NUM_TYPE num_neg_isolated =
    count_num_isolated_grid_vertices<false>(isovalue, scalar_grid);

  cout << "  # pos isolated grid vertices: "
       << num_pos_isolated << endl;
  cout << "  # neg isolated grid vertices: "
       << num_neg_isolated << endl;  
}


// *******************************************************************
// Misc Routines
// *******************************************************************

// Local namespace
namespace {


  /// @brief Usage message.
  void usage_msg(std::ostream & out)
  {
    out << "Usage: ijkisoinfo {isovalue} {nrrd file}" << endl;
  }

               
  /// @brief Get command line option.
  bool get_option(char * arg, COMMAND_LINE_OPTION_TYPE & opt)
  {
    const std::string s = arg;
    
    if (s == "-subsample") {
      opt = SUBSAMPLE_OPT;
      return true;
    }

    return false;
  }
  
  
  /// @brief Store isovalues and input filename in io_info.
  void process_isovalues_and_input_filename
  (const int argc, char **argv, const int iarg, INPUT_ARG & input_arg)
  {    
    // check for more parameter tokens
    for (int j = iarg; j+1 < argc; j++) {

      COMMAND_LINE_OPTION_TYPE optA;
      if (get_option(argv[j], optA) ||
          (argv[j][0] == '-' && !IJK::is_type<float>(argv[j]))) {
        // argv[iarg] is not an isovalue
        cerr << "Usage error. Illegal option: " << argv[iarg] << endl;
        cerr << endl;
        usage_error();
      }
    }
  
    if (iarg+2 > argc) {
      cerr << "Error.  Missing input isovalue or input file name." << endl;
      cerr << endl;
      usage_error();
    };

    // store isovalues
    for (int j = iarg; j+1 < argc; j++) {
      SCALAR_TYPE value;
      if (!IJK::string2val(argv[j], value)) {
        cerr << "Error. \"" << argv[j] << "\" is not a valid input isovalue." 
             << endl;
        usage_error();
      }

      input_arg.isovalue_string.push_back(argv[j]);
      input_arg.isovalue.push_back(value);
    }
    
    input_arg.nrrd_filename = argv[argc-1];
  }

  
};


void memory_exhaustion()
{
  cerr << "Error: Out of memory.  Terminating program." << endl;
  exit(10);
}


void parse_command_line
(int argc, char **argv, INPUT_ARG & input_arg)
{
  IJK::ERROR error;

  int iarg = 1;

  while (iarg < argc && argv[iarg][0] == '-') {

    COMMAND_LINE_OPTION_TYPE optA;
    if (!get_option(argv[iarg], optA)) {
      // Unknown parameter.  Possibly negative scalar value.
      break;
    }

    switch (optA) {

    case SUBSAMPLE_OPT:
      input_arg.subsample_resolution = 
        get_arg_int(iarg, argc, argv, error);
      input_arg.flag_subsample = true;
      iarg++;
      break;

    default:
      error.AddMessage
        ("Programming error. Option ", argv[iarg],
         " is not implemented.");
      throw error;
    }

    iarg++;
  }

  process_isovalues_and_input_filename(argc, argv, iarg, input_arg);
}


void usage_error()
{
  usage_msg(cerr);
  exit(-1);
}


void help()
{
  usage_msg(cout);
  cout << "  Print information about level set determined by {isovalue}." << endl;
  exit(0);
}


// *******************************************************************
// Class INPUT_ARG member functions
// *******************************************************************

void INPUT_ARG::Init()
{
  flag_subsample = false;
  subsample_resolution = 2;
}
