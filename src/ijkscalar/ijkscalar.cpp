/*!
 *  @file ijkscalar.cpp
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

#include <fstream>
#include <iostream>
#include <vector>

#include "ijkcommand_line.tpp"
#include "ijkcoord.tpp"
#include "ijkgrid_nrrd.tpp"
#include "ijkstring.tpp"

#include "ijkgrid_macros.h"

#include "ijkscalar.h"

using namespace std;
using namespace IJK;
using namespace IJKSCALAR;

// Commands
typedef enum
  { MERGE_MIN, MERGE_MAX, SCALAR_OP, HELP, UNDEFINED_COMMAND } 
  IJKSCALAR_COMMAND;

typedef enum
  { NEGATE_SCALAR, ADD_SCALAR, MULTIPLY_SCALAR, SET_SCALAR,
    SET_BOUNDARY_SCALAR,
    SUBSAMPLE_SCALAR, NONUNIFORM_SUBSAMPLE_SCALAR, 
    SUPERSAMPLE_SCALAR, NONUNIFORM_SUPERSAMPLE_SCALAR, 
    UNDEFINED_SCALAR_OP } SCALAR_OP_TYPE;

// global variables
const int TWO_INPUT_GRIDS = 2;
char * scalar_filename[TWO_INPUT_GRIDS] = {NULL};
char * gradient_filename[TWO_INPUT_GRIDS] = {NULL};
string output_scalar_filename;
IJKSCALAR_COMMAND command(UNDEFINED_COMMAND);
SCALAR_OP_TYPE scalar_op(UNDEFINED_SCALAR_OP);
SCALAR_TYPE scalar_operand(0);
int scalar_operandI(0);
std::vector<int> scalar_operandIV;
bool flag_gradients(false);

// compute routines
void merge_scalar
(const SGRID scalar_grid[2], const IJKSCALAR_COMMAND command,
 SGRID & output_scalar_grid);
void merge_scalar
(const SGRID scalar_grid[2], const GRADIENT_GRID gradient_grid[2],
 const IJKSCALAR_COMMAND command,
 SGRID & output_scalar_grid, GRADIENT_GRID & output_gradient_grid);
void apply_scalar_op
(const SGRID & scalar_grid, const SCALAR_OP_TYPE scalar_op,
 const SCALAR_TYPE scalar_operand, SGRID & output_scalar_grid);
void apply_scalar_op
(const SGRID & scalar_grid, const GRADIENT_GRID & gradient_grid,
 const SCALAR_OP_TYPE scalar_op, const SCALAR_TYPE scalar_operand, 
 SGRID & output_scalar_grid, GRADIENT_GRID & output_gradient_grid);
void set_boundary_scalar(const SCALAR_TYPE val, SGRID & scalar_grid);

// read/write routines
void read_scalar_grid(const char * filename, SGRID & scalar_grid);
void read_gradient_grid(const char * filename, GRADIENT_GRID & gradient_grid);
void read_gradient_grid(const string & filename, GRADIENT_GRID & gradient_grid);
void write_scalar_grid
(const string & output_filename, const SGRID & scalar_grid);
void write_gradient_grid
(const string & output_filename, const GRADIENT_GRID & gradient_grid);


// check routines.
bool check_scalar_operands
(const SGRID & scalar_grid, IJK::ERROR & error);

// misc routines
void memory_exhaustion();
void help(), usage_error();
void parse_command_line(int argc, char **argv);
void parse_filename
(const string & filename, string & prefix, string & suffix);
void split_string(const string & s, const char c,
                  string & prefix, string & suffix);
void get_gradient_filename
(const char * scalar_filename, const char * gradient_filename, 
 string & gradient_filename_str);
void get_gradient_filename
(const string & scalar_filename, const char * gradient_filename, 
 string & gradient_filename_str);


// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
  IJK::PROCEDURE_ERROR error("ijkscalar");

  SGRID output_scalar_grid;
  GRADIENT_GRID output_gradient_grid;
  string gradient_filename_str[2];
  string output_gradient_filename;

  try {

    std::set_new_handler(memory_exhaustion);

    parse_command_line(argc, argv);

    if (flag_gradients) {
      for (int i = 0; i < 2; i++) {
        get_gradient_filename
          (scalar_filename[i], gradient_filename[i], gradient_filename_str[i]);
      }
    }


    if (output_scalar_filename == "") 
      { output_scalar_filename = "out.nrrd"; }


    if (command == SCALAR_OP) {
      SGRID scalar_grid;
      GRADIENT_GRID gradient_grid;

      read_scalar_grid(scalar_filename[0], scalar_grid);

      if (!check_scalar_operands(scalar_grid, error)) 
        { throw error; }

      output_scalar_grid.SetSize(scalar_grid);

      if (flag_gradients) {

        read_gradient_grid(gradient_filename_str[0], gradient_grid);
        gradient_grid.SetSize(scalar_grid, scalar_grid.Dimension());

        output_gradient_grid.SetSize(gradient_grid);

        apply_scalar_op
          (scalar_grid, gradient_grid, scalar_op, scalar_operand, 
           output_scalar_grid, output_gradient_grid);

        write_scalar_grid(output_scalar_filename, output_scalar_grid);
        get_gradient_filename
          (output_scalar_filename, NULL, output_gradient_filename);
        write_gradient_grid(output_gradient_filename, output_gradient_grid);
      }
      else {
        apply_scalar_op
          (scalar_grid, scalar_op, scalar_operand, output_scalar_grid);
        write_scalar_grid(output_scalar_filename, output_scalar_grid);
      }

    }
    else {
      SGRID scalar_grid[TWO_INPUT_GRIDS];
      GRADIENT_GRID gradient_grid[TWO_INPUT_GRIDS];

      for (int i = 0; i < 2; i++) {
        read_scalar_grid(scalar_filename[i], scalar_grid[i]);
      }

      if (!scalar_grid[0].CheckDimension
          (scalar_grid[1], scalar_filename[0], scalar_filename[1],
           error)) { throw error; }

      output_scalar_grid.SetSize(scalar_grid[0]);

      if (flag_gradients) {

        for (int i = 0; i < 2; i++) {
          read_gradient_grid(gradient_filename_str[i], gradient_grid[i]);
          gradient_grid[i].SetSize(scalar_grid[i], scalar_grid[i].Dimension());
        }

        output_gradient_grid.SetSize(gradient_grid[0]);


        merge_scalar(scalar_grid, gradient_grid, command, 
                     output_scalar_grid, output_gradient_grid);
        write_scalar_grid(output_scalar_filename, output_scalar_grid);
        get_gradient_filename
          (output_scalar_filename, NULL, output_gradient_filename);
        write_gradient_grid(output_gradient_filename, output_gradient_grid);
      }
      else {

        merge_scalar(scalar_grid, command, output_scalar_grid);
        write_scalar_grid(output_scalar_filename, output_scalar_grid);
      }
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



// **************************************************
// Compute routines
// **************************************************

void compute_scalar_min
(const SGRID scalar_grid[2], SGRID & output_scalar_grid);
void compute_scalar_min
(const SGRID scalar_grid[2], const GRADIENT_GRID gradient_grid[2],
 SGRID & output_scalar_grid, GRADIENT_GRID & output_gradient_grid);
void compute_scalar_max
(const SGRID scalar_grid[2], SGRID & output_scalar_grid);
void compute_scalar_max
(const SGRID scalar_grid[2], const GRADIENT_GRID gradient_grid[2],
 SGRID & output_scalar_grid, GRADIENT_GRID & output_gradient_grid);

void merge_scalar
(const SGRID scalar_grid[2], const IJKSCALAR_COMMAND command,
 SGRID & output_scalar_grid)
{
  IJK::PROCEDURE_ERROR error("merge_scalar");

  switch(command) {

  case MERGE_MIN:
    compute_scalar_min(scalar_grid, output_scalar_grid);
    break;

  case MERGE_MAX:
    compute_scalar_max(scalar_grid, output_scalar_grid);
    break;

  case UNDEFINED_COMMAND:
    error.AddMessage("Programming error.  Undefined merge scalar command.");
    throw error;

  default:
    error.AddMessage("Programming error.  Illegal merge scalar command.");
    throw error;
  }
}

void merge_scalar
(const SGRID scalar_grid[2], const GRADIENT_GRID gradient_grid[2],
 const IJKSCALAR_COMMAND command,
 SGRID & output_scalar_grid, GRADIENT_GRID & output_gradient_grid)
{
  IJK::PROCEDURE_ERROR error("merge_scalar");

  switch(command) {

  case MERGE_MIN:
    compute_scalar_min(scalar_grid, gradient_grid,
                       output_scalar_grid, output_gradient_grid);
    break;

  case MERGE_MAX:
    compute_scalar_max(scalar_grid, gradient_grid,
                       output_scalar_grid, output_gradient_grid);
    break;

  case UNDEFINED_COMMAND:
    error.AddMessage("Programming error.  Undefined merge scalar command.");
    throw error;

  default:
    error.AddMessage("Programming error.  Illegal merge scalar command.");
    throw error;
  }
}

/// Apply scalar operation.
/// @param scalar_operand Scalar operand.  
///        For some operations, scalar_operand is ignored.
void apply_scalar_op
(const SGRID & scalar_grid, const SCALAR_OP_TYPE scalar_op,
 const SCALAR_TYPE scalar_operand, SGRID & output_scalar_grid)
{
  const int dimension = scalar_grid.Dimension();
  IJK::ARRAY<COORD_TYPE> spacing(dimension, 1);
  IJK::PROCEDURE_ERROR error("apply_scalar_op");

  output_scalar_grid.Copy(scalar_grid);

  // *** SHOULD BE AUTOMATIC IN Copy. ***
  output_scalar_grid.SetSpacing(scalar_grid.SpacingPtrConst());


  switch(scalar_op) {

  case NEGATE_SCALAR:
    output_scalar_grid.Multiply(-1);
    break;

  case ADD_SCALAR:
    output_scalar_grid.Add(scalar_operand);
    break;

  case MULTIPLY_SCALAR:
    output_scalar_grid.Multiply(scalar_operand);
    break;

  case SET_SCALAR:
    output_scalar_grid.SetAll(scalar_operand);
    break;

  case SET_BOUNDARY_SCALAR:
    set_boundary_scalar(scalar_operand, output_scalar_grid);
    break;

  case SUBSAMPLE_SCALAR:
    output_scalar_grid.Subsample(scalar_grid, scalar_operandI);
    output_scalar_grid.SetSpacing
      (scalar_operandI, scalar_grid.SpacingPtrConst());
    break;

  case NONUNIFORM_SUBSAMPLE_SCALAR:
    output_scalar_grid.SubsampleNonUniform
      (scalar_grid, IJK::vector2pointer(scalar_operandIV));
    for (int d = 0; d < dimension; d++)
      { spacing[d] = scalar_operandIV[d]*scalar_grid.Spacing(d); }
    output_scalar_grid.SetSpacing(spacing.PtrConst());
    break;

  case SUPERSAMPLE_SCALAR:
    output_scalar_grid.Supersample(scalar_grid, scalar_operandI);
    output_scalar_grid.SetSpacing
      (1.0/scalar_operandI, scalar_grid.SpacingPtrConst());
    break;

  case NONUNIFORM_SUPERSAMPLE_SCALAR:
    output_scalar_grid.SupersampleNonUniform
      (scalar_grid, IJK::vector2pointer(scalar_operandIV));
    for (int d = 0; d < dimension; d++)
      { spacing[d] = scalar_grid.Spacing(d)/scalar_operandIV[d]; }
    output_scalar_grid.SetSpacing(spacing.PtrConst());
    break;

  case UNDEFINED_SCALAR_OP:
    error.AddMessage("Programming error.  Undefined scalar_operation.");
    throw error;

  default:
    error.AddMessage("Programming error.  Illegal scalar operation.");
    throw error;
  }
}

/// Apply scalar operation.
/// @param scalar_operand Scalar operand.  
///        For some operations, scalar_operand is ignored.
void apply_scalar_op
(const SGRID & scalar_grid, const GRADIENT_GRID & gradient_grid,
 const SCALAR_OP_TYPE scalar_op, const SCALAR_TYPE scalar_operand, 
 SGRID & output_scalar_grid, GRADIENT_GRID & output_gradient_grid)
{
  IJK::PROCEDURE_ERROR error("apply_scalar_op");

  output_scalar_grid.Copy(scalar_grid);
  output_gradient_grid.Copy(gradient_grid);

  switch(scalar_op) {

  case NEGATE_SCALAR:
    output_scalar_grid.Multiply(-1);
    output_gradient_grid.ScalarMultiply(-1);
    break;

  case ADD_SCALAR:
    output_scalar_grid.Add(scalar_operand);
    // Adding scalar constant at each vertex does not change the gradient grid.
    break;

  case MULTIPLY_SCALAR:
    output_scalar_grid.Multiply(scalar_operand);
    output_gradient_grid.ScalarMultiply(scalar_operand);
    break;

  case SET_SCALAR:
    output_scalar_grid.SetAll(scalar_operand);
    output_gradient_grid.SetAllCoord(0);   // gradient is 0 everywhere.
    break;

  case SUBSAMPLE_SCALAR:
    error.AddMessage
      ("Error. Scalar operation -subsample not implemented with -grad.");
    throw error;

  case NONUNIFORM_SUBSAMPLE_SCALAR:
    error.AddMessage
      ("Error. Scalar operation -nonuniform_subsample not implemented with -grad.");
    throw error;

  case SUPERSAMPLE_SCALAR:
    error.AddMessage
      ("Error. Scalar operation -supersample not implemented with -grad.");
    throw error;

  case NONUNIFORM_SUPERSAMPLE_SCALAR:
    error.AddMessage
      ("Error. Scalar operation -nonuniform_supersample not implemented with -grad.");
    throw error;

  case UNDEFINED_SCALAR_OP:
    error.AddMessage("Programming error.  Undefined scalar_operation.");
    throw error;

  default:
    error.AddMessage("Programming error.  Illegal scalar operation.");
    throw error;
  }
    
}

void compute_scalar_min
(const SGRID scalar_grid[2], SGRID & output_scalar_grid)
{
  for (VERTEX_INDEX iv = 0; iv < scalar_grid[0].NumVertices(); iv++) {
    SCALAR_TYPE s = 
      std::min(scalar_grid[0].Scalar(iv), scalar_grid[1].Scalar(iv));
    output_scalar_grid.Set(iv, s);
  }
}

void compute_scalar_min
(const SGRID scalar_grid[2], const GRADIENT_GRID gradient_grid[2],
 SGRID & output_scalar_grid, GRADIENT_GRID & output_gradient_grid)
{
  for (VERTEX_INDEX iv = 0; iv < scalar_grid[0].NumVertices(); iv++) {
    SCALAR_TYPE s0 = scalar_grid[0].Scalar(iv);
    SCALAR_TYPE s1 = scalar_grid[1].Scalar(iv);
    if (s0 < s1) {
      output_scalar_grid.Set(iv,s0);
      output_gradient_grid.Set(iv, gradient_grid[0].VectorPtrConst(iv));
    }
    else if (s0 > s1) {
      output_scalar_grid.Set(iv,s1);
      output_gradient_grid.Set(iv, gradient_grid[1].VectorPtrConst(iv));
    }
    else {
      output_scalar_grid.Set(iv,s0);
      for (int ic = 0; ic < gradient_grid[0].VectorLength(); ic++) {
        output_gradient_grid.Set(iv, ic, 0);
      }
    }
  }
}

void compute_scalar_max
(const SGRID scalar_grid[2], SGRID & output_scalar_grid)
{
  for (VERTEX_INDEX iv = 0; iv < scalar_grid[0].NumVertices(); iv++) {
    SCALAR_TYPE s = 
      std::max(scalar_grid[0].Scalar(iv), scalar_grid[1].Scalar(iv));
    output_scalar_grid.Set(iv, s);
  }
}

void compute_scalar_max
(const SGRID scalar_grid[2], const GRADIENT_GRID gradient_grid[2],
 SGRID & output_scalar_grid, GRADIENT_GRID & output_gradient_grid)
{
  for (VERTEX_INDEX iv = 0; iv < scalar_grid[0].NumVertices(); iv++) {
    SCALAR_TYPE s0 = scalar_grid[0].Scalar(iv);
    SCALAR_TYPE s1 = scalar_grid[1].Scalar(iv);
    if (s0 > s1) {
      output_scalar_grid.Set(iv,s0);
      output_gradient_grid.Set(iv, gradient_grid[0].VectorPtrConst(iv));
    }
    else if (s0 < s1) {
      output_scalar_grid.Set(iv,s1);
      output_gradient_grid.Set(iv, gradient_grid[1].VectorPtrConst(iv));
    }
    else {
      output_scalar_grid.Set(iv,s0);
      for (int ic = 0; ic < gradient_grid[0].VectorLength(); ic++) {
        output_gradient_grid.Set(iv, ic, 0);
      }
    }
  }
}


// Set all boundary vertices to val.
void set_boundary_scalar(const SCALAR_TYPE val, SGRID & scalar_grid)
{
  IJK_FOR_EACH_BOUNDARY_GRID_VERTEX(iv, scalar_grid, VERTEX_INDEX) {
    scalar_grid.Set(iv, val);
  }
}

// **************************************************
// Read/write files
// **************************************************

void read_scalar_grid(const char * filename, SGRID & scalar_grid)
{
  NRRD_DATA<int, AXIS_SIZE_TYPE> nrrd_scalar_header;
  std::vector<COORD_TYPE> grid_spacing;
  IJK::PROCEDURE_ERROR error("read_scalar_grid");

  GRID_NRRD_IN<int,AXIS_SIZE_TYPE> nrrd_in;

  nrrd_in.ReadScalarGrid
    (filename, scalar_grid, nrrd_scalar_header, error);
  if (nrrd_in.ReadFailed()) { throw error; }

  nrrd_scalar_header.GetSpacing(grid_spacing);
  scalar_grid.SetSpacing(&(grid_spacing[0]));
}

void read_gradient_grid(const char * filename, GRADIENT_GRID & gradient_grid)
{
  IJK::PROCEDURE_ERROR error("read_gradient_grid");

  GRID_NRRD_IN<int,AXIS_SIZE_TYPE> nrrd_in;
  NRRD_DATA<int, AXIS_SIZE_TYPE> nrrd_header;

  nrrd_in.ReadVectorGrid
    (filename, gradient_grid, nrrd_header, error);
  if (nrrd_in.ReadFailed()) { throw error; }
}

void read_gradient_grid(const string & filename, GRADIENT_GRID & gradient_grid)
{
  read_gradient_grid(filename.c_str(), gradient_grid);
}

void write_scalar_grid
(const string & output_filename, const SGRID & scalar_grid)
{
  const int dimension = scalar_grid.Dimension();
  IJK::ARRAY<double> grid_spacing(dimension, 1);
  NRRD_DATA<int, AXIS_SIZE_TYPE> nrrd_header;

  for (int d = 0; d < dimension; d++) 
    { grid_spacing[d] = scalar_grid.Spacing(d); }

  nrrd_header.SetSize(scalar_grid.Dimension(), scalar_grid.AxisSize());
  nrrdAxisInfoSet_nva(nrrd_header.DataPtr(), nrrdAxisInfoSpacing, 
                      grid_spacing.PtrConst());

  cout << "Writing scalar field to " << output_filename << endl;
  write_scalar_grid_nrrd(output_filename, scalar_grid, nrrd_header);
}

void write_gradient_grid
(const string & output_filename, const GRADIENT_GRID & gradient_grid)
{
  cout << "Writing gradient field to " << output_filename << endl;
  write_vector_grid_nrrd(output_filename, gradient_grid);
}

// **************************************************
// Check Routines
// **************************************************

bool check_gradient_grid
(const GRADIENT_GRID & gradient_grid, IJK::ERROR & error)
{
  if (gradient_grid.VectorLength() != gradient_grid.Dimension()) {
    error.AddMessage
      ("Error in gradient grid (file ", gradient_filename, ").");
    error.AddMessage
      ("  Vector length of gradient grid should equal volume dimension.");
    error.AddMessage
      ("  Volume dimension = ", gradient_grid.Dimension(), ".");
    error.AddMessage
      ("  Vector length = ", gradient_grid.VectorLength(), ".");

    return(false);
  }

  return(true);
}

bool check_input_grids
(const SGRID & scalar_grid, const GRADIENT_GRID & gradient_grid,
 IJK::ERROR & error)
{
  if (!check_gradient_grid(gradient_grid, error)) 
    { return(false); }

  IJK::ERROR size_error;
  if (!gradient_grid.Check
      (scalar_grid, "Gradient grid", "Scalar grid", size_error)) {
    error.AddMessage
      ("Scalar grid (file ", scalar_filename, 
       ") and gradient grid (file ", gradient_filename, ") do not match.");
    for (int i = 0; i < size_error.NumMessages(); i++) 
      { error.AddMessage(size_error.Message(i)); }

    return(false);
  }

  return(true);
}

bool check_scalar_operands
(const SGRID & scalar_grid, IJK::ERROR & error)
{
  if (command == SCALAR_OP) {
    if (scalar_op == NONUNIFORM_SUBSAMPLE_SCALAR) {

      if (scalar_operandIV.size() != scalar_grid.Dimension()) {
        error.AddMessage("Number of values in -nonuniform_subsample does not match grid dimension.");
        error.AddMessage("  Number of values is: ", 
                         scalar_operandIV.size(), ".");
        error.AddMessage("  Grid dimension is: ", scalar_grid.Dimension(), ".");
        return(false);
      }
    }
  }

  return(true);
}


// **************************************************
// Misc Routines
// **************************************************

void memory_exhaustion()
{
  cerr << "Error: Out of memory.  Terminating program." << endl;
  exit(10);
}


IJKSCALAR_COMMAND get_ijkscalar_command(char * s_char)
{
  string s = s_char;
  if (s == "min") { return(MERGE_MIN); }
  else if (s == "max") { return(MERGE_MAX); }
  else if (s == "scalarop") { return(SCALAR_OP); }
  else if (s == "help") { return(HELP); }
  else { return(UNDEFINED_COMMAND); }
}

/// Get argument for option argv[iarg].
/// Does not modify iarg.
template <typename ETYPE>
void get_value
(const int iarg, const int argc, char **argv, ETYPE & x)
{
  if (iarg+1 >= argc) { 
    cerr << "Usage error. Missing argument for option " 
         << argv[iarg] << " and missing file name." << endl;
    usage_error();
  }

  if (!IJK::string2val(argv[iarg+1], x)) {
    cerr << "Error in argument for option: " << argv[iarg] << endl;
    cerr << "Non-numeric character in string: " << argv[iarg+1] << endl;
    exit(50);
  }
}

/// Prints error message if scalar_operand is already set
void check_operand_not_set(const SCALAR_OP_TYPE scalar_operand)
{
  if (scalar_operand != UNDEFINED_SCALAR_OP) {
    cerr << "Usage error.  Only one option -add, -negate, -mult, -set"
         << endl;
    cerr << "    or -set_boundary permitted." << endl;
    usage_error();
  }
}

void parse_command_line(int argc, char **argv)
{
  IJK::PROCEDURE_ERROR error("ijkscalar");

  if (argc < 2) { usage_error(); }

  if (string(argv[1]) == "-help" || string(argv[1]) == "help")
    { help(); }

  if (argv[1][0] == '-') {
    cerr << "Usage error.  First argument to ijkscalar must be a command."
         << endl;
    cerr << "  " << argv[1] << " is not a command." << endl;
    usage_error();
  }

  command = get_ijkscalar_command(argv[1]);

  if (command == HELP) { help(); }

  if (command == UNDEFINED_COMMAND) {
    cerr << "Usage error.  Command \"" << argv[1]
         << "\" not defined." << endl;
    usage_error();
  }

  int iarg = 2;
  while (iarg < argc && argv[iarg][0] == '-') {
    string s = argv[iarg];

    if (s == "-o") {
      iarg++;
      if (iarg >= argc) usage_error();
      output_scalar_filename = argv[iarg];
    }
    else if (s == "-grad") {
      flag_gradients = true;
    }
    else if (s == "-add") {
      check_operand_not_set(scalar_op);
      get_value(iarg, argc, argv, scalar_operand);
      iarg++;
      scalar_op = ADD_SCALAR;
    }
    else if (s == "-negate") {
      check_operand_not_set(scalar_op);
      scalar_op = NEGATE_SCALAR;
    }
    else if (s == "-mult") {
      check_operand_not_set(scalar_op);
      get_value(iarg, argc, argv, scalar_operand);
      iarg++;
      scalar_op = MULTIPLY_SCALAR;
    }
    else if (s == "-set") {
      check_operand_not_set(scalar_op);
      get_value(iarg, argc, argv, scalar_operand);
      iarg++;
      scalar_op = SET_SCALAR;
    }
    else if (s == "-set_boundary") {
      check_operand_not_set(scalar_op);
      get_value(iarg, argc, argv, scalar_operand);
      iarg++;
      scalar_op = SET_BOUNDARY_SCALAR;
    }
    else if (s == "-subsample") {
      check_operand_not_set(scalar_op);
      get_value(iarg, argc, argv, scalar_operandI);
      iarg++;
      scalar_op = SUBSAMPLE_SCALAR;
    }
    else if (s == "-nonuniform_subsample") {
      check_operand_not_set(scalar_op);
      IJK::get_arg_multiple_arguments
        (iarg, argc, argv, scalar_operandIV, error);
      iarg++;
      scalar_op = NONUNIFORM_SUBSAMPLE_SCALAR;
    }
    else if (s == "-supersample") {
      check_operand_not_set(scalar_op);
      get_value(iarg, argc, argv, scalar_operandI);
      iarg++;
      scalar_op = SUPERSAMPLE_SCALAR;
    }
    else if (s == "-nonuniform_supersample") {
      check_operand_not_set(scalar_op);
      IJK::get_arg_multiple_arguments
        (iarg, argc, argv, scalar_operandIV, error);
      iarg++;
      scalar_op = NONUNIFORM_SUPERSAMPLE_SCALAR;
    }

    else if (s == "-help") {
      help();
    }
    else
      { 
        cerr << "Illegal command option: " << s << endl;
        usage_error(); 
      }

    iarg++;
  }

  if (command == SCALAR_OP) {

    if (iarg+1 != argc) {

      if (iarg+1 > argc) {
        cerr << "Usage error.  Missing input filename." << endl;
      }
      else {
        cerr << "Usage error.  Parameter -scalarop uses only one input file."
             << endl;
      }
      usage_error();
    }

    if (scalar_op == UNDEFINED_SCALAR_OP) {
      cerr << "Usage error.  Option -add, -negate, -mult, -set or -set_boundary required"  << endl
           << "  with command scalarop." << endl;
      usage_error();
    }

    scalar_filename[0] = argv[iarg];
    scalar_filename[1] = NULL;
  }
  else {
    if (iarg+2 != argc) {

      if (iarg+2 > argc) {
        cerr << "Usage error.  Missing input filename." << endl;
      }
      else {
        cerr << "Usage error.  Only two input files required." << endl;
      }
      usage_error();
    }

    scalar_filename[0] = argv[iarg];
    scalar_filename[1] = argv[iarg+1];
  }

  if (command == SCALAR_OP) {

    if (scalar_op == SUBSAMPLE_SCALAR || scalar_op == SUPERSAMPLE_SCALAR) {
      if (scalar_operandI < 1) {
        cerr << "Input error. Parameter for -subsample must be an integer greater than 0."
             << endl;
        usage_error();
      }
    }
  }

}

void parse_filename
(const string & filename, string & prefix, string & suffix)
{
  // remove path from file name
  string prefix2, suffix2;
  split_string(filename, '/', prefix2, suffix2);
  if (suffix2 == "") 
    { split_string(filename, '.', prefix, suffix); }
  else 
    { split_string(suffix2, '.', prefix, suffix); }
  
  if (suffix != "nrrd" && suffix != "nhdr") {
    prefix = filename;
    suffix = "";
  }
}

void split_string(const string & s, const char c,
                  string & prefix, string & suffix)
// split string at last occurrence of character c into prefix and suffix
{
  string::size_type i = s.rfind(c);
  if (i == string::npos) {
    prefix = s;
    suffix = "";
  }
  else {
    if (i > 0) { prefix = s.substr(0,i); }
    else { prefix = ""; };

    if (i+1 < s.length()) { suffix = s.substr(i+1, s.length()-i-1); }
    else { suffix = ""; };
  }
}

void get_gradient_filename
(const char * scalar_filename, const char * gradient_filename, 
 string & gradient_filename_str)
{
  if (gradient_filename == NULL) {

    if (scalar_filename == NULL) {
      gradient_filename == "";
    }
    else {

      string prefix, suffix;

      parse_filename(string(scalar_filename), prefix, suffix);
      gradient_filename_str = prefix + ".grad." + suffix;
    }
  }
  else {
    gradient_filename_str = gradient_filename;
  }
}

void get_gradient_filename
(const string & scalar_filename, const char * gradient_filename, 
 string & gradient_filename_str)
{
  get_gradient_filename(scalar_filename.c_str(), gradient_filename,
                        gradient_filename_str);
}

void usage_error()
{
  cerr << "Usage: ijkscalar COMMAND [OPTIONS] <infile1> [<infile2>]"
       << endl;
  cerr << "COMMAND:" << endl;
  cerr << "  min || max || scalarop || help" << endl;
  cerr << "OPTIONS:" << endl;
  cerr << "  [-grad] [-add <x>|-negate|-mult <x>|-set <x>|-set_boundary<x>|" 
       << endl;
  cerr << "   -subsample <n>|-nonuniform_subsample \"<n0 n1...>\"|" << endl;
  cerr << "   -supersample <n>|-nonuniform_supersample \"<n0 n1...>\"]"
       << endl;
  cerr << "  [-o <output filename>]" << endl;
  exit(10);
}

void help()
{
  cerr << "Usage: ijkscalar COMMAND [OPTIONS] <infile1> [<infile2>] <outfile>"
       << endl;
  cerr << endl;
  cerr << "COMMAND:" << endl;
  cerr << "  min:  Compute pointwise minimum of two scalar grids." << endl
       << "        Requires two input files." << endl;
  cerr << "  max:  Compute pointwise maximum of two scalar grids." << endl
       << "        Requires two input files." << endl;
  cerr << "  scalarop: Apply a scalar operation to the scalar field." << endl
       << "        Requires one input file." << endl;
  cerr << "  help: Print this help message." << endl;
  cerr << endl;
  cerr << "OPTIONS:" << endl;
  cerr << "  -grad: Compute and output gradient file for computed scalar grid."
       << endl;
  cerr << "  -add <x>:  Add x to each scalar value." << endl;
  cerr << "  -negate:   Multiply each scalar value by -1." << endl;
  cerr << "  -mult <x>: Multiply each scalar value by x." << endl;
  cerr << "  -set <x>:  Set each scalar value to x." << endl;
  cerr << "  -set_boundary <x>: Set each boundary scalar value to x." << endl;
  cerr << "  -subsample <n>: Uniform subsample with period n." << endl;
  cerr << "  -nonuniform_subsample \"<n0 n1...>\": Non-uniform subsampling."
       << endl;
  cerr << "       Period for axis d is nd." << endl;
  cerr << "  -supersample <n>: Uniform supersample with period n." << endl;
  cerr << "  -nonuniform_supersample \"<n0 n1...>\": Non-uniform supersampling."
       << endl;
  cerr << "       Period for axis d is nd." << endl;
  cerr << "  -o <output filename>: Filename of output scalar file." << endl
       << "       Filename of gradient file is derived from scalar output filename."
       << endl;
  exit(10);
}
