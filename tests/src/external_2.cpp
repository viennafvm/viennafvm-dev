/* =======================================================================
   Copyright (c) 2013, Institute for Microelectronics,
                       Institute for Analysis and Scientific Computing,
                       TU Wien.

                            -----------------
                ViennaFVM - The Vienna Finite Volume Library
                            -----------------

   Authors:      Karl Rupp                           rupp@iue.tuwien.ac.at
                 Josef Weinbub                    weinbub@iue.tuwien.ac.at

   (A list of additional contributors can be found in the PDF manual)

   License:      MIT (X11), see file LICENSE in the base directory
======================================================================= */

#ifdef _MSC_VER      //Visual Studio complains about potentially dangerous things, which are perfectly legal in our context
  #pragma warning( disable : 4355 )     //use of this in member initializer list
  #pragma warning( disable : 4503 )     //truncated name decoration
#endif

//
// A check for the absence of external linkage (otherwise, library is not truly 'header-only')
//

#include <iostream>

#include "viennagrid/mesh/mesh.hpp"
#include "viennagrid/config/default_configs.hpp"

//
// *** ViennaFVM
//

#include "viennafvm/forwards.h"
#include "viennafvm/linear_assembler.hpp"
#include "viennafvm/io/vtk_writer.hpp"
#include "viennafvm/boundary.hpp"
#include "viennafvm/pde_solver.hpp"
#include "viennafvm/initial_guess.hpp"
#include "viennafvm/linear_solvers/viennacl.hpp"
#include "viennafvm/problem_description.hpp"


void other_func()
{
  typedef viennagrid::triangular_2d_mesh  MeshType;

  MeshType mesh;

  viennafvm::problem_description<MeshType> problem_desc(mesh);
  std::cout << "other_func(): Number of quantities in pde_solver per default: " << problem_desc.quantities().size() << std::endl;
}
