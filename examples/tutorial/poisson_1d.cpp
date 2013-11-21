/* =======================================================================
   Copyright (c) 2011, Institute for Microelectronics, TU Wien
   http://www.iue.tuwien.ac.at
                             -----------------
           ViennaFVM - The Vienna Finite Volume Method Library
                             -----------------

   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at
               (add your name here)

   license:    To be discussed, see file LICENSE in the ViennaFVM base directory
======================================================================= */

// include necessary system headers
#include <iostream>

#define VIENNAFVM_DEBUG

// ViennaFVM includes:
#include "viennafvm/forwards.h"
//#include "viennafvm/poisson_assembler.hpp"
#include "viennafvm/linear_assembler.hpp"
#include "viennafvm/io/vtk_writer.hpp"
#include "viennafvm/ncell_quantity.hpp"
#include "viennafvm/pde_solver.hpp"
#include "viennafvm/boundary.hpp"
#include "viennafvm/linear_solvers/viennacl.hpp"
#include "viennafvm/problem_description.hpp"

// ViennaGrid includes:
#include "viennagrid/mesh/mesh.hpp"
#include "viennagrid/config/default_configs.hpp"
#include "viennagrid/io/netgen_reader.hpp"
#include "viennagrid/io/vtk_writer.hpp"
#include "viennagrid/algorithm/voronoi.hpp"

// ViennaData includes:
#include "viennadata/api.hpp"

// ViennaMath includes:
#include "viennamath/expression.hpp"


int main()
{
  typedef viennagrid::line_1d_mesh   MeshType;
  typedef viennagrid::result_of::segmentation<MeshType>::type SegmentationType;

  typedef viennamath::function_symbol   FunctionSymbol;
  typedef viennamath::equation          Equation;

  //
  // Create a mesh from file
  //
  MeshType mesh;
  SegmentationType segmentation(mesh);

  try
  {
    viennagrid::io::netgen_reader my_netgen_reader;
    my_netgen_reader(mesh, segmentation, "../examples/data/line23.mesh");
  }
  catch (...)
  {
    std::cerr << "File-Reader failed. Aborting program..." << std::endl;
    return EXIT_FAILURE;
  }

  //
  // Create problem description
  //
  viennafvm::problem_description<MeshType> problem_desc(mesh);

  // Initialize unknowns:
  typedef viennafvm::problem_description<MeshType>::quantity_type   QuantityType;
  QuantityType & quan = problem_desc.add_quantity("u");
  viennafvm::set_dirichlet_boundary(quan, segmentation(1), 0.0);
  viennafvm::set_dirichlet_boundary(quan, segmentation(5), 1.0);

  viennafvm::set_unknown(quan, segmentation(2));
  viennafvm::set_unknown(quan, segmentation(3));
  viennafvm::set_unknown(quan, segmentation(4));

  QuantityType & quan_permittivity = problem_desc.add_quantity("permittivity", 3);
  viennafvm::set_initial_value(quan_permittivity, segmentation(3), 1.0);
  viennafvm::set_initial_value(quan_permittivity, segmentation(4), 1.0);
  viennafvm::set_initial_value(quan_permittivity, segmentation(5), 1.0);


  // Specify Poisson equation:
  FunctionSymbol u(0, viennamath::unknown_tag<>());   //an unknown function used for PDE specification
  FunctionSymbol permittivity(1);
  Equation poisson_eq = viennamath::make_equation( viennamath::div(permittivity * viennamath::grad(u)), 0);  // \Delta u = 0

  //
  // Setup Linear Solver
  //
  viennafvm::linsolv::viennacl  linear_solver;

  //
  // Pass system to solver:
  //
  viennafvm::pde_solver my_solver;
  my_solver(problem_desc,
            viennafvm::make_linear_pde_system(poisson_eq, u),  // PDE with associated unknown
            linear_solver);

  //
  // Writing solution back to mesh (discussion about proper way of returning a solution required...)
  //
  viennafvm::io::write_solution_to_VTK_file(problem_desc.quantities(), "poisson_1d", mesh, segmentation);

  std::cout << "*****************************************" << std::endl;
  std::cout << "* Poisson solver finished successfully! *" << std::endl;
  std::cout << "*****************************************" << std::endl;
  return EXIT_SUCCESS;
}
