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
  typedef double   numeric_type;

  typedef viennagrid::line_1d_mesh   MeshType;
  typedef viennagrid::result_of::segmentation<MeshType>::type SegmentationType;

  typedef viennagrid::result_of::cell_tag<MeshType>::type CellTag;

  typedef viennagrid::result_of::element<MeshType, CellTag>::type        CellType;

  typedef viennamath::function_symbol   FunctionSymbol;
  typedef viennamath::equation          Equation;

  typedef viennadata::storage<> StorageType;

  typedef viennafvm::boundary_key      BoundaryKey;

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
  // Create PDE solver instance
  //
  viennafvm::pde_solver<MeshType> pde_solver(mesh);

  // Initialize unknowns:
  typedef viennafvm::pde_solver<MeshType>::quantity_type   QuantityType;
  QuantityType quan("u", cells(mesh).size());
  viennafvm::set_dirichlet_boundary(quan, segmentation(1), 0.0);
  viennafvm::set_dirichlet_boundary(quan, segmentation(5), 1.0);

  viennafvm::set_unknown(quan, segmentation(2));
  viennafvm::set_unknown(quan, segmentation(3));
  viennafvm::set_unknown(quan, segmentation(4));

  pde_solver.add_quantity(quan);

  QuantityType quan_permittivity("permittivity", cells(mesh).size(), 3);
  viennafvm::set_initial_value(quan_permittivity, segmentation(3), 1.0);
  viennafvm::set_initial_value(quan_permittivity, segmentation(4), 1.0);
  viennafvm::set_initial_value(quan_permittivity, segmentation(5), 1.0);
  pde_solver.add_quantity(quan_permittivity);


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
  pde_solver(viennafvm::make_linear_pde_system(poisson_eq, u),  // PDE with associated unknown
             linear_solver);

  //
  // Writing solution back to mesh (discussion about proper way of returning a solution required...)
  //
  viennafvm::io::write_solution_to_VTK_file(pde_solver.quantities(), "poisson_1d", mesh, segmentation);

  std::cout << "*****************************************" << std::endl;
  std::cout << "* Poisson solver finished successfully! *" << std::endl;
  std::cout << "*****************************************" << std::endl;
  return EXIT_SUCCESS;
}
