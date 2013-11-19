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
  typedef double   numeric_type;

  typedef viennagrid::tetrahedral_3d_mesh   MeshType;
  typedef viennagrid::result_of::segmentation<MeshType>::type SegmentationType;

  typedef viennagrid::result_of::cell_tag<MeshType>::type CellTag;

  typedef viennagrid::result_of::element<MeshType, CellTag>::type        CellType;

  typedef viennagrid::result_of::element_range<MeshType, CellTag>::type  CellContainer;
  typedef viennagrid::result_of::iterator<CellContainer>::type                CellIterator;
  typedef viennagrid::result_of::vertex_range<CellType>::type               VertexOnCellContainer;
  typedef viennagrid::result_of::iterator<VertexOnCellContainer>::type        VertexOnCellIterator;

  typedef viennamath::function_symbol   FunctionSymbol;
  typedef viennamath::equation          Equation;

  typedef viennafvm::boundary_key      BoundaryKey;

  //
  // Create a mesh from file
  //
  MeshType mesh;
  SegmentationType segmentation(mesh);

  try
  {
    viennagrid::io::netgen_reader my_netgen_reader;
    my_netgen_reader(mesh, segmentation, "../examples/data/cube3072.mesh");
  }
  catch (...)
  {
    std::cerr << "File-Reader failed. Aborting program..." << std::endl;
    return EXIT_FAILURE;
  }

  //
  // Create PDE solver instance:
  //
  viennafvm::problem_description<MeshType> problem_desc(mesh);

  // Initialize unknowns:
  typedef viennafvm::problem_description<MeshType>::quantity_type   QuantityType;
  QuantityType & quan = problem_desc.add_quantity("u");

  //
  // Setting boundary information on mesh (this should come from device specification)
  //
  //setting some boundary flags:
  viennagrid::result_of::default_point_accessor<MeshType>::type point_accessor = viennagrid::default_point_accessor(mesh);

  CellContainer cells(mesh);
  for (CellIterator cit  = cells.begin();
                    cit != cells.end();
                  ++cit)
  {
    bool cell_on_boundary = false;

    VertexOnCellContainer vertices(*cit);
    for (VertexOnCellIterator vit  = vertices.begin();
                              vit != vertices.end();
                            ++vit)
    {
      //boundary for first equation: Homogeneous Dirichlet everywhere
      if (point_accessor(*vit)[0] == 0.0 || point_accessor(*vit)[0] == 1.0
        || point_accessor(*vit)[2] == 0.0 || point_accessor(*vit)[2] == 1.0 )
      {
        cell_on_boundary = true;
        break;
      }
    }
    if (cell_on_boundary)
    {
      quan.set_boundary_type(*cit, viennafvm::BOUNDARY_DIRICHLET);
      quan.set_boundary_value(*cit, 0.0);
    }
    else
      quan.set_unknown_mask(*cit, true);
  }


  problem_desc.add_quantity("permittivity");

  // Specify Poisson equation:
  FunctionSymbol u(0, viennamath::unknown_tag<>());   //an unknown function used for PDE specification
  FunctionSymbol permittivity(1);
  Equation poisson_eq = viennamath::make_equation( viennamath::div(permittivity * viennamath::grad(u)), -1);  // \Delta u = -1

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
  viennafvm::io::write_solution_to_VTK_file(problem_desc.quantities(), "poisson_3d", mesh, segmentation);

  std::cout << "*****************************************" << std::endl;
  std::cout << "* Poisson solver finished successfully! *" << std::endl;
  std::cout << "*****************************************" << std::endl;
  return EXIT_SUCCESS;
}
