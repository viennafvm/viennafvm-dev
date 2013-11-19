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

#define VIENNAFVM_DEBUG

// include necessary system headers
#include <iostream>

// ViennaFVM includes:
#include "viennafvm/forwards.h"
#include "viennafvm/linear_assembler.hpp"
#include "viennafvm/io/vtk_writer.hpp"
#include "viennafvm/boundary.hpp"
#include "viennafvm/pde_solver.hpp"
#include "viennafvm/linear_solvers/viennacl.hpp"
#include "viennafvm/problem_description.hpp"


// ViennaGrid includes:
#include "viennagrid/mesh/mesh.hpp"
#include <viennagrid/config/default_configs.hpp>
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

  typedef viennagrid::triangular_2d_mesh                        MeshType;
  typedef viennagrid::result_of::segmentation<MeshType>::type   SegmentationType;

  typedef viennagrid::result_of::cell_tag<MeshType>::type CellTag;

  typedef viennagrid::result_of::element<MeshType, CellTag>::type        CellType;

  typedef viennagrid::result_of::element_range<MeshType, CellTag>::type       CellContainer;
  typedef viennagrid::result_of::iterator<CellContainer>::type                CellIterator;
  typedef viennagrid::result_of::vertex_range<CellType>::type                 VertexOnCellContainer;
  typedef viennagrid::result_of::iterator<VertexOnCellContainer>::type        VertexOnCellIterator;

  typedef viennamath::function_symbol   FunctionSymbol;
  typedef viennamath::equation          Equation;

  //
  // Create a mesh from file
  //
  MeshType mesh;
  SegmentationType segmentation(mesh);

  try
  {
    viennagrid::io::netgen_reader my_reader;
    my_reader(mesh, segmentation, "../examples/data/square128.mesh");
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
  QuantityType & quan_u = problem_desc.add_quantity("u");
  QuantityType & quan_v = problem_desc.add_quantity("v");

  //
  // Setting boundary information on mesh (this should come from device specification)
  //

  viennagrid::result_of::default_point_accessor<MeshType>::type point_accessor = viennagrid::default_point_accessor(mesh);

  CellContainer cells(mesh);
  for (CellIterator cit  = cells.begin();
                    cit != cells.end();
                  ++cit)
  {
    VertexOnCellContainer vertices(*cit);
    for (VertexOnCellIterator vit  = vertices.begin();
                              vit != vertices.end();
                            ++vit)
    {
      //boundary for first equation: Homogeneous Dirichlet everywhere
      if ( point_accessor(*vit)[0] == 0.0 || point_accessor(*vit)[0] == 1.0
           || point_accessor(*vit)[1] == 0.0 || point_accessor(*vit)[1] == 1.0 )
      {
        //simulation with ID 0 uses homogeneous boundary data
        quan_u.set_boundary_type(*cit, viennafvm::BOUNDARY_DIRICHLET);
        quan_u.set_boundary_value(*cit, 0.0);
      }
      else
        quan_u.set_unknown_mask(*cit, true);

      //boundary for second equation (ID 1): 0 at left boundary, 1 at right boundary
      if ( point_accessor(*vit)[0] == 0.0)
      {
        quan_v.set_boundary_type(*cit, viennafvm::BOUNDARY_DIRICHLET);
        quan_v.set_boundary_value(*cit, 0.0);
      }
      else if ( point_accessor(*vit)[0] == 1.0)
      {
        quan_v.set_boundary_type(*cit, viennafvm::BOUNDARY_DIRICHLET);
        quan_v.set_boundary_value(*cit, 1.0);
      }
      else
        quan_v.set_unknown_mask(*cit, true);
    }
  }


  //
  // Specify two PDEs:
  //
  FunctionSymbol u(0, viennamath::unknown_tag<>());   //an unknown function used for PDE specification
  FunctionSymbol v(1, viennamath::unknown_tag<>());   //an unknown function used for PDE specification
  Equation poisson_equ_1 = viennamath::make_equation( viennamath::laplace(u), -1);
  Equation poisson_equ_2 = viennamath::make_equation( viennamath::laplace(v), 0);

  //
  // Setup Linear Solver
  //
  viennafvm::linsolv::viennacl  linear_solver;

  //
  // Pass system to solver:
  //
  viennafvm::pde_solver my_solver;
  my_solver(problem_desc,
            viennafvm::make_linear_pde_system(poisson_equ_1, u),  // PDE with associated unknown
            linear_solver);

  my_solver(problem_desc,
            viennafvm::make_linear_pde_system(poisson_equ_2, v, viennafvm::make_linear_pde_options(1, 1)),  // PDE with associated unknown
             linear_solver);


  //
  // Writing solution back to mesh
  //
  viennafvm::io::write_solution_to_VTK_file(problem_desc.quantities(), "poisson_2d", mesh, segmentation);

  std::cout << "*****************************************" << std::endl;
  std::cout << "* Poisson solver finished successfully! *" << std::endl;
  std::cout << "*****************************************" << std::endl;
  return EXIT_SUCCESS;
}

