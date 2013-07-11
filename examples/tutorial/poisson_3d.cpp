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
#include "viennafvm/linear_solve.hpp"
#include "viennafvm/pde_solver.hpp"

// ViennaGrid includes:
#include "viennagrid/domain/domain.hpp"
#include <viennagrid/config/default_configs.hpp>
#include "viennagrid/io/netgen_reader.hpp"
#include "viennagrid/io/vtk_writer.hpp"
#include "viennagrid/algorithm/voronoi.hpp"

// ViennaData includes:
#include "viennadata/api.hpp"

// ViennaMath includes:
#include "viennamath/expression.hpp"

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>


struct permittivity_key
{
  // Operator< is required for compatibility with std::map
  bool operator<(permittivity_key const & /*other*/) const { return false; }
};


int main()
{
  typedef double   numeric_type;

  typedef viennagrid::tetrahedral_3d_domain   DomainType;
  typedef viennagrid::result_of::segmentation<DomainType>::type SegmentationType;

  typedef viennagrid::result_of::cell_tag<DomainType>::type CellTag;

  typedef viennagrid::result_of::element<DomainType, CellTag>::type        CellType;

  typedef viennamath::function_symbol   FunctionSymbol;
  typedef viennamath::equation          Equation;

  typedef viennadata::storage<> StorageType;


  typedef viennagrid::result_of::element_range<DomainType, CellTag>::type  CellContainer;
  typedef viennagrid::result_of::iterator<CellContainer>::type                CellIterator;
  typedef viennagrid::result_of::vertex_range<CellType>::type               VertexOnCellContainer;
  typedef viennagrid::result_of::iterator<VertexOnCellContainer>::type        VertexOnCellIterator;

  typedef viennamath::function_symbol   FunctionSymbol;
  typedef viennamath::equation          Equation;

  typedef viennafvm::boundary_key      BoundaryKey;

  //
  // Create a domain from file
  //
  DomainType my_domain;
  SegmentationType my_segmentation(my_domain);
  StorageType my_storage;

  try
  {
    viennagrid::io::netgen_reader my_netgen_reader;
    my_netgen_reader(my_domain, my_segmentation, "../examples/data/cube3072.mesh");
  }
  catch (...)
  {
    std::cerr << "File-Reader failed. Aborting program..." << std::endl;
    return EXIT_FAILURE;
  }

  // Specify Poisson equation:
  viennafvm::ncell_quantity<CellType, viennamath::expr::interface_type>  permittivity; permittivity.wrap_constant( my_storage, permittivity_key() );

  FunctionSymbol u(0, viennamath::unknown_tag<>());   //an unknown function used for PDE specification
  Equation poisson_eq = viennamath::make_equation( viennamath::div(permittivity * viennamath::grad(u)), -1);  // \Delta u = -1

  //
  // Setting boundary information on domain (this should come from device specification)
  //
  //setting some boundary flags:

  viennagrid::result_of::default_point_accessor<DomainType>::type point_accessor = viennagrid::default_point_accessor(my_domain);
  
  viennadata::result_of::accessor<StorageType, permittivity_key, double, CellType>::type permittivity_accessor =
      viennadata::accessor<permittivity_key, double, CellType>(my_storage, permittivity_key());

  viennadata::result_of::accessor<StorageType, BoundaryKey, bool, CellType>::type boundary_accessor =
      viennadata::accessor<BoundaryKey, bool, CellType>(my_storage, BoundaryKey(0));
  
  CellContainer cells = viennagrid::elements(my_domain);
  for (CellIterator cit  = cells.begin();
                    cit != cells.end();
                  ++cit)
  {
    bool cell_on_boundary = false;

    // write dummy permittivity to cells:
    permittivity_accessor(*cit) = 1.0;

    VertexOnCellContainer vertices = viennagrid::elements(*cit);
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
    boundary_accessor(*cit) = cell_on_boundary;
  }


  //
  // Create PDE solver instance
  //
  viennafvm::pde_solver<> pde_solver;

  //
  // Pass system to solver:
  //
  pde_solver(my_storage,
             viennafvm::make_linear_pde_system(poisson_eq, u),  // PDE with associated unknown
             my_domain);

  //
  // Writing solution back to domain (discussion about proper way of returning a solution required...)
  //
  viennafvm::io::write_solution_to_VTK_file(my_storage, pde_solver.result(), "poisson_3d", my_domain, my_segmentation, 0);

  std::cout << "*****************************************" << std::endl;
  std::cout << "* Poisson solver finished successfully! *" << std::endl;
  std::cout << "*****************************************" << std::endl;
  return EXIT_SUCCESS;
}
