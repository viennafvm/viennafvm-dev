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
#include "viennafvm/boundary.hpp"

// ViennaGrid includes:
#include "viennagrid/domain.hpp"
#include <viennagrid/config/simplex.hpp>
#include "viennagrid/io/netgen_reader.hpp"
#include "viennagrid/io/vtk_writer.hpp"
#include "viennagrid/algorithm/voronoi.hpp"

// ViennaData includes:
#include "viennadata/api.hpp"

// ViennaMath includes:
#include "viennamath/expression.hpp"


struct permittivity_key
{
  // Operator< is required for compatibility with std::map
  bool operator<(permittivity_key const & other) const { return false; }
};


int main()
{
  typedef double   numeric_type;

  typedef viennagrid::config::line_1d                       ConfigType;
  typedef viennagrid::result_of::domain<ConfigType>::type   DomainType;
  typedef typename ConfigType::cell_tag                     CellTag;

  typedef viennagrid::result_of::ncell<ConfigType, CellTag::dim>::type        CellType;

  typedef viennamath::function_symbol   FunctionSymbol;
  typedef viennamath::equation          Equation;

  typedef viennafvm::boundary_key      BoundaryKey;

  //
  // Create a domain from file
  //
  DomainType my_domain;

  try
  {
    viennagrid::io::netgen_reader my_netgen_reader;
    my_netgen_reader(my_domain, "../examples/data/line10.mesh");
  }
  catch (...)
  {
    std::cerr << "File-Reader failed. Aborting program..." << std::endl;
    return EXIT_FAILURE;
  }

  // Specify Poisson equation:
  viennafvm::ncell_quantity<CellType, viennamath::expr::interface_type>  permittivity; permittivity.wrap_constant( permittivity_key() );
  std::cout << "No. of segments: " << my_domain.segments().size() << std::endl;
  viennafvm::set_quantity_region(permittivity_key(), my_domain, true);
  viennafvm::set_quantity_value(permittivity_key(), my_domain, 1);
  viennafvm::set_quantity_value(permittivity_key(), my_domain.segments()[2], 2);
  viennafvm::set_quantity_value(permittivity_key(), my_domain.segments()[3], 2);
  viennafvm::set_quantity_value(permittivity_key(), my_domain.segments()[4], 2);

  FunctionSymbol u(0, viennamath::unknown_tag<>());   //an unknown function used for PDE specification
  Equation poisson_eq = viennamath::make_equation( viennamath::div(permittivity * viennamath::grad(u)), 0);  // \Delta u = 0

  //
  // Setting boundary information on domain (this should come from device specification)
  //
  viennafvm::set_dirichlet_boundary(my_domain.segments()[0], 0.0, 0);
  viennafvm::set_dirichlet_boundary(my_domain.segments()[4], 1.0, 0);


  //
  // Create PDE solver instance
  //
  viennafvm::pde_solver<> pde_solver;

  //
  // Pass system to solver:
  //
  pde_solver(viennafvm::make_linear_pde_system(poisson_eq, u),  // PDE with associated unknown
             my_domain);

  //
  // Writing solution back to domain (discussion about proper way of returning a solution required...)
  //
  viennafvm::io::write_solution_to_VTK_file(pde_solver.result(), "poisson_1d", my_domain, 0);

  std::cout << "*****************************************" << std::endl;
  std::cout << "* Poisson solver finished successfully! *" << std::endl;
  std::cout << "*****************************************" << std::endl;
  return EXIT_SUCCESS;
}
