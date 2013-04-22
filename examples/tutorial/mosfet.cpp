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

// Define NDEBUG to get any reasonable performance with ublas:
#define NDEBUG

// include necessary system headers
#include <iostream>

// ViennaFVM includes:
#include "viennafvm/forwards.h"
#include "viennafvm/linear_assembler.hpp"
#include "viennafvm/io/vtk_writer.hpp"
#include "viennafvm/boundary.hpp"
#include "viennafvm/pde_solver.hpp"

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

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>



struct permittivity_key
{
  // Operator< is required for compatibility with std::map
  bool operator<(permittivity_key const & other) const { return false; }
};

struct net_doping_key
{
  // Operator< is required for compatibility with std::map
  bool operator<(net_doping_key const & other) const { return false; }
};


int main()
{
  typedef double   numeric_type;

  typedef viennagrid::config::triangular_2d                           ConfigType;
  typedef viennagrid::result_of::domain<ConfigType>::type             DomainType;
  typedef typename ConfigType::cell_tag                     CellTag;

  typedef viennagrid::result_of::ncell<ConfigType, CellTag::dim>::type        CellType;

  typedef viennamath::function_symbol   FunctionSymbol;
  typedef viennamath::equation          Equation;

  //
  // Create a domain from file
  //
  DomainType my_domain;

  try
  {
    viennagrid::io::netgen_reader my_reader;
    my_reader(my_domain, "../examples/data/mosfet.mesh");
  }
  catch (...)
  {
    std::cerr << "File-Reader failed. Aborting program..." << std::endl;
    return EXIT_FAILURE;
  }

  //
  // Initialize cell_quantity
  //
  viennafvm::ncell_quantity<CellType, viennamath::expr::interface_type>  permittivity; permittivity.wrap_constant( permittivity_key() );
  viennafvm::set_quantity_region(permittivity_key(), my_domain, true);               // permittivity is (for simplicity) defined everywhere
  viennafvm::set_quantity_value(permittivity_key(), my_domain, 11.7);                // relative permittivity of silicon
  viennafvm::set_quantity_value(permittivity_key(), my_domain.segments()[2], 15.6);  // relative permittivty of HfO2

  viennafvm::ncell_quantity<CellType, viennamath::expr::interface_type>  net_doping; net_doping.wrap_constant( net_doping_key() );

  //
  // Specify PDEs:
  //
  FunctionSymbol psi(0);   // potential, using id=0
  FunctionSymbol n(1);     // electron concentration, using id=1
  FunctionSymbol p(2);     // hole concentration, using id=2

  double q  = 1.6e-19;
  double kB = 1.38e-23; // Boltzmann constant
  double mu = 1;        // mobility (constant is fine for the moment)
  double T  = 300;
  double D  = mu * kB * T / q;  //diffusion constant

  // here is all the fun: specify DD system
  Equation poisson_eq = viennamath::make_equation( viennamath::div(permittivity * viennamath::grad(psi)),                     /* = */ q * (n - p - net_doping));
  Equation cont_eq_n  = viennamath::make_equation( viennamath::div(D * viennamath::grad(n) + mu * viennamath::grad(psi) * n), /* = */ 0);
  Equation cont_eq_p  = viennamath::make_equation( viennamath::div(D * viennamath::grad(p) - mu * viennamath::grad(psi) * p), /* = */ 0);

  //
  // Setting boundary information on domain (see mosfet.in2d for segment indices)
  //

  // potential:
  viennafvm::set_dirichlet_boundary(my_domain.segments()[0], 0.5, psi); // Gate contact
  viennafvm::set_dirichlet_boundary(my_domain.segments()[1], 0.0, psi); // Source contact
  viennafvm::set_dirichlet_boundary(my_domain.segments()[3], 0.5, psi); // Drain contact
  // using floating body, hence commented: viennafvm::set_dirichlet_boundary(my_domain.segments()[7], 0.0, psi); // Body contact

  // electron density
  viennafvm::set_dirichlet_boundary(my_domain.segments()[1], 1e20, n); // Source contact
  viennafvm::set_dirichlet_boundary(my_domain.segments()[3], 1e20, n); // Drain contact

  // hole density
  viennafvm::set_dirichlet_boundary(my_domain.segments()[1], 1e12, p); // Source contact
  viennafvm::set_dirichlet_boundary(my_domain.segments()[3], 1e12, p); // Drain contact


  //
  // Set quantity mask: By default, a quantity is defined on the entire domain.
  //                    Boundary equations are already specified above.
  //                    All that is left is to specify regions where a quantity 'does not make sense'
  //                    Here, we need to disable {n,p} in the gate oxide and the gate
  //

  viennafvm::disable_quantity(my_domain.segments()[0], n); // Gate contact
  viennafvm::disable_quantity(my_domain.segments()[2], n); // Gate oxide
  // add body contact if not using floating body

  viennafvm::disable_quantity(my_domain.segments()[0], p); // Gate contact
  viennafvm::disable_quantity(my_domain.segments()[2], p); // Gate oxide
  // add body contact if not using floating body


  //
  // TODO: Set initial guess and doping
  //



  //
  // Create PDE solver instance:
  //
  viennafvm::pde_solver<> pde_solver;


  //
  // Solve system and write solution vector to pde_result:
  // (discussion about proper interface required. Introduce a pde_result class?)
  //
  viennafvm::linear_pde_system<> pde_system;
  pde_system.add_pde(poisson_eq, psi); // equation and associated quantity
  //pde_system.add_pde(cont_eq_n, n);    // equation and associated quantity
  //pde_system.add_pde(cont_eq_p, p);    // equation and associated quantity

  pde_solver(pde_system, my_domain);   // weird math happening in here ;-)


  //
  // Writing all solution variables back to domain
  //
  std::vector<long> result_ids(pde_system.size());
  for (std::size_t i=0; i<pde_system.size(); ++i)
    result_ids[i] = pde_system.unknown(i)[0].id();

  viennafvm::io::write_solution_to_VTK_file(pde_solver.result(), "mosfet", my_domain, result_ids);

  std::cout << "********************************************" << std::endl;
  std::cout << "* MOSFET simulation finished successfully! *" << std::endl;
  std::cout << "********************************************" << std::endl;
  return EXIT_SUCCESS;
}

