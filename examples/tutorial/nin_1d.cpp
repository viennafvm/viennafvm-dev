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

//#define VIENNAFVM_DEBUG

// Define NDEBUG to get any reasonable performance with ublas:
//#define NDEBUG

// include necessary system headers
#include <iostream>

// ViennaFVM includes:
#include "viennafvm/forwards.h"
#include "viennafvm/linear_assembler.hpp"
#include "viennafvm/io/vtk_writer.hpp"
#include "viennafvm/boundary.hpp"
#include "viennafvm/pde_solver.hpp"
#include "viennafvm/initial_guess.hpp"

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

//
// Defining a bunch of accessor keys for physical quantities.
// These should be part of a semiconductor application based on ViennaFVM, not of ViennaFVM itself
// (which is just a generic finite volume solver and agnostic with respect to the actual physics).
//

struct permittivity_key
{
  // Operator< is required for compatibility with std::map
  bool operator<(permittivity_key const & /*other*/) const { return false; }
};

struct builtin_potential_key
{
  // Operator< is required for compatibility with std::map
  bool operator<(builtin_potential_key const & /*other*/) const { return false; }
};

// N_D
struct donator_doping_key
{
  // Operator< is required for compatibility with std::map
  bool operator<(donator_doping_key const & /*other*/) const { return false; }
};

// N_A
struct acceptor_doping_key
{
  // Operator< is required for compatibility with std::map
  bool operator<(acceptor_doping_key const & /*other*/) const { return false; }
};





double built_in_potential(double /*temperature*/, double doping_n, double doping_p)
{
  const double net_doping = doping_n - doping_p;
  const double x = std::abs(net_doping) / (2.0 * 1e16);

  double bpot = 0.026 * std::log(x + std::sqrt( 1.0 + x*x ) );
                              // V_T * arsinh( net_doping/(2 n_i))

  if ( net_doping < 0) //above formula does not yet consider the doping polarity
    bpot *= -1.0;

  return bpot;
}

template <typename StorageType, typename DomainType, typename SegmentationType>
void init_quantities(StorageType & storage, DomainType const & my_domain, SegmentationType const & my_segmentation, double n_plus, double p_plus)
{
  typedef typename viennagrid::result_of::cell<SegmentationType>::type CellType;
  
  //
  // Init permittivity
  //
  viennafvm::set_quantity_region( viennadata::accessor<permittivity_key, bool, CellType>(storage, permittivity_key()), my_domain, true);               // permittivity is (for simplicity) defined everywhere
  viennafvm::set_quantity_value( viennadata::accessor<permittivity_key, double, CellType>(storage, permittivity_key()), my_domain, 11.7 * 8.854e-12);                // relative permittivity of silicon

  //
  // Initialize doping
  //

  // donator doping
  viennafvm::set_quantity_region(viennadata::accessor<donator_doping_key, bool, CellType>(storage, donator_doping_key()), my_domain, true);

  viennafvm::set_quantity_value(viennadata::accessor<donator_doping_key, double, CellType>(storage, donator_doping_key()), my_domain,  n_plus);
  viennafvm::set_quantity_value(viennadata::accessor<donator_doping_key, double, CellType>(storage, donator_doping_key()), my_segmentation(3),  1e32/p_plus);

  // acceptor doping
  viennafvm::set_quantity_region(viennadata::accessor<acceptor_doping_key, bool, CellType>(storage, acceptor_doping_key()), my_domain, true);

  viennafvm::set_quantity_value(viennadata::accessor<acceptor_doping_key, double, CellType>(storage, acceptor_doping_key()), my_domain,  1e32/n_plus);
  viennafvm::set_quantity_value(viennadata::accessor<acceptor_doping_key, double, CellType>(storage, acceptor_doping_key()), my_segmentation(3),  p_plus);

  // built-in potential:
  viennafvm::set_quantity_region(viennadata::accessor<builtin_potential_key, bool, CellType>(storage, builtin_potential_key()), my_domain, true);   // defined everywhere

  viennafvm::set_quantity_value(viennadata::accessor<builtin_potential_key, double, CellType>(storage, builtin_potential_key()), my_segmentation(1), built_in_potential(300, n_plus, 1e32/n_plus));
  viennafvm::set_quantity_value(viennadata::accessor<builtin_potential_key, double, CellType>(storage, builtin_potential_key()), my_segmentation(2), built_in_potential(300, n_plus, 1e32/n_plus));
  viennafvm::set_quantity_value(viennadata::accessor<builtin_potential_key, double, CellType>(storage, builtin_potential_key()), my_segmentation(3), built_in_potential(300, 1e32/p_plus, p_plus));
  viennafvm::set_quantity_value(viennadata::accessor<builtin_potential_key, double, CellType>(storage, builtin_potential_key()), my_segmentation(4), built_in_potential(300, n_plus, 1e32/n_plus));
  viennafvm::set_quantity_value(viennadata::accessor<builtin_potential_key, double, CellType>(storage, builtin_potential_key()), my_segmentation(5), built_in_potential(300, n_plus, 1e32/n_plus));
}

/** @brief Scales the entire simulation domain (device) by the provided factor. This is accomplished by multiplying all point coordinates with this factor. */
template <typename DomainType>
void scale_domain(DomainType & domain, double factor)
{
  typedef typename viennagrid::result_of::vertex_range<DomainType> ::type VertexContainer;
  typedef typename viennagrid::result_of::iterator<VertexContainer>::type VertexIterator;

  typename viennagrid::result_of::default_point_accessor<DomainType>::type point_accessor = viennagrid::default_point_accessor(domain);

  VertexContainer vertices = viennagrid::elements(domain);
  for ( VertexIterator vit = vertices.begin();
        vit != vertices.end();
        ++vit )
  {
    point_accessor(*vit) *= factor; // scale
  }
}

int main()
{
  typedef double   numeric_type;

  typedef viennagrid::line_1d_domain   DomainType;
  typedef viennagrid::result_of::segmentation<DomainType>::type SegmentationType;

  typedef viennagrid::result_of::cell_tag<DomainType>::type CellTag;

  typedef viennagrid::result_of::element<DomainType, CellTag>::type        CellType;

  typedef viennamath::function_symbol   FunctionSymbol;
  typedef viennamath::equation          Equation;

  typedef viennadata::storage<> StorageType;

  //
  // Create a domain from file
  //
  DomainType my_domain;
  SegmentationType my_segmentation(my_domain);
  StorageType my_storage;

  try
  {
    viennagrid::io::netgen_reader my_reader;
    my_reader(my_domain, my_segmentation, "../examples/data/line23.mesh");
  }
  catch (...)
  {
    std::cerr << "File-Reader failed. Aborting program..." << std::endl;
    return EXIT_FAILURE;
  }

  scale_domain(my_domain, 1e-7);

  //
  // Set initial values
  //
  double n_plus = 1e24;
  double p_plus = 1e10;

  init_quantities(my_storage, my_domain, my_segmentation, n_plus, p_plus);

  //
  // Setting boundary information on domain (see mosfet.in2d for segment indices)
  //
  FunctionSymbol psi(0);   // potential, using id=0
  FunctionSymbol n(1);     // electron concentration, using id=1
  FunctionSymbol p(2);     // hole concentration, using id=2

  // potential:
  double built_in_pot = built_in_potential(300, n_plus, 1e32/n_plus); // should match specification in init_quantities()!
  viennafvm::set_dirichlet_boundary(my_storage, my_segmentation(1), 0.0 + built_in_pot, psi); // Gate contact
  viennafvm::set_dirichlet_boundary(my_storage, my_segmentation(5), 0.0 + built_in_pot, psi); // Source contact

  // electron density
  viennafvm::set_dirichlet_boundary(my_storage, my_segmentation(1), n_plus, n); // Source contact
  viennafvm::set_dirichlet_boundary(my_storage, my_segmentation(5), n_plus, n); // Drain contact

  // hole density
  viennafvm::set_dirichlet_boundary(my_storage, my_segmentation(1), 1e32/n_plus, p); // Source contact
  viennafvm::set_dirichlet_boundary(my_storage, my_segmentation(5), 1e32/n_plus, p); // Drain contact

  //
  // Initial conditions (required for nonlinear problems)
  //
  viennadata::result_of::accessor<StorageType, viennafvm::disable_quantity_key, double, CellType>::type psi_disabled_accessor =
      viennadata::accessor<viennafvm::disable_quantity_key, double, CellType>(my_storage, viennafvm::disable_quantity_key(psi.id()));

  viennadata::result_of::accessor<StorageType, builtin_potential_key, double, CellType>::type psi_source_accessor =
      viennadata::accessor<builtin_potential_key, double, CellType>(my_storage, builtin_potential_key());

  viennadata::result_of::accessor<StorageType, viennafvm::current_iterate_key, double, CellType>::type psi_iterate_accessor =
      viennadata::accessor<viennafvm::current_iterate_key, double, CellType>(my_storage, viennafvm::current_iterate_key(psi.id()));
  
  viennafvm::set_initial_guess(my_domain, psi, psi_disabled_accessor, psi_source_accessor, psi_iterate_accessor);
  //viennafvm::smooth_initial_guess(my_domain, psi, viennafvm::arithmetic_mean_smoother());
  //viennafvm::smooth_initial_guess(my_domain, psi, viennafvm::arithmetic_mean_smoother());


  viennadata::result_of::accessor<StorageType, viennafvm::disable_quantity_key, double, CellType>::type n_disabled_accessor =
      viennadata::accessor<viennafvm::disable_quantity_key, double, CellType>(my_storage, viennafvm::disable_quantity_key(n.id()));

  viennadata::result_of::accessor<StorageType, donator_doping_key, double, CellType>::type n_source_accessor =
      viennadata::accessor<donator_doping_key, double, CellType>(my_storage, donator_doping_key());

  viennadata::result_of::accessor<StorageType, viennafvm::current_iterate_key, double, CellType>::type n_iterate_accessor =
      viennadata::accessor<viennafvm::current_iterate_key, double, CellType>(my_storage, viennafvm::current_iterate_key(n.id()));
      
  viennafvm::set_initial_guess(my_domain, n, n_disabled_accessor, n_source_accessor, n_iterate_accessor);
  //viennafvm::smooth_initial_guess(my_domain, n, viennafvm::geometric_mean_smoother());
  //viennafvm::smooth_initial_guess(my_domain, n, viennafvm::geometric_mean_smoother());

  viennadata::result_of::accessor<StorageType, viennafvm::disable_quantity_key, double, CellType>::type p_disabled_accessor =
      viennadata::accessor<viennafvm::disable_quantity_key, double, CellType>(my_storage, viennafvm::disable_quantity_key(p.id()));

  viennadata::result_of::accessor<StorageType, acceptor_doping_key, double, CellType>::type p_source_accessor =
      viennadata::accessor<acceptor_doping_key, double, CellType>(my_storage, acceptor_doping_key());

  viennadata::result_of::accessor<StorageType, viennafvm::current_iterate_key, double, CellType>::type p_iterate_accessor =
      viennadata::accessor<viennafvm::current_iterate_key, double, CellType>(my_storage, viennafvm::current_iterate_key(p.id()));
  
  viennafvm::set_initial_guess(my_domain, p, p_disabled_accessor, p_source_accessor, p_iterate_accessor);
  //viennafvm::smooth_initial_guess(my_domain, p, viennafvm::geometric_mean_smoother());
  //viennafvm::smooth_initial_guess(my_domain, p, viennafvm::geometric_mean_smoother());



  //
  // Specify PDEs:
  //

  viennafvm::ncell_quantity<CellType, viennamath::expr::interface_type>  permittivity; permittivity.wrap_constant( my_storage, permittivity_key() );
  viennafvm::ncell_quantity<CellType, viennamath::expr::interface_type>  donator_doping; donator_doping.wrap_constant( my_storage, donator_doping_key() );
  viennafvm::ncell_quantity<CellType, viennamath::expr::interface_type>  acceptor_doping; acceptor_doping.wrap_constant( my_storage, acceptor_doping_key() );

  double q  = 1.6e-19;
  double kB = 1.38e-23; // Boltzmann constant
  double mu = 1;        // mobility (constant is fine for the moment)
  double T  = 300;
  double VT = kB * T / q;
  double D  = mu * VT;  //diffusion constant

  // here is all the fun: specify DD system
  Equation laplace_eq = viennamath::make_equation( viennamath::div(permittivity * viennamath::grad(psi)),                     /* = */ 0);
  Equation poisson_eq = viennamath::make_equation( viennamath::div(permittivity * viennamath::grad(psi)),                     /* = */ q * ((n - donator_doping) - (p - acceptor_doping)));
  Equation cont_eq_n  = viennamath::make_equation( viennamath::div(D * viennamath::grad(n) - mu * viennamath::grad(psi) * n), /* = */ 0);
  Equation cont_eq_p  = viennamath::make_equation( viennamath::div(D * viennamath::grad(p) + mu * viennamath::grad(psi) * p), /* = */ 0);

  // Specify the PDE system:
  //viennafvm::linear_pde_system<> pde_system_init;
  //pde_system_init.add_pde(poisson_eq, psi); // equation and associated quantity
  //viennafvm::pde_solver<> pde_solver_init;
  //pde_solver_init(pde_system_init, my_domain);   // weird math happening in here ;-)

  // actual DD system
  viennafvm::linear_pde_system<> pde_system;
  pde_system.add_pde(poisson_eq, psi); // equation and associated quantity
  pde_system.add_pde(cont_eq_n, n);    // equation and associated quantity
  pde_system.add_pde(cont_eq_p, p);    // equation and associated quantity

  pde_system.option(0).damping_term( (n + p) * (-q / VT) );
  pde_system.option(1).geometric_update(true);
  pde_system.option(2).geometric_update(true);

  pde_system.is_linear(false); // temporary solution up until automatic nonlinearity detection is running


  //
  // Create PDE solver instance and run the solver:
  //
  viennafvm::pde_solver<> pde_solver;

  pde_solver(my_storage, pde_system, my_domain);   // weird math happening in here ;-)


  //
  // Writing all solution variables back to domain
  //
  std::vector<long> result_ids(pde_system.size());
  for (std::size_t i=0; i<pde_system.size(); ++i)
    result_ids[i] = pde_system.unknown(i)[0].id();

  viennafvm::io::write_solution_to_VTK_file(my_storage, pde_solver.result(), "nin", my_domain, my_segmentation, result_ids);

  std::cout << "********************************************" << std::endl;
  std::cout << "* MOSFET simulation finished successfully! *" << std::endl;
  std::cout << "********************************************" << std::endl;
  return EXIT_SUCCESS;
}

