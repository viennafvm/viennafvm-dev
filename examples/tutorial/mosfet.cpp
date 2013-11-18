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
#define NDEBUG

// include necessary system headers
#include <iostream>

// ViennaFVM includes:
#define VIENNAFVM_VERBOSE
#include "viennafvm/forwards.h"
#include "viennafvm/linear_assembler.hpp"
#include "viennafvm/io/vtk_writer.hpp"
#include "viennafvm/boundary.hpp"
#include "viennafvm/pde_solver.hpp"
#include "viennafvm/initial_guess.hpp"
#include "viennafvm/linear_solvers/viennacl.hpp"

// ViennaGrid includes:
#include "viennagrid/mesh/mesh.hpp"
#include "viennagrid/config/default_configs.hpp"
#include "viennagrid/io/netgen_reader.hpp"
#include "viennagrid/io/vtk_writer.hpp"
#include "viennagrid/algorithm/voronoi.hpp"
#include "viennagrid/algorithm/scale.hpp"

// ViennaData includes:
#include "viennadata/api.hpp"

// ViennaMath includes:
#include "viennamath/expression.hpp"

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>


namespace names
{
  std::string permittivity()     { return "Permittivity"; }
  std::string donator_doping()   { return "N_D"; }
  std::string acceptor_doping()  { return "N_A"; }
  std::string potential()        { return "Potential"; }
  std::string electron_density() { return "Electron Density"; }
  std::string hole_density()     { return "Hole Density"; }
}



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

//
// Initialize all quantities involved
//
template <typename SegmentationType, typename ProblemDescriptionT>
void init_quantities(SegmentationType const & segmentation, ProblemDescriptionT & problem_description, double n_plus, double p_plus)
{
  typedef typename ProblemDescriptionT::quantity_type    QuantityType;

  std::size_t cell_num = viennagrid::cells(segmentation.mesh()).size();

  //
  // Electrostatic potential
  //
  QuantityType potential(names::potential(), cell_num);
  viennafvm::set_initial_value(potential, segmentation(1), built_in_potential(300, n_plus, 1e32/n_plus)); // gate
  viennafvm::set_initial_value(potential, segmentation(2), built_in_potential(300, n_plus, 1e32/n_plus)); // source contact
  viennafvm::set_initial_value(potential, segmentation(3), built_in_potential(300, n_plus, 1e32/n_plus)); // oxide (for simplicity set to same as gate)
  viennafvm::set_initial_value(potential, segmentation(4), built_in_potential(300, n_plus, 1e32/n_plus)); // drain contact
  viennafvm::set_initial_value(potential, segmentation(5), built_in_potential(300, n_plus, 1e32/n_plus)); // source
  viennafvm::set_initial_value(potential, segmentation(6), built_in_potential(300, n_plus, 1e32/n_plus)); // drain
  viennafvm::set_initial_value(potential, segmentation(7), built_in_potential(300, 1e32/p_plus, p_plus)); // body
  viennafvm::set_initial_value(potential, segmentation(8), built_in_potential(300, 1e32/p_plus, p_plus)); // body contact (floating body)

  viennafvm::set_dirichlet_boundary(potential, segmentation(1), 0.2 + built_in_potential(300, n_plus, 1e32/n_plus)); // Gate contact
  viennafvm::set_dirichlet_boundary(potential, segmentation(2), 0.0 + built_in_potential(300, n_plus, 1e32/n_plus)); // Source contact
  viennafvm::set_dirichlet_boundary(potential, segmentation(4), 0.2 + built_in_potential(300, n_plus, 1e32/n_plus)); // Drain contact
  viennafvm::set_dirichlet_boundary(potential, segmentation(8), 0.0 + built_in_potential(300, 1e32/p_plus, p_plus)); // Body contact

  viennafvm::set_unknown(potential, segmentation(3));
  viennafvm::set_unknown(potential, segmentation(5));
  viennafvm::set_unknown(potential, segmentation(6));
  viennafvm::set_unknown(potential, segmentation(7));

  problem_description.add_quantity(potential);


  //
  // Electron Density
  //
  QuantityType electron_density(names::electron_density(), cell_num, 1e16);   // initialize with intrinsic concentration
  viennafvm::set_initial_value(electron_density, segmentation(5), n_plus);      // source
  viennafvm::set_initial_value(electron_density, segmentation(6), n_plus);      // drain
  viennafvm::set_initial_value(electron_density, segmentation(7), 1e32/p_plus); // body
  viennafvm::set_initial_value(electron_density, segmentation(8), 1e32/p_plus); // body contact (floating body)

  viennafvm::set_dirichlet_boundary(electron_density, segmentation(2), n_plus);      // Source contact
  viennafvm::set_dirichlet_boundary(electron_density, segmentation(4), n_plus);      // Drain contact
  viennafvm::set_dirichlet_boundary(electron_density, segmentation(8), 1e32/p_plus); // Body contact

  viennafvm::set_unknown(electron_density, segmentation(5));
  viennafvm::set_unknown(electron_density, segmentation(6));
  viennafvm::set_unknown(electron_density, segmentation(7));

  problem_description.add_quantity(electron_density);


  //
  // Hole Density
  //
  QuantityType hole_density(names::hole_density(), cell_num, 1e16);   // initialize with intrinsic concentration
  viennafvm::set_initial_value(hole_density, segmentation(5), 1e32/n_plus); // source
  viennafvm::set_initial_value(hole_density, segmentation(6), 1e32/n_plus); // drain
  viennafvm::set_initial_value(hole_density, segmentation(7), p_plus); // body
  viennafvm::set_initial_value(hole_density, segmentation(8), p_plus); // body contact (floating body)

  viennafvm::set_dirichlet_boundary(hole_density, segmentation(2), 1e32/n_plus); // Source contact
  viennafvm::set_dirichlet_boundary(hole_density, segmentation(4), 1e32/n_plus); // Drain contact
  viennafvm::set_dirichlet_boundary(hole_density, segmentation(8), p_plus); // Body contact

  viennafvm::set_unknown(hole_density, segmentation(5));
  viennafvm::set_unknown(hole_density, segmentation(6));
  viennafvm::set_unknown(hole_density, segmentation(7));

  problem_description.add_quantity(hole_density);

  //
  // Init permittivity
  //
  QuantityType permittivity(names::permittivity(), cell_num, 11.7 * 8.854e-12);   // initialize with permittivity of silicon
  viennafvm::set_initial_value(permittivity, segmentation(3), 15.6 * 8.854e-12);   // permittivty of HfO2
  problem_description.add_quantity(permittivity);

  //
  // Initialize doping
  //

  // donator doping
  QuantityType donator_doping(names::donator_doping(), cell_num, 1e16);   // initialize with intrinsic concentration
  viennafvm::set_initial_value(donator_doping, segmentation(5), n_plus);      // source
  viennafvm::set_initial_value(donator_doping, segmentation(6), n_plus);      // drain
  viennafvm::set_initial_value(donator_doping, segmentation(7), 1e32/p_plus); // body
  viennafvm::set_initial_value(donator_doping, segmentation(8), 1e32/p_plus); // body contact (floating body)
  problem_description.add_quantity(donator_doping);

  // acceptor doping
  QuantityType acceptor_doping(names::acceptor_doping(), cell_num, 1e16);   // initialize with intrinsic concentration
  viennafvm::set_initial_value(acceptor_doping, segmentation(5), 1e32/n_plus); // source
  viennafvm::set_initial_value(acceptor_doping, segmentation(6), 1e32/n_plus); // drain
  viennafvm::set_initial_value(acceptor_doping, segmentation(7), p_plus);      // body
  viennafvm::set_initial_value(acceptor_doping, segmentation(8), p_plus);      // body contact (floating body)
  problem_description.add_quantity(acceptor_doping);

}

/*
template<typename MeshT, typename SegmentationT, typename StorageT>
void write_device_doping(MeshT& mesh, SegmentationT& segments, StorageT& storage)
{
  typedef typename viennagrid::result_of::cell_tag<MeshT>::type            CellTag;
  typedef typename viennagrid::result_of::element<MeshT, CellTag>::type    CellType;

  typedef typename viennadata::result_of::accessor<StorageT, donator_doping_key, double, CellType>::type  DonatorAccessor;
  typedef typename viennadata::result_of::accessor<StorageT, acceptor_doping_key, double, CellType>::type AcceptorAccessor;

  DonatorAccessor  donator_acc  = viennadata::make_accessor(storage, donator_doping_key());
  AcceptorAccessor acceptor_acc = viennadata::make_accessor(storage, acceptor_doping_key());

  viennagrid::io::vtk_writer<MeshT> my_vtk_writer;
  my_vtk_writer.add_scalar_data_on_cells( donator_acc , "donators" );
  my_vtk_writer.add_scalar_data_on_cells( acceptor_acc , "acceptors" );
  my_vtk_writer(mesh, segments, "mosfet_doping");
}

template<typename MeshT, typename SegmentationT, typename StorageT>
void write_device_initial_guesses(MeshT& mesh, SegmentationT& segments, StorageT& storage)
{
  typedef viennafvm::boundary_key                                                                         BoundaryKey;
  typedef viennafvm::current_iterate_key                                                                  IterateKey;

  typedef typename viennagrid::result_of::cell_tag<MeshT>::type                                         CellTag;
  typedef typename viennagrid::result_of::element<MeshT, CellTag>::type                                 CellType;

  typedef typename viennadata::result_of::accessor<StorageT, BoundaryKey, double, CellType>::type         BoundaryAccessor;
  typedef typename viennadata::result_of::accessor<StorageT, IterateKey, double, CellType>::type          InitGuessAccessor;

  BoundaryAccessor  bnd_pot_acc  = viennadata::make_accessor(storage, BoundaryKey(0));
  InitGuessAccessor init_pot_acc = viennadata::make_accessor(storage, IterateKey(0));

  BoundaryAccessor  bnd_n_acc  = viennadata::make_accessor(storage, BoundaryKey(1));
  InitGuessAccessor init_n_acc = viennadata::make_accessor(storage, IterateKey(1));

  BoundaryAccessor  bnd_p_acc  = viennadata::make_accessor(storage, BoundaryKey(2));
  InitGuessAccessor init_p_acc = viennadata::make_accessor(storage, IterateKey(2));

  viennagrid::io::vtk_writer<MeshT> bnd_vtk_writer;
  bnd_vtk_writer.add_scalar_data_on_cells( bnd_pot_acc , "potential" );
  bnd_vtk_writer.add_scalar_data_on_cells( bnd_n_acc ,   "electrons" );
  bnd_vtk_writer.add_scalar_data_on_cells( bnd_p_acc ,   "holes" );
  bnd_vtk_writer(mesh, segments, "mosfet_boundary_conditions");

  viennagrid::io::vtk_writer<MeshT> init_vtk_writer;
  init_vtk_writer.add_scalar_data_on_cells( init_pot_acc , "potential" );
  init_vtk_writer.add_scalar_data_on_cells( init_n_acc ,   "electrons" );
  init_vtk_writer.add_scalar_data_on_cells( init_p_acc ,   "holes" );
  init_vtk_writer(mesh, segments, "mosfet_initial_conditions");
}
*/


int main()
{
  typedef double   numeric_type;

  typedef viennagrid::triangular_2d_mesh                           MeshType;
  typedef viennagrid::result_of::segmentation<MeshType>::type    SegmentationType;

  typedef viennagrid::result_of::cell_tag<MeshType>::type            CellTag;
  typedef viennagrid::result_of::element<MeshType, CellTag>::type    CellType;

  typedef viennamath::function_symbol   FunctionSymbol;
  typedef viennamath::equation          Equation;

  //
  // Create a domain from file
  //
  MeshType mesh;
  SegmentationType segmentation(mesh);

  try
  {
    viennagrid::io::netgen_reader my_reader;
    my_reader(mesh, segmentation, "../examples/data/mosfet.mesh");
  }
  catch (...)
  {
    std::cerr << "File-Reader failed. Aborting program..." << std::endl;
    return EXIT_FAILURE;
  }

  viennagrid::scale(mesh, 1e-9); // scale to nanometer

  viennafvm::pde_solver<MeshType> pde_solver(mesh);


  //
  // Set initial values
  //
  double n_plus = 1e24;
  double p_plus = 1e20;

  init_quantities(segmentation, pde_solver, n_plus, p_plus);

  //
  // Setting boundary information on mesh (see mosfet.in2d for segment indices)
  //
  FunctionSymbol psi(0);   // potential, using id=0
  FunctionSymbol n(1);     // electron concentration, using id=1
  FunctionSymbol p(2);     // hole concentration, using id=2
  FunctionSymbol permittivity(3);        // permittivity
  FunctionSymbol donator_doping(4);      // donator doping
  FunctionSymbol acceptor_doping(5); // acceptor doping

  //
  // Write Doping and initial guesses/boundary conditions to VTK output files for inspection
  //
  //write_device_doping(mesh, segmentation, pde_solver.storage());
  //write_device_initial_guesses(mesh, segmentation, pde_solver.storage());

  //
  // Specify PDEs:
  //

  double q  = 1.6e-19;
  double kB = 1.38e-23; // Boltzmann constant
  double mu = 1;        // mobility (constant is fine for the moment)
  double T  = 300;
  double VT = kB * T / q;
  double D  = mu * VT;  //diffusion constant

  // here is all the fun: specify DD system
  Equation poisson_eq = viennamath::make_equation( viennamath::div(permittivity * viennamath::grad(psi)),                     /* = */ q * ((n - donator_doping) - (p - acceptor_doping)));
  Equation cont_eq_n  = viennamath::make_equation( viennamath::div(D * viennamath::grad(n) - mu * viennamath::grad(psi) * n), /* = */ 0);
  Equation cont_eq_p  = viennamath::make_equation( viennamath::div(D * viennamath::grad(p) + mu * viennamath::grad(psi) * p), /* = */ 0);

  // Specify the PDE system:
  viennafvm::linear_pde_system<> pde_system;
  pde_system.add_pde(poisson_eq, psi); pde_system.option(0).damping_term( (n + p) * (-q / VT) );
  pde_system.add_pde(cont_eq_n,  n);   pde_system.option(1).geometric_update(true);
  pde_system.add_pde(cont_eq_p,  p);   pde_system.option(2).geometric_update(true);

  pde_system.is_linear(false); // temporary solution up until automatic nonlinearity detection is running

  //
  // Setup Linear Solver
  //
  viennafvm::linsolv::viennacl  linear_solver;

  //
  // Create PDE solver instance and run the solver:
  //

  pde_solver(pde_system, linear_solver);   // weird math happening in here ;-)


  //
  // Writing all solution variables back to mesh
  //
  viennafvm::io::write_solution_to_VTK_file(pde_solver.quantities(), "mosfet", mesh, segmentation);

  std::cout << "********************************************" << std::endl;
  std::cout << "* MOSFET simulation finished successfully! *" << std::endl;
  std::cout << "********************************************" << std::endl;
  return EXIT_SUCCESS;
}

