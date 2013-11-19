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
#include "viennafvm/problem_description.hpp"

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


const unsigned int Source = 1;
const unsigned int Channel = 2;
const unsigned int Drain = 3;
const unsigned int Oxide = 4;
const unsigned int Gate = 5;
const unsigned int Body = 6;
const unsigned int BodyContact = 7;
const unsigned int SourceContact = 8;
const unsigned int DrainContact = 9;



template <typename SegmentationType, typename ProblemDescriptionT>
void init_quantities(SegmentationType const & segmentation, ProblemDescriptionT & problem_description)
{
  typedef typename ProblemDescriptionT::quantity_type    QuantityType;

  //
  // Electrostatic potential
  //
  QuantityType & potential = problem_description.add_quantity(names::potential());
  viennafvm::set_initial_value(potential, segmentation(Gate),          built_in_potential(300, 1e24, 1e8)); // gate
  viennafvm::set_initial_value(potential, segmentation(SourceContact), built_in_potential(300, 1e24, 1e8)); // source contact
  viennafvm::set_initial_value(potential, segmentation(Oxide),         built_in_potential(300, 1e24, 1e8)); // oxide (for simplicity set to same as gate)
  viennafvm::set_initial_value(potential, segmentation(DrainContact),  built_in_potential(300, 1e24, 1e8)); // drain contact
  viennafvm::set_initial_value(potential, segmentation(Source),        built_in_potential(300, 1e24, 1e8)); // source
  viennafvm::set_initial_value(potential, segmentation(Drain),         built_in_potential(300, 1e24, 1e8)); // drain
  viennafvm::set_initial_value(potential, segmentation(Channel),       built_in_potential(300, 1e8,  1e8)); // Channel
  viennafvm::set_initial_value(potential, segmentation(Body),          built_in_potential(300, 1e8,  1e8)); // body
  viennafvm::set_initial_value(potential, segmentation(BodyContact),   built_in_potential(300, 1e8,  1e8)); // body contact (floating body)

  viennafvm::set_dirichlet_boundary(potential, segmentation(Gate),          0.2 + built_in_potential(300, 1e24, 1e8)); // Gate contact
  viennafvm::set_dirichlet_boundary(potential, segmentation(SourceContact), 0.0 + built_in_potential(300, 1e24, 1e8)); // Source contact
  viennafvm::set_dirichlet_boundary(potential, segmentation(DrainContact),  0.2 + built_in_potential(300, 1e24, 1e8)); // Drain contact
  viennafvm::set_dirichlet_boundary(potential, segmentation(BodyContact),   0.0 + built_in_potential(300, 1e8,  1e24)); // Body contact

  viennafvm::set_unknown(potential, segmentation(Oxide));
  viennafvm::set_unknown(potential, segmentation(Source));
  viennafvm::set_unknown(potential, segmentation(Drain));
  viennafvm::set_unknown(potential, segmentation(Channel));
  viennafvm::set_unknown(potential, segmentation(Body));


  //
  // Electron Density
  //
  QuantityType & electron_density = problem_description.add_quantity(names::electron_density(), 1e16);   // initialize with intrinsic concentration
  viennafvm::set_initial_value(electron_density, segmentation(Source),        1e24);  // source
  viennafvm::set_initial_value(electron_density, segmentation(SourceContact), 1e24);  // source contact
  viennafvm::set_initial_value(electron_density, segmentation(Drain),         1e24);  // drain
  viennafvm::set_initial_value(electron_density, segmentation(DrainContact),  1e24);  // drain contact
  viennafvm::set_initial_value(electron_density, segmentation(Channel),       1e8);   // channel
  viennafvm::set_initial_value(electron_density, segmentation(Oxide),         1e8);   // oxide (this is for averaging purposes)
  viennafvm::set_initial_value(electron_density, segmentation(Body),          1e8);   // body
  viennafvm::set_initial_value(electron_density, segmentation(BodyContact),   1e8);   // body contact (floating body)

  viennafvm::set_dirichlet_boundary(electron_density, segmentation(SourceContact), 1e24);   // Source contact
  viennafvm::set_dirichlet_boundary(electron_density, segmentation(DrainContact),  1e24);   // Drain contact
  viennafvm::set_dirichlet_boundary(electron_density, segmentation(BodyContact),   1e8);    // Body contact

  viennafvm::set_unknown(electron_density, segmentation(Source));
  viennafvm::set_unknown(electron_density, segmentation(Drain));
  viennafvm::set_unknown(electron_density, segmentation(Channel));
  viennafvm::set_unknown(electron_density, segmentation(Body));

  //
  // Hole Density
  //
  QuantityType & hole_density = problem_description.add_quantity(names::hole_density(), 1e16);   // initialize with intrinsic concentration

  viennafvm::set_initial_value(hole_density, segmentation(Source),        1e8);  // source
  viennafvm::set_initial_value(hole_density, segmentation(SourceContact), 1e8);  // source
  viennafvm::set_initial_value(hole_density, segmentation(Drain),         1e8);  // drain
  viennafvm::set_initial_value(hole_density, segmentation(DrainContact),  1e8);  // drain
  viennafvm::set_initial_value(hole_density, segmentation(Channel),       1e24); // channel
  viennafvm::set_initial_value(hole_density, segmentation(Oxide),         1e24); // oxide (this is for averaging purposes)
  viennafvm::set_initial_value(hole_density, segmentation(Body),          1e24); // body
  viennafvm::set_initial_value(hole_density, segmentation(BodyContact),   1e24); // body contact (floating body)

  viennafvm::set_dirichlet_boundary(hole_density, segmentation(SourceContact), 1e8);   // Source contact
  viennafvm::set_dirichlet_boundary(hole_density, segmentation(DrainContact),  1e8);   // Drain contact
  viennafvm::set_dirichlet_boundary(hole_density, segmentation(BodyContact),   1e24);    // Body contact

  viennafvm::set_unknown(hole_density, segmentation(Source));
  viennafvm::set_unknown(hole_density, segmentation(Drain));
  viennafvm::set_unknown(hole_density, segmentation(Channel));
  viennafvm::set_unknown(hole_density, segmentation(Body));

  //
  // Init permittivity
  //
  QuantityType & permittivity = problem_description.add_quantity(names::permittivity(), 11.7 * 8.854e-12);   // initialize with permittivity of silicon
  viennafvm::set_initial_value(permittivity, segmentation(Oxide), 15.6 * 8.854e-12);   // permittivty of HfO2

  //
  // Initialize doping
  //

  // donator doping
  QuantityType & donator_doping = problem_description.add_quantity(names::donator_doping(), 1e16);   // initialize with intrinsic concentration
  viennafvm::set_initial_value(donator_doping, segmentation(Source),      1e24); // source
  viennafvm::set_initial_value(donator_doping, segmentation(Drain),       1e24); // drain
  viennafvm::set_initial_value(donator_doping, segmentation(Body),        1e18); // body
  viennafvm::set_initial_value(donator_doping, segmentation(BodyContact), 1e18); // body contact (floating body)
  viennafvm::set_initial_value(donator_doping, segmentation(Channel),     1e18); // channel

  // acceptor doping
  QuantityType & acceptor_doping = problem_description.add_quantity(names::acceptor_doping(), 1e16);   // initialize with intrinsic concentration
  viennafvm::set_initial_value(acceptor_doping, segmentation(Source),      1e8); // source
  viennafvm::set_initial_value(acceptor_doping, segmentation(Drain),       1e8); // drain
  viennafvm::set_initial_value(acceptor_doping, segmentation(Body),        1e14);      // body
  viennafvm::set_initial_value(acceptor_doping, segmentation(BodyContact), 1e14);      // body contact (floating body)
  viennafvm::set_initial_value(acceptor_doping, segmentation(Channel),     1e14);      // body contact (floating body)

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
  my_vtk_writer(mesh, segments, "mosfet_3d_doping");
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
  bnd_vtk_writer(mesh, segments, "mosfet_3d_boundary_conditions");

  viennagrid::io::vtk_writer<MeshT> init_vtk_writer;
  init_vtk_writer.add_scalar_data_on_cells( init_pot_acc , "potential" );
  init_vtk_writer.add_scalar_data_on_cells( init_n_acc ,   "electrons" );
  init_vtk_writer.add_scalar_data_on_cells( init_p_acc ,   "holes" );
  init_vtk_writer(mesh, segments, "mosfet_3d_initial_conditions");
}
*/



int main(int argc, char* argv[])
{
  if(argc != 2)
  {
      std::cerr << "Missing parameters - Usage: " << argv[0] << " path/to/trigate.mesh" << std::endl;
      return -1;
  }

  typedef double   numeric_type;

  typedef viennagrid::tetrahedral_3d_mesh                       MeshType;
  typedef viennagrid::result_of::segmentation<MeshType>::type   SegmentationType;

  typedef viennagrid::result_of::cell_tag<MeshType>::type            CellTag;
  typedef viennagrid::result_of::element<MeshType, CellTag>::type    CellType;

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
    my_reader(mesh, segmentation, argv[1]);
  }
  catch (...)
  {
    std::cerr << "File-Reader failed. Aborting program..." << std::endl;
    return EXIT_FAILURE;
  }

  viennagrid::scale(mesh, 1e-9); // scale to nanometer

  viennafvm::problem_description<MeshType> problem_desc(mesh);


  //
  // Set initial values
  //
  init_quantities(segmentation, problem_desc);

  //
  // Setting boundary information on mesh (see mosfet.in2d for segment indices)
  //
  FunctionSymbol psi(0);   // potential, using id=0
  FunctionSymbol n(1);     // electron concentration, using id=1
  FunctionSymbol p(2);     // hole concentration, using id=2
  FunctionSymbol permittivity(3);    // permittivity
  FunctionSymbol donator_doping(4);  // donator doping
  FunctionSymbol acceptor_doping(5); // acceptor doping


  //
  // Smooth initial guesses
  //
  for(int si = 0; si < 5; si++)
  {
    viennafvm::smooth_initial_guess(mesh, problem_desc.quantities().at(psi.id()), viennafvm::arithmetic_mean_smoother());
    viennafvm::smooth_initial_guess(mesh, problem_desc.quantities().at(n.id()),   viennafvm::geometric_mean_smoother());
    viennafvm::smooth_initial_guess(mesh, problem_desc.quantities().at(p.id()),   viennafvm::geometric_mean_smoother());
  }

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
  pde_system.add_pde(poisson_eq, psi); // equation and associated quantity
  pde_system.add_pde(cont_eq_n,  n);   // equation and associated quantity
  pde_system.add_pde(cont_eq_p,  p);   // equation and associated quantity

  pde_system.option(0).damping_term( (n + p) * (-q / VT) );
  pde_system.option(1).geometric_update(true);
  pde_system.option(2).geometric_update(true);

  pde_system.is_linear(false); // temporary solution up until automatic nonlinearity detection is running

  //
  // Setup Linear Solver
  //
  viennafvm::linsolv::viennacl  linear_solver;
  linear_solver.solver()         = viennafvm::linsolv::viennacl::solver_ids::bicgstab;
  linear_solver.preconditioner() = viennafvm::linsolv::viennacl::preconditioner_ids::ilu0;

  //
  // Create PDE solver instance and run the solver:
  //
  viennafvm::pde_solver my_solver;
  my_solver.set_damping(1.0);
  my_solver.set_nonlinear_iterations(100);
  my_solver.set_nonlinear_breaktol(1.0E-2);

  my_solver(problem_desc, pde_system, linear_solver);


  //
  // Writing all solution variables back to mesh
  //
  viennafvm::io::write_solution_to_VTK_file(problem_desc.quantities(), "mosfet_3d", mesh, segmentation);

  std::cout << "********************************************" << std::endl;
  std::cout << "* MOSFET simulation finished successfully! *" << std::endl;
  std::cout << "********************************************" << std::endl;
  return EXIT_SUCCESS;
}

