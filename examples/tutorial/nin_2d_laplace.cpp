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
#include "viennagrid/config/default_configs.hpp"
#include "viennagrid/io/netgen_reader.hpp"
#include "viennagrid/io/vtk_writer.hpp"
#include "viennagrid/algorithm/geometric_transform.hpp"

// ViennaData includes:
#include "viennadata/api.hpp"

// ViennaMath includes:
#include "viennamath/expression.hpp"


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

template <typename SegmentationType, typename ProblemDescriptionT>
void init_quantities( SegmentationType const & segmentation, ProblemDescriptionT & problem_description )
{
  typedef typename ProblemDescriptionT::quantity_type    QuantityType;

  //
  // Electrostatic potential
  //
  QuantityType & potential = problem_description.add_quantity(names::potential());

  viennafvm::set_dirichlet_boundary(potential, segmentation(1), 0.0); // Left contact
  viennafvm::set_dirichlet_boundary(potential, segmentation(5), 0.2); // Right contact

  viennafvm::set_unknown(potential, segmentation(2));
  viennafvm::set_unknown(potential, segmentation(3));
  viennafvm::set_unknown(potential, segmentation(4));

  //
  // Init permittivity
  //
  problem_description.add_quantity(names::permittivity(), 11.7 * 8.854e-12);   // initialize with permittivity of silicon
}


int main()
{
  typedef viennagrid::triangular_2d_mesh   MeshType;
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
    viennagrid::io::netgen_reader my_reader;
    my_reader(mesh, segmentation, "../examples/data/nin2d.mesh");
  }
  catch (...)
  {
    std::cerr << "File-Reader failed. Aborting program..." << std::endl;
    return EXIT_FAILURE;
  }

  viennagrid::scale(mesh, 1e-9); // scale to nanometer

  //
  // Create PDE solver instance:
  //
  viennafvm::problem_description<MeshType> problem_desc(mesh);

  //
  // Assign doping and set initial values
  //
  init_quantities(segmentation, problem_desc);

  //
  // Setting boundary information on mesh (see mosfet.in2d for segment indices)
  //
  FunctionSymbol psi(0);   // potential, using id=0
  FunctionSymbol permittivity(1);   // potential, using id=0

  //
  // Specify PDEs:
  //

  // here is all the fun: specify DD system
  Equation laplace_eq = viennamath::make_equation( viennamath::div(permittivity * viennamath::grad(psi)),                     /* = */ 0);

  viennafvm::linear_pde_system<> pde_system;
  pde_system.add_pde(laplace_eq, psi); // equation and associated quantity
  pde_system.is_linear(true);


  //
  // Setup Linear Solver
  //
  viennafvm::linsolv::viennacl  linear_solver;

  //
  // Run the solver:
  //
  viennafvm::pde_solver my_solver;
  my_solver(problem_desc, pde_system, linear_solver);   // weird math happening in here ;-)


  //
  // Writing all solution variables back to mesh
  //
  viennafvm::io::write_solution_to_VTK_file(problem_desc.quantities(), "nin_2d_laplace", mesh, segmentation);

  std::cout << "*************************************" << std::endl;
  std::cout << "* Simulation finished successfully! *" << std::endl;
  std::cout << "*************************************" << std::endl;
  return EXIT_SUCCESS;
}

