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
#define VIENNAFVM_VERBOSE
#include "viennafvm/forwards.h"
#include "viennafvm/linear_assembler.hpp"
#include "viennafvm/io/vtk_writer.hpp"
#include "viennafvm/boundary.hpp"
#include "viennafvm/pde_solver.hpp"
#include "viennafvm/initial_guess.hpp"
#include "viennafvm/linear_solvers/viennacl.hpp"

// ViennaGrid includes:
#include "viennagrid/config/default_configs.hpp"
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

template <typename SegmentationType, typename StorageType>
void init_quantities( SegmentationType const & segmentation, StorageType & storage )
{
  typedef typename viennagrid::result_of::cell<SegmentationType>::type CellType;

  //
  // Init permittivity
  //
  double eps0 = 8.854e-12;
  double epsr_silicon = 11.7; // silicon!
  double eps_silicon = eps0 * epsr_silicon;


  // donator doping
  viennafvm::set_quantity_region( segmentation(1), storage, permittivity_key(), false);

  viennafvm::set_quantity_region( segmentation(2), storage, permittivity_key(), true);
  viennafvm::set_quantity_value( segmentation(2),  storage, permittivity_key(), eps_silicon);

  viennafvm::set_quantity_region( segmentation(3), storage, permittivity_key(), true);
  viennafvm::set_quantity_value( segmentation(3),  storage, permittivity_key(), eps_silicon);

  viennafvm::set_quantity_region( segmentation(4), storage, permittivity_key(), true);
  viennafvm::set_quantity_value( segmentation(4),  storage, permittivity_key(), eps_silicon);

  viennafvm::set_quantity_region( segmentation(5), storage, permittivity_key(), false);
}

/** @brief Scales the entire simulation mesh (device) by the provided factor. This is accomplished by multiplying all point coordinates with this factor. */
template <typename MeshType>
void scale_mesh(MeshType & mesh, double factor)
{
  typedef typename viennagrid::result_of::vertex_range<MeshType> ::type VertexContainer;
  typedef typename viennagrid::result_of::iterator<VertexContainer>::type VertexIterator;

  typename viennagrid::result_of::default_point_accessor<MeshType>::type point_accessor = viennagrid::default_point_accessor(mesh);

  VertexContainer vertices(mesh);
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

  typedef viennagrid::triangular_2d_mesh   MeshType;
  typedef viennagrid::result_of::segmentation<MeshType>::type SegmentationType;

  typedef viennagrid::result_of::cell_tag<MeshType>::type CellTag;

  typedef viennagrid::result_of::element<MeshType, CellTag>::type        CellType;

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

//   for (int i = 0; i < 10; ++i)
//     std::cout << "Segment " << i << " - " << segmentation.segment_present(i) << std::endl;

  scale_mesh(mesh, 1e-9);

  //
  // Create PDE solver instance:
  //
  viennafvm::pde_solver<MeshType> pde_solver(mesh);

  //
  // Assign doping and set initial values
  //
  init_quantities(segmentation, pde_solver.storage());

  //
  // Setting boundary information on mesh (see mosfet.in2d for segment indices)
  //
  FunctionSymbol psi(0);   // potential, using id=0

  // potential:
  viennafvm::set_dirichlet_boundary( segmentation(1), pde_solver.storage(), psi, 0.0);
  viennafvm::set_dirichlet_boundary( segmentation(5), pde_solver.storage(), psi, 0.8);

  //
  // Specify PDEs:
  //

  viennafvm::ncell_quantity<CellType, viennamath::expr::interface_type>  permittivity; permittivity.wrap_constant( pde_solver.storage(), permittivity_key() );


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

  pde_solver(pde_system, linear_solver);   // weird math happening in here ;-)


  //
  // Writing all solution variables back to mesh
  //
  std::vector<long> result_ids(pde_system.size());
  for (std::size_t i=0; i<pde_system.size(); ++i)
    result_ids[i] = pde_system.unknown(i)[0].id();

  viennafvm::io::write_solution_to_VTK_file(pde_solver.result(), "nin_2d_laplace", mesh, segmentation, pde_solver.storage(), result_ids);

  std::cout << "*************************************" << std::endl;
  std::cout << "* Simulation finished successfully! *" << std::endl;
  std::cout << "*************************************" << std::endl;
  return EXIT_SUCCESS;
}

