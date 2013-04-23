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
#include "viennafvm/linear_solve.hpp"


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


int main()
{
  typedef double   numeric_type;

  typedef viennagrid::config::triangular_2d                           ConfigType;
  typedef viennagrid::result_of::domain<ConfigType>::type             DomainType;
  typedef typename ConfigType::cell_tag                     CellTag;

  typedef viennagrid::result_of::ncell<ConfigType, CellTag::dim>::type        CellType;
  typedef viennagrid::result_of::ncell_range<DomainType, CellTag::dim>::type  CellContainer;
  typedef viennagrid::result_of::iterator<CellContainer>::type                CellIterator;
  typedef viennagrid::result_of::ncell_range<CellType, 0>::type               VertexOnCellContainer;
  typedef viennagrid::result_of::iterator<VertexOnCellContainer>::type        VertexOnCellIterator;

  typedef boost::numeric::ublas::compressed_matrix<numeric_type>  MatrixType;
  typedef boost::numeric::ublas::vector<numeric_type>             VectorType;

  typedef viennamath::function_symbol   FunctionSymbol;
  typedef viennamath::equation          Equation;

  //
  // Create a domain from file
  //
  DomainType my_domain;

  try
  {
    viennagrid::io::netgen_reader my_reader;
    my_reader(my_domain, "../examples/data/square128.mesh");
  }
  catch (...)
  {
    std::cerr << "File-Reader failed. Aborting program..." << std::endl;
    return EXIT_FAILURE;
  }

  //
  // Specify two PDEs:
  //
  FunctionSymbol u(0, viennamath::unknown_tag<>());   //an unknown function used for PDE specification
  FunctionSymbol v(1, viennamath::unknown_tag<>());   //an unknown function used for PDE specification
  Equation poisson_equ_1 = viennamath::make_equation( viennamath::laplace(u), -1);
  Equation poisson_equ_2 = viennamath::make_equation( viennamath::laplace(v), 0);

  MatrixType system_matrix_1, system_matrix_2;
  VectorType load_vector_1, load_vector_2;

  //
  // Setting boundary information on domain (this should come from device specification)
  //
  //setting some boundary flags:
  CellContainer cells = viennagrid::ncells(my_domain);
  for (CellIterator cit  = cells.begin();
                    cit != cells.end();
                  ++cit)
  {
    VertexOnCellContainer vertices = viennagrid::ncells<0>(*cit);
    for (VertexOnCellIterator vit  = vertices.begin();
                              vit != vertices.end();
                            ++vit)
    {
      //boundary for first equation: Homogeneous Dirichlet everywhere
      if ( (*vit)[0] == 0.0 || (*vit)[0] == 1.0
           || (*vit)[1] == 0.0 || (*vit)[1] == 1.0 )
        viennafvm::set_dirichlet_boundary(*cit, 0.0, 0);  //simulation with ID 0 uses homogeneous boundary data

      //boundary for second equation (ID 1): 0 at left boundary, 1 at right boundary
      if ( (*vit)[0] == 0.0)
        viennafvm::set_dirichlet_boundary(*cit, 0.0, 1);
      else if ( (*vit)[0] == 1.0)
        viennafvm::set_dirichlet_boundary(*cit, 1.0, 1);
    }
  }


  //
  // Create PDE solver functors: (discussion about proper interface required)
  //
  viennafvm::linear_assembler fvm_assembler;


  //
  // Solve system and write solution vector to pde_result:
  // (discussion about proper interface required. Introduce a pde_result class?)
  //
  fvm_assembler(viennafvm::make_linear_pde_system(poisson_equ_1, u),
                my_domain,
                system_matrix_1,
                load_vector_1
              );

  VectorType pde_result_1 = viennafvm::solve(system_matrix_1, load_vector_1);
  viennafvm::io::write_solution_to_VTK_file(pde_result_1, "poisson_2d_1", my_domain, 0);

  fvm_assembler(viennafvm::make_linear_pde_system(poisson_equ_2, v, viennafvm::make_linear_pde_options(1, 1)),
                my_domain,
                system_matrix_2,
                load_vector_2
              );

  //std::cout << system_matrix_2 << std::endl;
  //std::cout << load_vector_2 << std::endl;

  VectorType pde_result_2 = viennafvm::solve(system_matrix_2, load_vector_2);


  //
  // Writing solution back to domain (discussion about proper way of returning a solution required...)
  //
  viennafvm::io::write_solution_to_VTK_file(pde_result_2, "poisson_2d_2", my_domain, 1);

  std::cout << "*****************************************" << std::endl;
  std::cout << "* Poisson solver finished successfully! *" << std::endl;
  std::cout << "*****************************************" << std::endl;
  return EXIT_SUCCESS;
}

