#ifndef VIENNAFVM_LINEAR_SOLVE_HPP
#define VIENNAFVM_LINEAR_SOLVE_HPP

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

//ViennaCL includes:
#ifndef VIENNACL_HAVE_UBLAS
 #define VIENNACL_HAVE_UBLAS
#endif

#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/norm_2.hpp"
#include "viennacl/linalg/prod.hpp"

namespace viennafvm
{
  //
  // Solve system of linear equations:
  //
  template <typename MatrixType, typename VectorType>
  VectorType solve(MatrixType const & system_matrix,
                   VectorType const & load_vector)
  {
    typedef typename VectorType::value_type        numeric_type;
    VectorType result(load_vector.size());

    std::cout << "* solve(): Solving linear system" << std::endl;

    result = viennacl::linalg::solve(system_matrix, load_vector, viennacl::linalg::bicgstab_tag());
    std::cout << "* solve(): Residual: " << norm_2(prod(system_matrix, result) - load_vector) / norm_2(load_vector) << std::endl;

    return result;
  }
}

#endif // LINEAR_SOLVE_HPP
