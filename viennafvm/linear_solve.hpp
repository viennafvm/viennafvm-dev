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
#include "viennacl/linalg/ilu.hpp"

namespace viennafvm
{

  template <typename MatrixType, typename VectorType>
  VectorType row_normalize_system(MatrixType & system_matrix, VectorType & rhs)
  {
      typedef typename MatrixType::iterator1    Iterator1;
      typedef typename MatrixType::iterator2    Iterator2;

      VectorType scale_factors(rhs.size());

      //preprocessing: scale rows such that ||a_i|| = 1
      for (Iterator1 iter1  = system_matrix.begin1();
          iter1 != system_matrix.end1();
        ++iter1)
      {
        double row_norm = 0.0;
        for (Iterator2 iter2  = iter1.begin();
          iter2 != iter1.end();
        ++iter2)
        {
          row_norm += *iter2 * *iter2;
        }

        row_norm = sqrt(row_norm);

        if (system_matrix(iter1.index1(), iter1.index1()) < 0.0)
          row_norm *= -1.0;

        for (Iterator2 iter2  = iter1.begin();
          iter2 != iter1.end();
        ++iter2)
        {
          *iter2 /= row_norm;
        }
        rhs(iter1.index1()) /= row_norm;
        scale_factors(iter1.index1()) = 1.0 / row_norm;
      }

      return scale_factors;
  }

  //
  // Solve system of linear equations:
  //
  template <typename MatrixType, typename VectorType>
  VectorType solve(MatrixType const & system_matrix,
                   VectorType const & load_vector)
  {
    typedef typename VectorType::value_type        numeric_type;
    VectorType result(load_vector.size());

    MatrixType system_matrix_2(system_matrix);
    VectorType load_vector_2(load_vector);

    row_normalize_system(system_matrix_2, load_vector_2);

    std::cout << "* solve(): Solving linear system" << std::endl;

    viennacl::linalg::ilu0_tag precond_tag;
    viennacl::linalg::ilu0_precond<MatrixType> preconditioner(system_matrix, precond_tag);

    result = viennacl::linalg::solve(system_matrix_2, load_vector_2, viennacl::linalg::bicgstab_tag(), preconditioner);
    std::cout << "* solve(): Residual: " << norm_2(prod(system_matrix, result) - load_vector) / norm_2(load_vector) << std::endl;

    return result;
  }
}

#endif // LINEAR_SOLVE_HPP
