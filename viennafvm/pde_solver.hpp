#ifndef VIENNAFVM_PDE_SOLVER_HPP
#define VIENNAFVM_PDE_SOLVER_HPP

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

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>

#include "viennafvm/forwards.h"
#include "viennafvm/linear_assembler.hpp"
#include "viennafvm/linear_solve.hpp"


namespace viennafvm
{
  namespace detail
  {
    template <typename PDESystem>
    inline bool is_linear(PDESystem const &) { return true; }
  }


  template <typename MatrixType = boost::numeric::ublas::compressed_matrix<viennafvm::numeric_type>,
            typename VectorType = boost::numeric::ublas::vector<viennafvm::numeric_type> >
  class pde_solver
  {
    public:

      template <typename PDESystemType, typename DomainType>
      void operator()(PDESystemType const & pde_system, DomainType const & domain)
      {
        bool is_linear = detail::is_linear(pde_system);

        if (is_linear)
        {
          MatrixType system_matrix;
          VectorType load_vector;

          // do something
          viennafvm::linear_assembler fvm_assembler;

          fvm_assembler(pde_system, domain, system_matrix, load_vector);


          //std::cout << system_matrix << std::endl;
          //std::cout << load_vector << std::endl;

          result_.resize(load_vector.size());
          result_ = viennafvm::solve(system_matrix, load_vector);

        }
        else // nonlinear
        {
          // more work to do. Newton or Picard iteration required
        }

      }

      VectorType const & result() { return result_; }

    private:
      VectorType result_;
  };

}

#endif // VIENNAFVM_PDE_SOLVER_HPP
