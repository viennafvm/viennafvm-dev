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

  template <typename PDESystemType, typename DomainType, typename VectorType>
  void apply_update(PDESystemType const & pde_system, DomainType const & domain, VectorType const & update)
  {
    typedef typename DomainType::config_type              Config;
    typedef typename Config::cell_tag                     CellTag;

    typedef typename viennagrid::result_of::point<Config>::type                  PointType;
    typedef typename viennagrid::result_of::ncell<Config, CellTag::dim>::type    CellType;

    typedef typename viennagrid::result_of::const_ncell_range<DomainType, CellTag::dim>::type   CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

    typedef typename PDESystemType::mapping_key_type   MappingKeyType;
    typedef typename PDESystemType::boundary_key_type  BoundaryKeyType;

    numeric_type alpha = 0.5;

    for (std::size_t pde_index = 0; pde_index < pde_system.size(); ++pde_index)
    {
      long unknown_id = pde_system.option(pde_index).data_id();

      BoundaryKeyType bnd_key(unknown_id);
      MappingKeyType  map_key(unknown_id);

      viennamath::function_symbol const & u = pde_system.unknown(pde_index)[0];
      numeric_type l2_update_norm = 0;

      CellContainer cells = viennagrid::ncells(domain);
      for (CellIterator cit = cells.begin(); cit != cells.end(); ++cit)
      {
        if (viennadata::access<BoundaryKeyType, bool>(bnd_key)(*cit))  // boundary cell
        {
          //nothing
        }
        else if (viennafvm::is_quantity_enabled(*cit, unknown_id))
        {
          double current_value = get_current_iterate(*cit, u);
          double update_value = update(viennadata::access<MappingKeyType, long>(map_key)(*cit));

          numeric_type new_value = pde_system.option(pde_index).geometric_update()
                                    ? std::pow(current_value, 1.0 - alpha) * std::pow( std::max(current_value + update_value, 0.1), alpha)
                                    : current_value + alpha * update_value;

          l2_update_norm += new_value * new_value;
          set_current_iterate(*cit, u, new_value);
        }
      }

      std::cout << "Update norm for quantity " << pde_index << ": " << std::sqrt(l2_update_norm) << std::endl;
    }
  }

  template <typename PDESystemType, typename DomainType, typename VectorType>
  void transfer_to_solution_vector(PDESystemType const & pde_system, DomainType const & domain, VectorType & result)
  {
    typedef typename DomainType::config_type              Config;
    typedef typename Config::cell_tag                     CellTag;

    typedef typename viennagrid::result_of::point<Config>::type                  PointType;
    typedef typename viennagrid::result_of::ncell<Config, CellTag::dim>::type    CellType;

    typedef typename viennagrid::result_of::const_ncell_range<DomainType, CellTag::dim>::type   CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

    typedef typename PDESystemType::mapping_key_type   MappingKeyType;
    typedef typename PDESystemType::boundary_key_type  BoundaryKeyType;

    for (std::size_t pde_index = 0; pde_index < pde_system.size(); ++pde_index)
    {
      long unknown_id = pde_system.option(pde_index).data_id();

      BoundaryKeyType bnd_key(unknown_id);
      MappingKeyType  map_key(unknown_id);

      CellContainer cells = viennagrid::ncells(domain);
      for (CellIterator cit = cells.begin(); cit != cells.end(); ++cit)
      {
        if (viennadata::access<BoundaryKeyType, bool>(bnd_key)(*cit))  // boundary cell
        {
          //nothing
        }
        else if (viennafvm::is_quantity_enabled(*cit, unknown_id))
        {
          result(viennadata::access<MappingKeyType, long>(map_key)(*cit)) = get_current_iterate(*cit, pde_system.unknown(pde_index)[0]);
        }
      }
    }
  }


  template <typename MatrixType = boost::numeric::ublas::compressed_matrix<viennafvm::numeric_type>,
            typename VectorType = boost::numeric::ublas::vector<viennafvm::numeric_type> >
  class pde_solver
  {
    public:

      template <typename PDESystemType, typename DomainType>
      void operator()(PDESystemType const & pde_system, DomainType const & domain)
      {
        bool is_linear = pde_system.is_linear(); //TODO: Replace with an automatic detection

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
          VectorType update = viennafvm::solve(system_matrix, load_vector);

          apply_update(pde_system, domain, update);
          transfer_to_solution_vector(pde_system, domain, result_);
        }
        else // nonlinear
        {
          std::size_t iter_max = 50;
          for (std::size_t iter=0; iter < iter_max; ++iter)
          {
            MatrixType system_matrix;
            VectorType load_vector;

            // assemble linearized systems
            viennafvm::linear_assembler fvm_assembler;

            fvm_assembler(pde_system, domain, system_matrix, load_vector);

            result_.resize(load_vector.size());
            VectorType update = viennafvm::solve(system_matrix, load_vector);

            apply_update(pde_system, domain, update);
          }
          transfer_to_solution_vector(pde_system, domain, result_);
        }

      }

      VectorType const & result() { return result_; }

    private:
      VectorType result_;
  };

}

#endif // VIENNAFVM_PDE_SOLVER_HPP
