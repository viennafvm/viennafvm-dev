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
  double apply_update(PDESystemType const & pde_system, std::size_t pde_index, DomainType const & domain, VectorType const & update, numeric_type alpha = 0.3)
  {
    typedef typename DomainType::config_type              Config;
    typedef typename Config::cell_tag                     CellTag;

    typedef typename viennagrid::result_of::point<Config>::type                  PointType;
    typedef typename viennagrid::result_of::ncell<Config, CellTag::dim>::type    CellType;

    typedef typename viennagrid::result_of::const_ncell_range<DomainType, CellTag::dim>::type   CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

    typedef typename PDESystemType::mapping_key_type   MappingKeyType;
    typedef typename PDESystemType::boundary_key_type  BoundaryKeyType;

    long unknown_id = pde_system.unknown(pde_index)[0].id();

    BoundaryKeyType bnd_key(unknown_id);
    MappingKeyType  map_key(unknown_id);

    viennamath::function_symbol const & u = pde_system.unknown(pde_index)[0];
    numeric_type l2_update_norm = 0;

    CellContainer cells = viennagrid::ncells(domain);

    // get damping term
    numeric_type A_n = 0.0;
    if (pde_system.option(pde_index).geometric_update())
    {
      for (CellIterator cit = cells.begin(); cit != cells.end(); ++cit)
      {
        if (viennafvm::is_quantity_enabled(*cit, unknown_id))
        {
          double current_value = get_current_iterate(*cit, u);
          double update_value = viennadata::access<BoundaryKeyType, bool>(bnd_key)(*cit)
                                 ? viennadata::access<BoundaryKeyType, numeric_type>(bnd_key)(*cit) - current_value
                                 : update(viennadata::access<MappingKeyType, long>(map_key)(*cit));

          if (current_value != 0)
            A_n = std::max(A_n, std::abs(update_value / current_value));
        }
      }
    }

    // apply update:
    for (CellIterator cit = cells.begin(); cit != cells.end(); ++cit)
    {
      if (viennafvm::is_quantity_enabled(*cit, unknown_id))
      {
        double current_value = get_current_iterate(*cit, u);
        double update_value = viennadata::access<BoundaryKeyType, bool>(bnd_key)(*cit)
                               ? viennadata::access<BoundaryKeyType, numeric_type>(bnd_key)(*cit) - current_value
                               : update(viennadata::access<MappingKeyType, long>(map_key)(*cit));

        numeric_type new_value = current_value + alpha * update_value;

        if (pde_system.option(pde_index).geometric_update())
        {
          if (update_value < 0)
            new_value = current_value + alpha * ( update_value / (1.0 - A_n * ( update_value / current_value ) ));
          else
            new_value = std::pow(current_value, 1.0 - alpha) * std::pow( current_value + update_value, alpha);
        }

        l2_update_norm += (new_value - current_value) * (new_value - current_value);
        set_current_iterate(*cit, u, new_value);
      }
    }

//    std::cout << "* Update norm for quantity " << pde_index << ": " << std::sqrt(l2_update_norm) << std::endl;
    return std::sqrt(l2_update_norm);
  }


  // [JW] do we really need this? it seems we can do the pde_index loop in the calling instance, as we are doing it
  // with the picard/non-linear approach anyway. that should thus also fly with the linear approach..
  //
//  template <typename PDESystemType, typename DomainType, typename VectorType>
//  void apply_update(PDESystemType const & pde_system, DomainType const & domain, VectorType const & update, numeric_type alpha = 0.5)
//  {
//    typedef typename DomainType::config_type              Config;
//    typedef typename Config::cell_tag                     CellTag;

//    typedef typename viennagrid::result_of::point<Config>::type                  PointType;
//    typedef typename viennagrid::result_of::ncell<Config, CellTag::dim>::type    CellType;

//    typedef typename viennagrid::result_of::const_ncell_range<DomainType, CellTag::dim>::type   CellContainer;
//    typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

//    typedef typename PDESystemType::mapping_key_type   MappingKeyType;
//    typedef typename PDESystemType::boundary_key_type  BoundaryKeyType;

//    for (std::size_t pde_index = 0; pde_index < pde_system.size(); ++pde_index)
//    {
//      apply_update(pde_system, pde_index, domain, update, alpha);
//    }
//  }

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
      long unknown_id = pde_system.unknown(pde_index)[0].id();

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

      typedef viennafvm::numeric_type numeric_type;

      pde_solver()
      {
        nonlinear_iterations = 40;
        linear_iterations = 500;
        linear_breaktol = 1e-14;
        damping = 1.0;
      }

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
          VectorType update = viennafvm::solve(system_matrix, load_vector, linear_iterations, linear_breaktol);
          //std::cout << update << std::endl;

          for (std::size_t pde_index = 0; pde_index < pde_system.size(); ++pde_index)
          {
            numeric_type update_norm = apply_update(pde_system, pde_index, domain, update, damping);
            std::cout << "* Update norm for quantity " << pde_index << ": " << update_norm << std::endl;
          }
          transfer_to_solution_vector(pde_system, domain, result_);
        }
        else // nonlinear
        {
          picard_iteration_ = true;

          bool converged = false;
          std::size_t required_nonlinear_iterations = 0;
          for (std::size_t iter=0; iter < nonlinear_iterations; ++iter)
          {
            required_nonlinear_iterations++;
            std::cout << " --- Iteration " << iter << " --- " << std::endl;

            if (picard_iteration_)
            {
              for (std::size_t pde_index = 0; pde_index < pde_system.size(); ++pde_index)
              {
                MatrixType system_matrix;
                VectorType load_vector;

                // assemble linearized systems
                viennafvm::linear_assembler fvm_assembler;

                fvm_assembler(pde_system, pde_index, domain, system_matrix, load_vector);

                //std::cout << system_matrix << std::endl;
                //std::cout << load_vector << std::endl;

                result_.resize(load_vector.size());
                VectorType update = viennafvm::solve(system_matrix, load_vector, linear_iterations, linear_breaktol);
                //std::cout << update << std::endl;

                numeric_type update_norm = apply_update(pde_system, pde_index, domain, update, damping);
                std::cout << "* Update norm for quantity " << pde_index << ": " << update_norm << std::endl;

                if(pde_index == 0) // check if the potential update has converged ..
                {
                    if(update_norm <= nonlinear_breaktol) converged = true;
                }
              }
            }
            else
            {
              throw "not implemented!";
            }

            if(converged) break; // .. the nonlinear for-loop

          } // nonlinear for-loop

          if(converged)
          {
              std::cout << std::endl;
              std::cout << "--------" << std::endl;
              std::cout << "Success: Simulation converged successfully!" << std::endl;
              std::cout << "  Potential update reached the break-tolerance of " << nonlinear_breaktol
                        << " in " << required_nonlinear_iterations << " iterations" << std::endl;
              std::cout << "--------" << std::endl;
          }
          else
          {
              std::cout << std::endl;
              std::cout << "--------" << std::endl;
              std::cout << "Error: Simulation did not converge" << std::endl;
              std::cout << "  Potential update did not reach the break-tolerance of " << nonlinear_breaktol
                        << " in " << nonlinear_iterations << " iterations" << std::endl;
              std::cout << "--------" << std::endl;
          }


          // need to pack all approximations into a single vector:
          std::size_t map_index = create_mapping(pde_system, domain);
          result_.resize(map_index);
          transfer_to_solution_vector(pde_system, domain, result_);
        }

      }

      VectorType const & result() { return result_; }

      std::size_t get_nonlinear_iterations() { return nonlinear_iterations; }
      void set_nonlinear_iterations(std::size_t max_iters) { nonlinear_iterations = max_iters; }
      
      numeric_type get_nonlinear_breaktol() { return nonlinear_iterations; }
      void set_nonlinear_breaktol(numeric_type value) { nonlinear_breaktol = value; }

      std::size_t get_linear_iterations() { return linear_iterations; }
      void set_linear_iterations(std::size_t max_iters) { linear_iterations = max_iters; }            
      
      numeric_type get_linear_breaktol() { return linear_iterations; }
      void set_linear_breaktol(numeric_type value) { linear_breaktol = value; }

      numeric_type get_damping() { return damping; }
      void set_damping(numeric_type value) { damping = value; }

    private:
      VectorType result_;
      bool picard_iteration_;
      std::size_t nonlinear_iterations;
      std::size_t linear_iterations;
      numeric_type      linear_breaktol;
      numeric_type      nonlinear_breaktol;
      numeric_type      damping;
  };

}

#endif // VIENNAFVM_PDE_SOLVER_HPP
