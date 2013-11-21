#ifndef VIENNAFVM_PDE_SOLVER_HPP
#define VIENNAFVM_PDE_SOLVER_HPP

/* =======================================================================
   Copyright (c) 2011, Institute for Microelectronics, TU Wien
   http://www.iue.tuwien.ac.at
                             -----------------
           ViennaFVM - The Vienna Finite Volume Method Library
                             -----------------

   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at
               Josef Weinbub                   weinbub@iue.tuwien.ac.at
               (add your name here)

   license:    see file LICENSE in the ViennaFVM base directory
======================================================================= */

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>

#ifdef VIENNAFVM_VERBOSE
#include "viennafvm/timer.hpp"
#endif
#include "viennafvm/forwards.h"
#include "viennafvm/quantity.hpp"
#include "viennafvm/linear_assembler.hpp"
#include "viennafvm/linear_solvers/viennacl.hpp"

namespace viennafvm
{

  template <typename PDESystemType, typename DomainType, typename QuantityContainer, typename VectorType>
  double apply_update(PDESystemType const & pde_system, std::size_t pde_index,
                      DomainType const & domain,
                      QuantityContainer & quantities,
                      VectorType const & update, numeric_type alpha = 0.3)
  {
    typedef typename viennagrid::result_of::cell_tag<DomainType>::type CellTag;

    typedef typename viennagrid::result_of::const_element_range<DomainType, CellTag>::type   CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

    typedef typename QuantityContainer::value_type QuantityType;

    long quantity_id = pde_system.unknown(pde_index)[0].id();

    numeric_type l2_update_norm = 0;

    QuantityType & quan = quantities.at(quantity_id);

    CellContainer cells(domain);

    // get damping term
    numeric_type A_n = 0.0;
    if (pde_system.option(pde_index).geometric_update())
    {
      for (CellIterator cit = cells.begin(); cit != cells.end(); ++cit)
      {

        double current_value = quan.get_value(*cit);
        double update_value  = quan.get_value(*cit);

        if (quan.get_unknown_index(*cit) >= 0)
          update_value = update(quan.get_unknown_index(*cit));
        else
          update_value = quan.get_boundary_value(*cit) - current_value;


        if (current_value != 0)
          A_n = std::max(A_n, std::abs(update_value / current_value));
      }
    }

    // apply update:
    for (CellIterator cit = cells.begin(); cit != cells.end(); ++cit)
    {
      double current_value = quan.get_value(*cit);
      double update_value  = quan.get_value(*cit);

      if (quan.get_unknown_index(*cit) >= 0)
        update_value = update(quan.get_unknown_index(*cit));
      else
        update_value = quan.get_boundary_value(*cit) - current_value;

      numeric_type new_value = current_value + alpha * update_value;

      if (pde_system.option(pde_index).geometric_update())
      {
        if (update_value < 0)
          new_value = current_value + alpha * ( update_value / (1.0 - A_n * ( update_value / current_value ) ));
        else
          new_value = std::pow(current_value, 1.0 - alpha) * std::pow( current_value + update_value, alpha);
      }
      else
      {
        new_value = current_value + alpha * update_value;
      }

      l2_update_norm += (new_value - current_value) * (new_value - current_value);
      quan.set_value(*cit, new_value);
    }

//    std::cout << "* Update norm for quantity " << pde_index << ": " << std::sqrt(l2_update_norm) << std::endl;
    return std::sqrt(l2_update_norm);
  }


  class pde_solver
  {
      typedef boost::numeric::ublas::compressed_matrix<viennafvm::numeric_type>    MatrixType;
      typedef boost::numeric::ublas::vector<viennafvm::numeric_type>               VectorType;

    public:
      typedef viennafvm::numeric_type numeric_type;

      pde_solver()
      {
        nonlinear_iterations  = 100;
        nonlinear_breaktol    = 1.0e-3;
        damping               = 1.0;
      }

      template<typename ProblemDescriptionT, typename PDESystemT, typename LinearSolverT>
      void operator()(ProblemDescriptionT & problem_description,
                      PDESystemT const & pde_system,
                      LinearSolverT& linear_solver,
                      std::size_t break_pde = 0)
      {
      #ifdef VIENNAFVM_VERBOSE
        std::streamsize cout_precision = std::cout.precision();
      #endif

        bool is_linear = pde_system.is_linear(); //TODO: Replace with an automatic detection

        if (is_linear)
        {
          for (std::size_t pde_index = 0; pde_index < pde_system.size(); ++pde_index)
          {
          #ifdef VIENNAFVM_VERBOSE
            viennafvm::Timer timer;
            timer.start();
            std::cout << " * Quantity " << pde_index << " : " << std::endl;
            std::cout << " ------------------------------------------" << std::endl;
          #endif

            MatrixType system_matrix;
            VectorType load_vector;

          #ifdef VIENNAFVM_VERBOSE
            viennafvm::Timer subtimer;
            subtimer.start();
          #endif
            viennafvm::linear_assembler fvm_assembler;
            fvm_assembler(pde_system, pde_index, problem_description.mesh(), problem_description.quantities(), system_matrix, load_vector);
          #ifdef VIENNAFVM_VERBOSE
            std::cout.precision(3);
            subtimer.get();
            std::cout << "   Assembly time : " << std::fixed << subtimer.get() << " s" << std::endl;
          #endif

            VectorType update;
            linear_solver(system_matrix, load_vector, update);
          #ifdef VIENNAFVM_VERBOSE
            std::cout << "   Precond time  : " << std::fixed << linear_solver.last_pc_time() << " s" << std::endl;
            std::cout << "   Solver time   : " << std::fixed << linear_solver.last_solver_time() << " s" << std::endl;
          #endif

          #ifdef VIENNAFVM_VERBOSE
            subtimer.start();
            numeric_type update_norm = apply_update(pde_system, pde_index, problem_description.mesh(), problem_description.quantities(), update, damping);
          #else
            apply_update(pde_system, pde_index, problem_description.mesh(), problem_description.quantities(), update, 1.0);  //this is linear, so there's no need for any damping
          #endif

          #ifdef VIENNAFVM_VERBOSE
            subtimer.get();
            std::cout << "   Update time   : " << std::fixed << subtimer.get() << " s" << std::endl;
          #endif

          #ifdef VIENNAFVM_VERBOSE
            timer.get();
            std::cout << "   Total time    : " << std::fixed << timer.get() << " s" << std::endl;

            std::cout.precision(cout_precision);
            std::cout.unsetf(std::ios_base::floatfield);

            std::cout << "   Solver iters  : " << linear_solver.last_iterations();
            if(linear_solver.last_iterations() == linear_solver.max_iterations())
              std::cout << " ( not converged ) " << std::endl;
            else std::cout << std::endl;

            std::cout << "   Solver error  : " << linear_solver.last_error() << std::endl;

            std::cout << "   Update norm   : "  << update_norm << std::endl;

            std::cout << std::endl;
          #endif
          }
        }
        else // nonlinear
        {
          picard_iteration_ = true;

        #ifdef VIENNAFVM_VERBOSE
          std::vector<double> previous_update_norms(pde_system.size());
        #endif

          bool converged = false;
          std::size_t required_nonlinear_iterations = 0;
          for (std::size_t iter=0; iter < nonlinear_iterations; ++iter)
          {
            required_nonlinear_iterations++;
          #ifdef VIENNAFVM_VERBOSE
            std::cout << " --- Nonlinear iteration " << iter << " --- " << std::endl;
          #endif
            if (picard_iteration_)
            {
              for (std::size_t pde_index = 0; pde_index < pde_system.size(); ++pde_index)
              {
              #ifdef VIENNAFVM_VERBOSE
                viennafvm::Timer timer;
                timer.start();
                std::cout << " * Quantity " << pde_index << " : " << std::endl;
                std::cout << "   ------------------------------------" << std::endl;
              #endif


                MatrixType system_matrix;
                VectorType load_vector;

              #ifdef VIENNAFVM_VERBOSE
                viennafvm::Timer subtimer;
                subtimer.start();
              #endif
                // assemble linearized systems
                viennafvm::linear_assembler fvm_assembler;
                fvm_assembler(pde_system, pde_index, problem_description.mesh(), problem_description.quantities(), system_matrix, load_vector);
              #ifdef VIENNAFVM_VERBOSE
                std::cout.precision(3);
                subtimer.get();
                std::cout << "   Assembly time : " << std::fixed << subtimer.get() << " s" << std::endl;
              #endif

                VectorType update;
                linear_solver(system_matrix, load_vector, update);
              #ifdef VIENNAFVM_VERBOSE
                std::cout << "   Precond time  : " << std::fixed << linear_solver.last_pc_time() << " s" << std::endl;
                std::cout << "   Solver time   : " << std::fixed << linear_solver.last_solver_time() << " s" << std::endl;
              #endif

              #ifdef VIENNAFVM_VERBOSE
                subtimer.start();
              #endif
                numeric_type update_norm = apply_update(pde_system, pde_index, problem_description.mesh(), problem_description.quantities(), update, damping);
              #ifdef VIENNAFVM_VERBOSE
                subtimer.get();
                std::cout << "   Update time   : " << std::fixed << subtimer.get() << " s" << std::endl;
              #endif

              #ifdef VIENNAFVM_VERBOSE
                timer.get();
                std::cout << "   Total time    : " << std::fixed << timer.get() << " s" << std::endl;

                std::cout.precision(cout_precision);
                std::cout.unsetf(std::ios_base::floatfield);

                std::cout << "   Solver iters  : " << linear_solver.last_iterations();
                if(linear_solver.last_iterations() == linear_solver.max_iterations())
                  std::cout << " ( not converged ) " << std::endl;
                else std::cout << std::endl;

                std::cout << "   Solver error  : " << linear_solver.last_error() << std::endl;

                std::string norm_tendency_indicator;
                if(iter == 0)
                {
                  previous_update_norms[pde_index] = update_norm;
                  norm_tendency_indicator = "";
                }
                else
                {
                  if(update_norm > previous_update_norms[pde_index])
                    norm_tendency_indicator = "<up>";
                  else
                  if(update_norm < previous_update_norms[pde_index])
                    norm_tendency_indicator = "<down>";
                  else
                    norm_tendency_indicator = "<=>";
                  previous_update_norms[pde_index] = update_norm;
                }

                std::cout << "   Update norm   : "  << update_norm << " " << norm_tendency_indicator;
                if(pde_index == break_pde)
                  std::cout << " ( **** )" << std::endl;
                else
                  std::cout << std::endl;

                std::cout << std::endl;
              #endif

                if(pde_index == break_pde) // check if the potential update has converged ..
                {
                    if(update_norm <= nonlinear_breaktol) converged = true;
                }
              }
            }
            else
            {
              throw "not implemented!";
            }
          #ifdef VIENNAFVM_VERBOSE
            std::cout << std::endl;
          #endif
            if(converged) break; // .. the nonlinear for-loop

          } // nonlinear for-loop

        #ifdef VIENNAFVM_VERBOSE
          if(converged)
          {
              std::cout << std::endl;
              std::cout << "--------" << std::endl;
              std::cout << "Success: Simulation converged successfully!" << std::endl;
              std::cout << "  Update norm of observed variable reached the break-tolerance of " << nonlinear_breaktol
                        << " in " << required_nonlinear_iterations << " iterations" << std::endl;
              std::cout << "--------" << std::endl;
          }
          else
          {
              std::cout << std::endl;
              std::cout << "--------" << std::endl;
              std::cout << "Warning: Simulation did not converge!" << std::endl;
              std::cout << "  Update norm of observed variable did not reach the break-tolerance of " << nonlinear_breaktol
                        << " in " << nonlinear_iterations << " iterations" << std::endl;
              std::cout << "--------" << std::endl;
          }
        #endif

        }

      }

      std::size_t get_nonlinear_iterations() { return nonlinear_iterations; }
      void set_nonlinear_iterations(std::size_t max_iters) { nonlinear_iterations = max_iters; }

      numeric_type get_nonlinear_breaktol() { return nonlinear_iterations; }
      void set_nonlinear_breaktol(numeric_type value) { nonlinear_breaktol = value; }

      numeric_type get_damping() { return damping; }
      void set_damping(numeric_type value) { damping = value; }

    private:

      bool picard_iteration_;
      std::size_t     nonlinear_iterations;
      numeric_type    nonlinear_breaktol;
      numeric_type    damping;
  };

}

#endif // VIENNAFVM_PDE_SOLVER_HPP

