#ifndef VIENNAFVM_FLUX_HPP
#define VIENNAFVM_FLUX_HPP

/* =======================================================================
   Copyright (c) 2011, Institute for Microelectronics, TU Wien
   http://www.iue.tuwien.ac.at
                             -----------------
           ViennaFVM - The Vienna Finite Volume Method Library
                             -----------------

   authors:    Karl Rupp                             rupp@iue.tuwien.ac.at
               Josef Weinbub                      weinbub@iue.tuwien.ac.at

               (add your name here)

   license:    To be discussed, see file LICENSE in the ViennaFVM base directory
======================================================================= */

// *** system includes
//
#include <vector>

// *** local includes
//

// *** vienna includes
//
#include "viennamath/expression.hpp"
#include "viennamath/manipulation/substitute.hpp"
#include "viennafvm/cell_quan.hpp"
#include "viennamath/manipulation/diff.hpp"
#include "viennamath/manipulation/eval.hpp"

namespace viennafvm
{

  namespace detail
  {

    /** @brief Transforms a strong formulation of an equation to a weak form, assuming homogeneous Neumann boundary conditions */
    template <typename InterfaceType>
    struct func_symbol_scanner : public viennamath::rt_traversal_interface<InterfaceType>
    {
      public:
        func_symbol_scanner(viennamath::rt_function_symbol<InterfaceType> const & u) : u_(u), has_func_symbol_(false) {}

        void operator()(InterfaceType const * e) const
        {
            viennamath::callback_if_castable< viennamath::rt_function_symbol<InterfaceType> >::apply(e, *this);
            viennamath::callback_if_castable< viennamath::rt_binary_expr<InterfaceType> >::apply(e, *this);
        }

        void operator()(viennamath::rt_function_symbol<InterfaceType> const & fs) const
        {
          typedef typename InterfaceType::numeric_type   NumericType;

          if (fs.deep_equal(&u_)) //this is a surface integral
            has_func_symbol_ = true;
        }

        void operator()(viennamath::rt_binary_expr<InterfaceType> const & bin) const
        {
          typedef typename InterfaceType::numeric_type   NumericType;
          typedef viennamath::op_binary<viennamath::op_plus<NumericType>, InterfaceType>   PlusOperatorType;
          typedef viennamath::op_binary<viennamath::op_minus<NumericType>, InterfaceType>  MinusOperatorType;
          typedef viennamath::op_binary<viennamath::op_mult<NumericType>, InterfaceType>   ProductOperatorType;
          typedef viennamath::op_binary<viennamath::op_div<NumericType>, InterfaceType>    DivisionOperatorType;

          if (    dynamic_cast<const PlusOperatorType *>(bin.op()) != NULL
              || dynamic_cast<const MinusOperatorType *>(bin.op()) != NULL
              || dynamic_cast<const ProductOperatorType *>(bin.op()) != NULL
              || dynamic_cast<const DivisionOperatorType *>(bin.op()) != NULL
                  )
          {
            func_symbol_scanner<InterfaceType> scanner(u_);

            scanner(bin.lhs());
            has_func_symbol_ |= scanner.found();

            scanner.reset();

            scanner(bin.rhs());
            has_func_symbol_ |= scanner.found();
          }
          else //TODO: Add checks!
          {
            throw "Gradient extraction failed";
          }
        }

        bool found() const { return has_func_symbol_; }
        void reset() const { has_func_symbol_ = false; }

      private:
        viennamath::rt_function_symbol<InterfaceType> const & u_;
        mutable bool has_func_symbol_;
    };


    /** @brief Scans for grad(u) a given function u */
    template <typename InterfaceType>
    struct gradient_scanner : public viennamath::rt_traversal_interface<InterfaceType>
    {
      public:
        gradient_scanner(viennamath::rt_function_symbol<InterfaceType> const & u) : u_(u), has_grad_(false) {}

        void operator()(InterfaceType const * e) const
        {
          viennamath::callback_if_castable< viennamath::rt_unary_expr<InterfaceType> >::apply(e, *this);
          viennamath::callback_if_castable< viennamath::rt_binary_expr<InterfaceType> >::apply(e, *this);
        }

        void operator()(viennamath::rt_unary_expr<InterfaceType> const & unary_expr) const
        {
          typedef typename InterfaceType::numeric_type   NumericType;
          typedef viennamath::op_unary<viennamath::op_gradient<NumericType>, InterfaceType>  GradientOperator;

          const GradientOperator * symb_op = dynamic_cast<const GradientOperator *>(unary_expr.op());
          if (symb_op != NULL) //this is a gradient operator
          {
            if (unary_expr.lhs()->deep_equal(&u_))
              has_grad_ = true;
          }
        }

        void operator()(viennamath::rt_binary_expr<InterfaceType> const & bin) const
        {
          typedef typename InterfaceType::numeric_type   NumericType;
          typedef viennamath::op_binary<viennamath::op_plus<NumericType>, InterfaceType>   PlusOperatorType;
          typedef viennamath::op_binary<viennamath::op_minus<NumericType>, InterfaceType>  MinusOperatorType;
          typedef viennamath::op_binary<viennamath::op_mult<NumericType>, InterfaceType>   ProductOperatorType;
          typedef viennamath::op_binary<viennamath::op_div<NumericType>, InterfaceType>    DivisionOperatorType;

          if (    dynamic_cast<const PlusOperatorType *>(bin.op()) != NULL
              || dynamic_cast<const MinusOperatorType *>(bin.op()) != NULL
              || dynamic_cast<const ProductOperatorType *>(bin.op()) != NULL
              || dynamic_cast<const DivisionOperatorType *>(bin.op()) != NULL
                  )
          {
            gradient_scanner<InterfaceType> * scanner = new gradient_scanner<InterfaceType>(u_);
            viennamath::rt_traversal_wrapper<InterfaceType> wrapped_scanner(scanner);

            bin.lhs()->recursive_traversal(wrapped_scanner);
            has_grad_ |= scanner->found();

            scanner->reset();

            bin.rhs()->recursive_traversal(wrapped_scanner);
            has_grad_ |= scanner->found();

          }
          else
          {
            std::cout << "Gradient scanner encountered unexpected expression: " << bin << std::endl;
            throw "Gradient scanning failed";
          }
        }

        bool found() const { return has_grad_; }
        void reset() const { has_grad_ = false; }

      private:
        viennamath::rt_function_symbol<InterfaceType> const & u_;
        mutable bool has_grad_;
    };


    /** @brief Scans for grad(u) a given function u */
    template <typename InterfaceType>
    struct gradient_argument_extractor : public viennamath::rt_traversal_interface<InterfaceType>
    {
      public:
        gradient_argument_extractor(viennamath::rt_function_symbol<InterfaceType> const & u) : u_(u), arg_() {}

        void operator()(InterfaceType const * e) const
        {
          viennamath::callback_if_castable< viennamath::rt_unary_expr<InterfaceType> >::apply(e, *this);
          viennamath::callback_if_castable< viennamath::rt_binary_expr<InterfaceType> >::apply(e, *this);
        }

        void operator()(viennamath::rt_unary_expr<InterfaceType> const & unary_expr) const
        {
          typedef typename InterfaceType::numeric_type   NumericType;
          typedef viennamath::op_unary<viennamath::op_gradient<NumericType>, InterfaceType>  GradientOperator;

          const GradientOperator * symb_op = dynamic_cast<const GradientOperator *>(unary_expr.op());
          if (symb_op != NULL) //this is a gradient operator
          {
            if (unary_expr.lhs()->deep_equal(&u_))
              arg_ = viennamath::rt_expr<InterfaceType>(unary_expr.lhs()->clone());
          }
        }

        void operator()(viennamath::rt_binary_expr<InterfaceType> const & bin) const
        {
          typedef typename InterfaceType::numeric_type   NumericType;
          typedef viennamath::op_binary<viennamath::op_plus<NumericType>, InterfaceType>   PlusOperatorType;
          typedef viennamath::op_binary<viennamath::op_minus<NumericType>, InterfaceType>  MinusOperatorType;
          typedef viennamath::op_binary<viennamath::op_mult<NumericType>, InterfaceType>   ProductOperatorType;
          typedef viennamath::op_binary<viennamath::op_div<NumericType>, InterfaceType>    DivisionOperatorType;

          if (    dynamic_cast<const PlusOperatorType *>(bin.op()) != NULL
              || dynamic_cast<const MinusOperatorType *>(bin.op()) != NULL
              || dynamic_cast<const ProductOperatorType *>(bin.op()) != NULL
              || dynamic_cast<const DivisionOperatorType *>(bin.op()) != NULL
                  )
          {
            // continue iteration on operands (assuming there is only one gradient of the respective function symbol)
            (*this)(bin.lhs());
            (*this)(bin.rhs());
          }
          else
          {
            std::cout << "Gradient extraction encountered unexpected expression: " << bin << std::endl;
            throw "Gradient extraction failed";
          }
        }

        viennamath::rt_expr<InterfaceType> get() const { return arg_; }

      private:
        viennamath::rt_function_symbol<InterfaceType> const & u_;
        mutable viennamath::rt_expr<InterfaceType> arg_;
    };

    /** @brief Scans for the factor in front of grad(u) */
    template <typename InterfaceType>
    struct gradient_prefactor_extractor : public viennamath::rt_traversal_interface<InterfaceType>
    {
      public:
        gradient_prefactor_extractor(viennamath::rt_function_symbol<InterfaceType> const & u) : u_(u), arg_() {}

        void operator()(InterfaceType const * e) const
        {
          viennamath::callback_if_castable< viennamath::rt_unary_expr<InterfaceType> >::apply(e, *this);
          viennamath::callback_if_castable< viennamath::rt_binary_expr<InterfaceType> >::apply(e, *this);
        }

        void operator()(viennamath::rt_unary_expr<InterfaceType> const & unary_expr) const
        {
          typedef typename InterfaceType::numeric_type   NumericType;
          typedef viennamath::op_unary<viennamath::op_gradient<NumericType>, InterfaceType>  GradientOperator;

          const GradientOperator * symb_op = dynamic_cast<const GradientOperator *>(unary_expr.op());
          if (symb_op != NULL) //this is a gradient operator
          {
            if (unary_expr.lhs()->deep_equal(&u_))
              arg_ = viennamath::rt_constant<typename InterfaceType::numeric_type, InterfaceType>(1);
          }
        }

        void operator()(viennamath::rt_binary_expr<InterfaceType> const & bin) const
        {
          typedef typename InterfaceType::numeric_type   NumericType;
          typedef viennamath::op_binary<viennamath::op_plus<NumericType>, InterfaceType>   PlusOperatorType;
          typedef viennamath::op_binary<viennamath::op_minus<NumericType>, InterfaceType>  MinusOperatorType;
          typedef viennamath::op_binary<viennamath::op_mult<NumericType>, InterfaceType>   ProductOperatorType;
          typedef viennamath::op_binary<viennamath::op_div<NumericType>, InterfaceType>    DivisionOperatorType;

          if (   dynamic_cast<const PlusOperatorType *>(bin.op()) != NULL
              || dynamic_cast<const MinusOperatorType *>(bin.op()) != NULL)
          {
            gradient_scanner<InterfaceType> * scanner = new gradient_scanner<InterfaceType>(u_);
            viennamath::rt_traversal_wrapper<InterfaceType> wrapped_scanner(scanner);

            bin.lhs()->recursive_traversal(wrapped_scanner);
            if (scanner->found())
            {
              gradient_argument_extractor<InterfaceType> extractor(u_);
              extractor(bin.lhs());
              arg_ = viennamath::rt_expr<InterfaceType>(extractor.get().get()->clone());
            }
            else
            {
              scanner->reset();

              bin.rhs()->recursive_traversal(wrapped_scanner);
              if (scanner->found())
              {
                gradient_argument_extractor<InterfaceType> extractor(u_);
                extractor(bin.rhs());
                if (dynamic_cast<const MinusOperatorType *>(bin.op()) != NULL)
                {
                  arg_ = viennamath::rt_binary_expr<InterfaceType>(new viennamath::rt_constant<NumericType, InterfaceType>(-1),
                                                                   new viennamath::op_binary<viennamath::op_mult<NumericType>, InterfaceType>(),
                                                                   extractor.get().get()->clone());
                }
                else
                {
                  arg_ = viennamath::rt_expr<InterfaceType>(extractor.get().get()->clone());
                }
              }
            }
          }
          else if (dynamic_cast<const ProductOperatorType *>(bin.op()) != NULL
                || dynamic_cast<const DivisionOperatorType *>(bin.op()) != NULL
                  )
          {
            gradient_scanner<InterfaceType> * scanner = new gradient_scanner<InterfaceType>(u_);
            viennamath::rt_traversal_wrapper<InterfaceType> wrapped_scanner(scanner);

            bin.lhs()->recursive_traversal(wrapped_scanner);
            if (scanner->found())
            {
              gradient_argument_extractor<InterfaceType> extractor(u_);
              extractor(bin.lhs());
              arg_ = viennamath::rt_binary_expr<InterfaceType>(extractor.get().get()->clone(),
                                                               bin.op()->clone(),
                                                               bin.rhs()->clone());
            }
            else
            {
              scanner->reset();

              bin.rhs()->recursive_traversal(wrapped_scanner);
              if (scanner->found())
              {
                gradient_argument_extractor<InterfaceType> extractor(u_);
                extractor(bin.rhs());
                arg_ = viennamath::rt_binary_expr<InterfaceType>(bin.lhs()->clone(),
                                                                 bin.op()->clone(),
                                                                 extractor.get().get()->clone());
              }
            }
          }
          else
          {
            std::cout << "Gradient extraction encountered unexpected expression: " << bin << std::endl;
            throw "Gradient extraction failed";
          }
        }

        viennamath::rt_expr<InterfaceType> get() const { return arg_; }

      private:
        viennamath::rt_function_symbol<InterfaceType> const & u_;
        mutable viennamath::rt_expr<InterfaceType> arg_;
    };


    /** @brief A helper functor for updating the cell_quan tokens in a ViennaMath expression */
    template <typename NCellType, typename InterfaceType>
    struct ncell_updater : public viennamath::rt_traversal_interface<>
    {
      public:
        ncell_updater(NCellType const & ncell) : nc_(ncell) {}

        void operator()(InterfaceType const * e) const
        {
          if (viennamath::callback_if_castable< viennafvm::ncell_quan<NCellType, InterfaceType> >::apply(e, *this))
            return;
        }

        void operator()(viennafvm::ncell_quan<NCellType, InterfaceType> const & cq) const
        {
          cq.update(nc_);
          //std::cout << "cell_quan updated!" << std::endl;
        }

      private:
        NCellType const & nc_;
    };

  } //namespace detail


  template <typename CellType, typename FacetType, typename InterfaceType>
  class flux_handler
  {
    public:
      flux_handler(viennamath::rt_expr<InterfaceType> const & integrand, viennamath::rt_function_symbol<InterfaceType> const & u)
      {
        detail::gradient_scanner<InterfaceType> gradient_scanner(u);
        detail::func_symbol_scanner<InterfaceType> fsymbol_scanner(u);

        gradient_scanner(integrand.get());
        if (gradient_scanner.found())
        {
          std::cout << "Gradient found in flux expression!" << std::endl;
        }
        else
        {
          std::cout << "Gradient NOT found in flux expression: " << integrand << std::endl;
          throw "No gradient in flux!";
        }

        fsymbol_scanner(integrand.get());
        if (fsymbol_scanner.found())
        {
          std::cout << "Function symbol found!" << std::endl;
        }
        else
          std::cout << "Function symbol NOT found!" << std::endl;

        if (gradient_scanner.found() && fsymbol_scanner.found()) //advection-diffusion
        {
          std::cout << "Advection-Diffusion detected!" << std::endl;
        }
        else //pure diffusion
        {
          std::cout << "Purely diffusive problem detected!" << std::endl;

          // replace grad() by 1/distance in expression, where 1/distance is a cell quantity

          detail::gradient_argument_extractor<InterfaceType>  arg_extractor(u);
          detail::gradient_prefactor_extractor<InterfaceType> prefactor_extractor(u);

          arg_extractor(integrand.get());
          viennamath::rt_expr<InterfaceType> gradient_argument = arg_extractor.get();
          std::cout << "Gradient argument: " << gradient_argument << std::endl;

          prefactor_extractor(integrand.get());
          viennamath::rt_expr<InterfaceType> gradient_prefactor  = prefactor_extractor.get();
          std::cout << "Gradient prefactor: " << gradient_prefactor << std::endl;

          // Instantiate residual accessor for nonlinear terms:
          viennafvm::ncell_quan<CellType, InterfaceType> current_iterate;
          current_iterate.wrap_constant(viennafvm::current_iterate_key(u.id()));

          viennamath::rt_expr<InterfaceType> replaced_gradient_prefactor = viennamath::substitute(u, viennamath::rt_expr<InterfaceType>(current_iterate.clone()), gradient_prefactor);

          // replace grad() by (outer - inner) / distance.
          // Note that because of separate in() and out() member functions for the evaluation,
          // this is just the derivative of the gradient argument divided by the distance.
          //
          // Examples:
          //   # grad(u^2) leads to 2u/d for both inner and outer flux, where u is the current iterate.
          //   # grad(u)   leads to 1/d, presumably the most common case
          //
          // Todo: Deeper thought about stabilization in nonlinear case.
          //
          viennafvm::ncell_quan<FacetType, InterfaceType> distance;
          distance.wrap_constant(viennafvm::facet_distance_key());
          viennamath::rt_expr<InterfaceType> modified_gradient = viennamath::substitute(u, viennamath::rt_expr<InterfaceType>(current_iterate.clone()), viennamath::diff(gradient_argument, u)) / distance;

          in_integrand_  = replaced_gradient_prefactor * modified_gradient;
          out_integrand_ = replaced_gradient_prefactor * modified_gradient;

          std::cout << "Expression for in-flux:  " << in_integrand_ << std::endl;
          std::cout << "Expression for out-flux: " << out_integrand_ << std::endl;
        }

      }

      double in(CellType const & inner_cell, FacetType const & facet, CellType const & outer_cell) const
      {
        std::vector<double> p(3); //dummy point

        viennamath::rt_traversal_wrapper<InterfaceType> cell_updater( new detail::ncell_updater<CellType, InterfaceType>(inner_cell) );
        in_integrand_.get()->recursive_traversal(cell_updater);

        viennamath::rt_traversal_wrapper<InterfaceType> facet_updater( new detail::ncell_updater<FacetType, InterfaceType>(facet) );
        in_integrand_.get()->recursive_traversal(facet_updater);

        //return 1.0 / viennadata::access<viennafvm::facet_distance_key, double>()(facet);  //TODO: Unhack!
        return viennamath::eval(in_integrand_, p);
      }

      double out(CellType const & inner_cell, FacetType const & facet, CellType const & outer_cell) const
      {
        std::vector<double> p(3); //dummy point

        viennamath::rt_traversal_wrapper<InterfaceType> cell_updater( new detail::ncell_updater<CellType, InterfaceType>(outer_cell) );
        out_integrand_.get()->recursive_traversal(cell_updater);

        viennamath::rt_traversal_wrapper<InterfaceType> facet_updater( new detail::ncell_updater<FacetType, InterfaceType>(facet) );
        out_integrand_.get()->recursive_traversal(facet_updater);

        return viennamath::eval(out_integrand_, p);
      }

    private:
      viennamath::rt_expr<InterfaceType> in_integrand_;
      viennamath::rt_expr<InterfaceType> out_integrand_;
  };



} // end namespace viennafvm

#endif
