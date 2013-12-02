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
#include <cmath>

// *** local includes
//

// *** vienna includes
//
#include "viennamath/expression.hpp"
#include "viennamath/manipulation/substitute.hpp"
#include "viennafvm/ncell_quantity.hpp"
#include "viennamath/manipulation/diff.hpp"
#include "viennamath/manipulation/eval.hpp"

#include "viennafvm/util.hpp"

#include "viennagrid/algorithm/interface.hpp"
#include "viennagrid/algorithm/centroid.hpp"

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
          if (fs.deep_equal(&u_)) //this is the same function symbol
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
              gradient_prefactor_extractor<InterfaceType> extractor(u_);
              extractor(bin.lhs());
              arg_ = viennamath::rt_expr<InterfaceType>(extractor.get().get()->clone());
              viennamath::inplace_simplify(arg_);
            }
            else
            {
              scanner->reset();

              bin.rhs()->recursive_traversal(wrapped_scanner);
              if (scanner->found())
              {
                gradient_prefactor_extractor<InterfaceType> extractor(u_);
                extractor(bin.rhs());
                if (dynamic_cast<const MinusOperatorType *>(bin.op()) != NULL)
                {
                  arg_ = viennamath::rt_binary_expr<InterfaceType>(new viennamath::rt_constant<NumericType, InterfaceType>(-1),
                                                                   new viennamath::op_binary<viennamath::op_mult<NumericType>, InterfaceType>(),
                                                                   extractor.get().get()->clone());
                  viennamath::inplace_simplify(arg_);
                }
                else
                {
                  arg_ = viennamath::rt_expr<InterfaceType>(extractor.get().get()->clone());
                  viennamath::inplace_simplify(arg_);
                }
              }
            }
          }
          else if (dynamic_cast<const ProductOperatorType *>(bin.op()) != NULL
                || dynamic_cast<const DivisionOperatorType *>(bin.op()) != NULL
                  )
          {
            gradient_scanner<InterfaceType> scanner(u_);

            scanner(bin.lhs());
            if (scanner.found())
            {
              gradient_prefactor_extractor<InterfaceType> extractor(u_);
              extractor(bin.lhs());
              arg_ = viennamath::rt_binary_expr<InterfaceType>(extractor.get().get()->clone(),
                                                               bin.op()->clone(),
                                                               bin.rhs()->clone());
              viennamath::inplace_simplify(arg_);
            }
            else
            {
              scanner.reset();

              scanner(bin.rhs());
              if (scanner.found())
              {
                gradient_prefactor_extractor<InterfaceType> extractor(u_);
                extractor(bin.rhs());
                arg_ = viennamath::rt_binary_expr<InterfaceType>(bin.lhs()->clone(),
                                                                 bin.op()->clone(),
                                                                 extractor.get().get()->clone());
                viennamath::inplace_simplify(arg_);
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


    /** @brief Scans for the factor in front of grad(u) */
    template <typename InterfaceType>
    struct func_symbol_prefactor_extractor : public viennamath::rt_traversal_interface<InterfaceType>
    {
      public:
        func_symbol_prefactor_extractor(viennamath::rt_function_symbol<InterfaceType> const & u) : u_(u), arg_() {}

        void operator()(InterfaceType const * e) const
        {
          typedef typename InterfaceType::numeric_type     NumericType;

          arg_ = viennamath::rt_constant<NumericType, InterfaceType>(0);
          viennamath::callback_if_castable< viennamath::rt_binary_expr<InterfaceType> >::apply(e, *this);
          viennamath::callback_if_castable< viennamath::rt_function_symbol<InterfaceType> >::apply(e, *this);
        }

        void operator()(viennamath::rt_function_symbol<InterfaceType> const & /*fs*/) const
        {
          typedef typename InterfaceType::numeric_type     NumericType;

          arg_ = viennamath::rt_constant<NumericType, InterfaceType>(1);
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
            func_symbol_scanner<InterfaceType> fs_scanner(u_);

            fs_scanner(bin.lhs());
            if (fs_scanner.found())
            {
              func_symbol_prefactor_extractor<InterfaceType> extractor(u_);
              extractor(bin.lhs());
              arg_ = viennamath::rt_expr<InterfaceType>(extractor.get().get()->clone());
              viennamath::inplace_simplify(arg_);
            }
            else
            {
              fs_scanner.reset();

              fs_scanner(bin.rhs());
              if (fs_scanner.found())
              {
                func_symbol_prefactor_extractor<InterfaceType> extractor(u_);
                extractor(bin.rhs());
                if (dynamic_cast<const MinusOperatorType *>(bin.op()) != NULL)
                {
                  arg_ = viennamath::rt_binary_expr<InterfaceType>(new viennamath::rt_constant<NumericType, InterfaceType>(-1),
                                                                   new viennamath::op_binary<viennamath::op_mult<NumericType>, InterfaceType>(),
                                                                   extractor.get().get()->clone());
                  viennamath::inplace_simplify(arg_);
                }
                else
                {
                  arg_ = viennamath::rt_expr<InterfaceType>(extractor.get().get()->clone());
                  viennamath::inplace_simplify(arg_);
                }
              }
            }
          }
          else if (dynamic_cast<const ProductOperatorType *>(bin.op()) != NULL
                || dynamic_cast<const DivisionOperatorType *>(bin.op()) != NULL
                  )
          {
            func_symbol_scanner<InterfaceType> * scanner = new func_symbol_scanner<InterfaceType>(u_);
            viennamath::rt_traversal_wrapper<InterfaceType> wrapped_scanner(scanner);

            bin.lhs()->recursive_traversal(wrapped_scanner);
            if (scanner->found())
            {
              func_symbol_prefactor_extractor<InterfaceType> extractor(u_);
              extractor(bin.lhs());
              arg_ = viennamath::rt_binary_expr<InterfaceType>(extractor.get().get()->clone(),
                                                               bin.op()->clone(),
                                                               bin.rhs()->clone());
              viennamath::inplace_simplify(arg_);
            }
            else
            {
              scanner->reset();

              bin.rhs()->recursive_traversal(wrapped_scanner);
              if (scanner->found())
              {
                func_symbol_prefactor_extractor<InterfaceType> extractor(u_);
                extractor(bin.rhs());
                arg_ = viennamath::rt_binary_expr<InterfaceType>(bin.lhs()->clone(),
                                                                 bin.op()->clone(),
                                                                 extractor.get().get()->clone());
                viennamath::inplace_simplify(arg_);
              }
            }
          }
          else
          {
            std::cout << "Function symbol prefactor extraction encountered unexpected expression: " << bin << std::endl;
            throw "Function symbol prefactor extraction failed";
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
        ncell_updater(NCellType const & ncell) : nc_(&ncell), other_cell_(NULL) {}
        ncell_updater(NCellType const & c1, NCellType const & c2) : nc_(&c1), other_cell_(&c2) {}

        void operator()(InterfaceType const * e) const
        {
          if (viennamath::callback_if_castable< viennafvm::ncell_quantity<NCellType, InterfaceType> >::apply(e, *this))
            return;
        }

        void operator()(viennafvm::ncell_quantity<NCellType, InterfaceType> const & cq) const
        {
          if (other_cell_ != NULL)
            cq.update(*nc_, *other_cell_);
          else
            cq.update(*nc_);
          //std::cout << "cell_quan updated!" << std::endl;
        }

      private:
        NCellType const * nc_;
        NCellType const * other_cell_;
    };

  } //namespace detail


  template <typename QuantityContainerT, typename CellType, typename FacetType, typename InterfaceType>
  class flux_handler
  {
    public:
      flux_handler(QuantityContainerT const & quantities,
                   viennamath::rt_expr<InterfaceType> const & integrand,
                   viennamath::rt_function_symbol<InterfaceType> const & u) : quantities_(quantities), has_advection_(false)
      {
        detail::gradient_scanner<InterfaceType> gradient_scanner(u);
        detail::func_symbol_scanner<InterfaceType> fsymbol_scanner(u);

        gradient_scanner(integrand.get());
        if (gradient_scanner.found())
        {
          //std::cout << "Gradient found in flux expression!" << std::endl;
        }
        else
        {
          std::cout << "Gradient NOT found in flux expression: " << integrand << std::endl;
          throw "No gradient in flux!";
        }

        fsymbol_scanner(integrand.get());

        //
        // Now prepare expression
        //

        detail::gradient_argument_extractor<InterfaceType>  arg_extractor(u);
        detail::gradient_prefactor_extractor<InterfaceType> grad_prefactor_extractor(u);

        arg_extractor(integrand.get());
        viennamath::rt_expr<InterfaceType> gradient_argument = arg_extractor.get();
#ifdef VIENNAFVM_DEBUG
        std::cout << " - Gradient argument: " << gradient_argument << std::endl;
#endif
        grad_prefactor_extractor(integrand.get());
        viennamath::rt_expr<InterfaceType> gradient_prefactor  = grad_prefactor_extractor.get();
#ifdef VIENNAFVM_DEBUG
        std::cout << " - Gradient prefactor: " << gradient_prefactor << std::endl;
#endif
        // Instantiate residual accessor for nonlinear terms:
        viennafvm::ncell_quantity<CellType, InterfaceType> current_iterate;
        current_iterate.wrap_constant(quantities.at(u.id()), true);

        viennamath::rt_expr<InterfaceType> replaced_gradient_prefactor = viennamath::substitute(u,
                                                                                                viennamath::rt_expr<InterfaceType>(current_iterate.clone()),
                                                                                                gradient_prefactor); // Note: requires additional thoughts on whether this makes sense


        if (gradient_scanner.found() && fsymbol_scanner.found()) //advection-diffusion
        {
#ifdef VIENNAFVM_DEBUG
          std::cout << " - Detected type of equation: Advection-Diffusion" << std::endl;
#endif
          has_advection_ = true;

          // extract prefactor of 'u':
          detail::func_symbol_prefactor_extractor<InterfaceType> fs_prefactor_extractor(u);

          fs_prefactor_extractor(integrand.get());
          viennamath::rt_expr<InterfaceType> fs_prefactor  = fs_prefactor_extractor.get();
#ifdef VIENNAFVM_DEBUG
          std::cout << " - Prefactor of " << u << ": " << fs_prefactor << std::endl;
#endif

          //
          // Canonical form: A * grad(n) + B * n   along straight line of length d
          //
          // The ratio |B/(A*d)| determines the discretization. If close to zero, we use a standard differencing scheme
          //

          A_ = replaced_gradient_prefactor;
          B_ = fs_prefactor;

#ifdef VIENNAFVM_DEBUG
          std::cout << " - Expression for stabilization term B/A (without distance d): " << B_ / A_ << std::endl;
#endif
        }
        else //pure diffusion
        {
#ifdef VIENNAFVM_DEBUG
          std::cout << " - Detected type of equation: Purely Diffusive" << std::endl;
#endif

          // replace grad() by 1/distance in expression, where 1/distance is a cell quantity

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
          viennamath::rt_expr<InterfaceType> modified_gradient = viennamath::substitute(u,
                                                                                        viennamath::rt_constant<double, InterfaceType>(1.0),
                                                                                        //viennamath::diff(gradient_argument, u)) / distance;
                                                                                        gradient_argument);

          in_integrand_  = modified_gradient;
          out_integrand_ = modified_gradient;
          integrand_prefactor_ = replaced_gradient_prefactor;

#ifdef VIENNAFVM_DEBUG
          std::cout << " - Expression for in-flux:  " << in_integrand_ << std::endl;
          std::cout << " - Expression for out-flux: " << out_integrand_ << std::endl;
#endif
        }

      }

      double in(CellType const & inner_cell, FacetType const & /*facet*/, CellType const & outer_cell, double distance) const
      {
        std::vector<double> p(3); //dummy point

        viennamath::rt_traversal_wrapper<InterfaceType> cell_updater_inner( new detail::ncell_updater<CellType, InterfaceType>(inner_cell) );
        viennamath::rt_traversal_wrapper<InterfaceType> cell_updater_outer( new detail::ncell_updater<CellType, InterfaceType>(outer_cell) );
        viennamath::rt_traversal_wrapper<InterfaceType> facet_updater( new detail::ncell_updater<CellType, InterfaceType>(inner_cell, outer_cell) );

        if (has_advection_)
        {
          A_.get()->recursive_traversal(cell_updater_inner);
          A_.get()->recursive_traversal(facet_updater);
          B_.get()->recursive_traversal(cell_updater_inner);
          B_.get()->recursive_traversal(facet_updater);

          double val_A = viennamath::eval(A_, p);
          double val_B = viennamath::eval(B_, p);
          double exponent = val_B / (val_A / distance);

          if ( std::abs(exponent) > 0.01) // Actual tolerance is not critical - this is for stabilization purposes only
            return val_B / (std::exp(exponent) - 1);
          else
            return val_A / distance;  // Note: Can be obtained from tailor expansion of the equation above
        }

        // pure diffusion:
        in_integrand_.get()->recursive_traversal(cell_updater_inner);
        in_integrand_.get()->recursive_traversal(facet_updater);

        integrand_prefactor_.get()->recursive_traversal(facet_updater);
        integrand_prefactor_.get()->recursive_traversal(cell_updater_inner);
        double eps_inner = viennamath::eval(integrand_prefactor_, p);
        integrand_prefactor_.get()->recursive_traversal(cell_updater_outer);
        double eps_outer = viennamath::eval(integrand_prefactor_, p);

        return viennamath::eval(in_integrand_, p) / distance * 2.0 * eps_inner * eps_outer / (eps_inner + eps_outer);
      }

      double out(CellType const & inner_cell, FacetType const & /*facet*/, CellType const & outer_cell, double distance) const
      {
        std::vector<double> p(3); //dummy point

        viennamath::rt_traversal_wrapper<InterfaceType> cell_updater_inner( new detail::ncell_updater<CellType, InterfaceType>(inner_cell) );
        viennamath::rt_traversal_wrapper<InterfaceType> cell_updater_outer( new detail::ncell_updater<CellType, InterfaceType>(outer_cell) );
        viennamath::rt_traversal_wrapper<InterfaceType> facet_updater( new detail::ncell_updater<CellType, InterfaceType>(inner_cell, outer_cell) );

        if (has_advection_)
        {
          A_.get()->recursive_traversal(cell_updater_inner);
          A_.get()->recursive_traversal(facet_updater);
          B_.get()->recursive_traversal(cell_updater_inner);
          B_.get()->recursive_traversal(facet_updater);

          double val_A = viennamath::eval(A_, p);
          double val_B = viennamath::eval(B_, p);
          double exponent = val_B / (val_A / distance);

          if ( std::abs(exponent) > 0.01) // Actual tolerance is not critical - this is for stabilization purposes only
            return val_B / (1.0 - std::exp(-exponent));
          else
            return val_A / distance;  // Note: Can be obtained from tailor expansion of the equation above
        }

        // pure diffusion:
        out_integrand_.get()->recursive_traversal(facet_updater);
        out_integrand_.get()->recursive_traversal(cell_updater_inner);

        integrand_prefactor_.get()->recursive_traversal(facet_updater);
        integrand_prefactor_.get()->recursive_traversal(cell_updater_inner);
        double eps_inner = viennamath::eval(integrand_prefactor_, p);
        integrand_prefactor_.get()->recursive_traversal(cell_updater_outer);
        double eps_outer = viennamath::eval(integrand_prefactor_, p);

        return viennamath::eval(out_integrand_, p) / distance * 2.0 * eps_inner * eps_outer / (eps_inner + eps_outer);
      }

    private:

      QuantityContainerT const & quantities_;

      bool has_advection_;

      // diffusive case:
      viennamath::rt_expr<InterfaceType> in_integrand_;
      viennamath::rt_expr<InterfaceType> out_integrand_;
      viennamath::rt_expr<InterfaceType> integrand_prefactor_;

      // diffusion-advection:
      viennamath::rt_expr<InterfaceType> A_;
      viennamath::rt_expr<InterfaceType> B_;
  };



  //
  // Accessor for fluxes:
  //
  template <typename ProblemDescriptionT>
  class flux_accessor
  {
    typedef typename ProblemDescriptionT::mesh_type     MeshType;

    typedef typename viennagrid::result_of::facet<MeshType>::type       FacetType;
    typedef typename viennagrid::result_of::cell<MeshType>::type        CellType;

    typedef typename ProblemDescriptionT::quantity_container_type       QuantityContainerType;
    typedef viennamath::rt_expr<>::interface_type                       InterfaceType;

  public:

    flux_accessor(ProblemDescriptionT const & problem_description, viennamath::rt_expr<> const & flux_term, viennamath::rt_function_symbol<> const & unknown)
      : problem_description_(problem_description),
        flux_(problem_description_.quantities(),
              prepare_for_evaluation<CellType>(problem_description_.quantities(), flux_term, unknown),
              unknown),
        unknown_index_(unknown.id()) {}

    /** @brief Evaluates the global flux vector on a cell
      *
      * Since FVM uses a cell-centered approach, only the projections of the flux onto the facet normals are easily computable.
      * This function computes the global flux vector out of these projections.
      */
    std::vector<double> operator()(CellType const & cell) const
    {
      //
      // Algorithm: Iterate over first k facets (k ... spatial dimension), extract unit normals n_i,
      //            setup linear system of equations for the project, then solve it.
      //

      std::vector<double> flux(viennagrid::result_of::topologic_cell_dimension<MeshType>::value);

      if (flux.size() == 1) // 1d is trivial: Just take the flux from the facet in positive orientation
      {
        FacetType const & facet1 = viennagrid::facets(cell)[0];
        FacetType const & facet2 = viennagrid::facets(cell)[1];

        if (viennagrid::default_point_accessor(problem_description_.mesh())(facet1)[0] < viennagrid::default_point_accessor(problem_description_.mesh())(facet2)[0])
        {
          flux[0] = operator()(cell, facet2);
          return flux;
        }

        flux[0] = operator()(cell, facet1);
      }
      else // solve projection equations:
      {
        throw "To be implemented!";
      }

      return flux;

    }

    /** @brief Evaluates the normal projection of the flux density onto the facet out of the provided cell.
      *
      * For example:
      *    --------
      *   |        |
      *   |  Cell  | -->
      *   |        | Flux out of cell w.r.t. the facet on the right
      *    --------
      *
      *  @return The flux density projection. Multiply by facet volume if you want to have the total flux (rather than the density) through that facet.
      */
    double operator()(CellType const & cell, FacetType const & facet) const
    {
      typedef typename viennagrid::result_of::point<MeshType>::type     PointType;
      typedef typename ProblemDescriptionT::quantity_type               QuantityType;

      QuantityType const & quan = problem_description_.quantities()[unknown_index_];

      CellType const * outer_cell = util::other_cell_of_facet(facet, cell, problem_description_.mesh());

      if (outer_cell == NULL) // this is the mesh boundary, so no outflux here:
        return 0;

      long index_inner_cell = quan.get_unknown_index(cell);
      long index_outer_cell = quan.get_unknown_index(*outer_cell);

      if (index_inner_cell < 0 && index_outer_cell < 0) // We don't have any unknowns here, hence no flux
        return 0;

      // compute distance:
      PointType centroid_1         = viennagrid::centroid(cell);
      PointType centroid_2         = viennagrid::centroid(*outer_cell);
      PointType center_connection  = centroid_2 - centroid_1;
      double distance              = viennagrid::norm(center_connection);

      // compute flux (same way as in matrix assembly in order to obtain consistent results)
      double quan_outer = quan.get_value(*outer_cell);
      double quan_inner = quan.get_value(cell);

      return   flux_.out(cell, facet, *outer_cell, distance) * quan_outer
             - flux_.in (cell, facet, *outer_cell, distance) * quan_inner;
    }


  private:
    ProblemDescriptionT const & problem_description_;
    viennafvm::flux_handler<QuantityContainerType, CellType, FacetType, InterfaceType> flux_;
    long unknown_index_;
  };



  template <typename SegmentationT, typename FluxEvaluatorT>
  double flux_between_segments(SegmentationT const & seg_src, SegmentationT const & seg_dest, FluxEvaluatorT const & flux_evaluator)
  {
    typedef typename viennagrid::result_of::point<SegmentationT>::type                 PointType;

    typedef typename viennagrid::result_of::cell_tag<SegmentationT>::type                       CellTag;
    typedef typename viennagrid::result_of::const_element_range<SegmentationT, CellTag>::type   CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

    typedef typename viennagrid::result_of::cell<SegmentationT>::type                       CellType;
    typedef typename viennagrid::result_of::facet_tag<SegmentationT>::type                  FacetTag;
    typedef typename viennagrid::result_of::const_element_range<CellType, FacetTag>::type   FacetOnCellContainer;
    typedef typename viennagrid::result_of::iterator<FacetOnCellContainer>::type            FacetOnCellIterator;

    double total_flux = 0.0;

    CellContainer cells(seg_src);
    for (CellIterator cit = cells.begin(); cit != cells.end(); ++cit)
    {
      FacetOnCellContainer facets_on_cell(*cit);
      for (FacetOnCellIterator focit = facets_on_cell.begin(); focit != facets_on_cell.end(); ++focit)
      {
        if (viennagrid::is_interface(seg_src, seg_dest, *focit))
        {
          CellType const * other_cell = util::other_cell_of_facet(*focit, *cit, seg_src.parent().mesh());

          assert(other_cell != NULL && bool("Logic error: Interface facet only attached to one cell!"));

          PointType centroid_1         = viennagrid::centroid(*cit);
          PointType centroid_2         = viennagrid::centroid(*other_cell);
          PointType center_connection  = centroid_1 - centroid_2;
          PointType outer_normal       = util::unit_outer_normal(*focit, *cit, viennagrid::default_point_accessor(seg_src)); //note: consistent orientation of center_connection and outer_normal is important here!

          double center_connection_len = viennagrid::norm(center_connection);
          double effective_facet_ratio = viennagrid::inner_prod(center_connection, outer_normal) / center_connection_len;  // inner product of unit vectors
          double effective_facet_area  = viennagrid::volume(*focit) * effective_facet_ratio;

          total_flux += flux_evaluator(*cit, *focit) * effective_facet_area;
        }
      }
    }

    return total_flux;
  }


} // end namespace viennafvm

#endif

