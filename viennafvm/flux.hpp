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

#ifndef VIENNAFVM_FLUX_HPP
#define VIENNAFVM_FLUX_HPP

// *** system includes
//
#include <vector>

// *** local includes
//

// *** vienna includes
//
#include "viennamath/expression.hpp"
#include "viennamath/manipulation/substitute.hpp"

namespace viennafvm {

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
            std::cout << "Gradient extraction encountered unexpected expression: " << bin << std::endl;
            throw "Gradient extraction failed";
          }
        }

        bool found() const { return has_grad_; }
        void reset() const { has_grad_ = false; }

      private:
        viennamath::rt_function_symbol<InterfaceType> const & u_;
        mutable bool has_grad_;
    };

  } //namespace detail


  template <typename CellType, typename FacetType>
  class flux_handler
  {
    public:
      template <typename InterfaceType>
      flux_handler(viennamath::rt_expr<InterfaceType> integrand, viennamath::rt_function_symbol<InterfaceType> const & u)
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
        }
      }

      double in(CellType const & inner_cell, FacetType const & facet, CellType const & outer_cell) const
      {
        return 1.0 / viennadata::access<viennafvm::facet_distance_key, double>()(facet);  //TODO: Unhack!
      }

      double out(CellType const & inner_cell, FacetType const & facet, CellType const & outer_cell) const
      {
        return 1.0 / viennadata::access<viennafvm::facet_distance_key, double>()(facet);  //TODO: Unhack!
      }

    private:

  };



} // end namespace viennafvm

#endif

