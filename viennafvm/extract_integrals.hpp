#ifndef VIENNAFVM_EXTRACTINTEGRALS_HPP
#define VIENNAFVM_EXTRACTINTEGRALS_HPP

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
#include "viennafvm/ncell_quantity.hpp"

namespace viennafvm {




  namespace detail
  {
    /** @brief Transforms a strong formulation of an equation to a weak form, assuming homogeneous Neumann boundary conditions */
    template <typename InterfaceType>
    struct integrand_extractor : public viennamath::rt_manipulation_interface<InterfaceType>
    {
      public:
        integrand_extractor(viennamath::id_type interval_id = 0) : id_(interval_id) {}

        InterfaceType * operator()(InterfaceType const * e) const
        {
          if (   !viennamath::callback_if_castable< viennamath::rt_unary_expr<InterfaceType> >::apply(e, *this)
              && !viennamath::callback_if_castable< viennamath::rt_binary_expr<InterfaceType> >::apply(e, *this))
          {
            std::cout << "Integrand extraction stalled at e=" << e->deep_str() << std::endl;
            throw "Cannot derive weak form!";
          }

          return integrated_expr.get()->clone();
        }

        void operator()(viennamath::rt_unary_expr<InterfaceType> const & unary_expr) const
        {
          typedef typename InterfaceType::numeric_type   NumericType;
          typedef viennamath::op_unary<viennamath::op_rt_symbolic_integral<InterfaceType>, InterfaceType>  SymbolicIntegralOperator;

          const SymbolicIntegralOperator * symb_op = dynamic_cast<const SymbolicIntegralOperator *>(unary_expr.op());
          if (symb_op != NULL) //this is a surface integral
          {
            if (symb_op->op().interval().id() == id_)
              integrated_expr = viennamath::rt_expr<InterfaceType>(unary_expr.lhs()->clone());
            else
              integrated_expr = viennamath::rt_constant<NumericType, InterfaceType>(0);
          }
          else
            throw "Cannot extract integrand!";
        }

        void operator()(viennamath::rt_binary_expr<InterfaceType> const & bin) const
        {
          typedef typename InterfaceType::numeric_type   NumericType;
          typedef viennamath::op_binary<viennamath::op_plus<NumericType>, InterfaceType>   PlusOperatorType;
          typedef viennamath::op_binary<viennamath::op_minus<NumericType>, InterfaceType>  MinusOperatorType;

          if (    dynamic_cast<const PlusOperatorType *>(bin.op()) != NULL
              || dynamic_cast<const MinusOperatorType *>(bin.op()) != NULL) //integration is additive :-)
          {
            viennamath::rt_manipulation_wrapper<InterfaceType> manipulator(new integrand_extractor<InterfaceType>(id_));
            //Note: In the following, directly passing *this is not possible due to the need for a wrapper...
            integrated_expr = new viennamath::rt_binary_expr<InterfaceType>(bin.lhs()->recursive_manipulation(manipulator),
                                                                            bin.op()->clone(),
                                                                            bin.rhs()->recursive_manipulation(manipulator));
          }
          else //TODO: Add checks!
          {
            throw "Integrand extraction failed";
          }
        }

        bool modifies(InterfaceType const * /*e*/) const { return true; }

      private:
        viennamath::id_type id_;
        mutable viennamath::rt_expr<InterfaceType> integrated_expr;
    };



    template <typename QuantityContainerT, typename CellType, typename InterfaceType>
    struct function_symbol_replacer : public viennamath::rt_manipulation_interface<InterfaceType>
    {
      public:
        function_symbol_replacer(QuantityContainerT const & quantities, viennamath::rt_function_symbol<InterfaceType> u) : quantities_(quantities), u_(u) {}

        InterfaceType * operator()(InterfaceType const * e) const
        {
          manipulated_expression_ = viennamath::rt_expr<InterfaceType>(e->clone());

          viennamath::callback_if_castable< viennamath::rt_function_symbol<InterfaceType> >::apply(e, *this);
          viennamath::callback_if_castable< viennamath::rt_unary_expr<InterfaceType> >::apply(e, *this);
          viennamath::callback_if_castable< viennamath::rt_binary_expr<InterfaceType> >::apply(e, *this);

          return manipulated_expression_.get()->clone();
        }

        void operator()(viennamath::rt_function_symbol<InterfaceType> const & func_symb) const
        {
          if (func_symb.id() != u_.id())
          {
            viennafvm::ncell_quantity<CellType, InterfaceType> current_iterate;
            current_iterate.wrap_constant(quantities_.at(func_symb.id()), false);

            manipulated_expression_ = viennamath::rt_expr<InterfaceType>(current_iterate.clone());
          }

        }

        void operator()(viennamath::rt_unary_expr<InterfaceType> const & unary_expr) const
        {
          typedef typename InterfaceType::numeric_type   NumericType;
          typedef viennamath::op_unary<viennamath::op_gradient<NumericType>, InterfaceType>  GradientOperator;

          const GradientOperator * grad_op = dynamic_cast<const GradientOperator *>(unary_expr.op());
          if (grad_op != NULL) //this is a gradient
          {
            typedef viennamath::rt_function_symbol<InterfaceType> FunctionSymbol;

            const FunctionSymbol * fs = dynamic_cast<const FunctionSymbol *>(unary_expr.lhs());
            if (fs != NULL)
            {
              if (fs->id() != u_.id())
              {
                viennafvm::ncell_quantity<CellType, InterfaceType> current_iterate;
                current_iterate.wrap_constant(quantities_.at(fs->id()), true);

                //manipulated_expression_ = viennamath::rt_unary_expr<InterfaceType>(current_iterate.clone(), grad_op->clone());
                manipulated_expression_ = viennamath::rt_expr<InterfaceType>(current_iterate.clone()); // Shortcut: Gradient values are assumed to be accessible at facets
              }
            }
          }
          else
            throw "Cannot derive weak form!";
        }

        void operator()(viennamath::rt_binary_expr<InterfaceType> const & bin) const
        {
          typedef typename InterfaceType::numeric_type   NumericType;
          typedef viennamath::op_binary<viennamath::op_plus<NumericType>, InterfaceType>   PlusOperatorType;
          typedef viennamath::op_binary<viennamath::op_minus<NumericType>, InterfaceType>  MinusOperatorType;
          typedef viennamath::op_binary<viennamath::op_mult<NumericType>, InterfaceType>   MultOperatorType;
          typedef viennamath::op_binary<viennamath::op_div<NumericType>, InterfaceType>    DivOperatorType;

          if (   dynamic_cast<const PlusOperatorType *>(bin.op()) != NULL
              || dynamic_cast<const MinusOperatorType *>(bin.op()) != NULL
              || dynamic_cast<const MultOperatorType *>(bin.op()) != NULL
              || dynamic_cast<const DivOperatorType *>(bin.op()) != NULL)
          {
            viennamath::rt_manipulation_wrapper<InterfaceType> manipulator(new function_symbol_replacer<QuantityContainerT, CellType, InterfaceType>(quantities_, u_));
            //Note: In the following, directly passing *this is not possible due to the need for a wrapper...
            manipulated_expression_ = new viennamath::rt_binary_expr<InterfaceType>(bin.lhs()->recursive_manipulation(manipulator),
                                                                                    bin.op()->clone(),
                                                                                    bin.rhs()->recursive_manipulation(manipulator));
          }
          else //TODO: Add checks!
          {
            throw "Integrand extraction failed";
          }
        }

        bool modifies(InterfaceType const * /*e*/) const { return true; }

      private:
        QuantityContainerT const & quantities_;
        viennamath::rt_function_symbol<InterfaceType> u_;
        mutable viennamath::rt_expr<InterfaceType> manipulated_expression_;
    };

  } //namespace detail


  /** @brief Replaces all function symbols with their respective cell accessors */
  template <typename CellType, typename QuantityContainerT, typename InterfaceType>
  viennamath::rt_expr<InterfaceType> prepare_for_evaluation(QuantityContainerT & quantities, viennamath::rt_expr<InterfaceType> const & ex, viennamath::rt_function_symbol<InterfaceType> const & u)
  {
    // replace 'other' function symbols (not the same as u) with ncell_quantity:
    viennamath::rt_manipulation_wrapper<InterfaceType> wrapped_function_symbol_replacer( new detail::function_symbol_replacer<QuantityContainerT, CellType, InterfaceType>(quantities, u) );
    viennamath::rt_expr<InterfaceType> replaced_ex(ex.get()->recursive_manipulation( wrapped_function_symbol_replacer ));

    return viennamath::simplify(replaced_ex);
  }


  /** @brief Discards all terms other than surface integrals from an expression */
  template <typename CellType, typename QuantityContainerT, typename InterfaceType>
  viennamath::rt_expr<InterfaceType> extract_surface_integrand(QuantityContainerT & quantities,
                                                               viennamath::rt_expr<InterfaceType> const & ex,
                                                               viennamath::rt_function_symbol<InterfaceType> const & u)
  {
    viennamath::rt_manipulation_wrapper<InterfaceType> wrapped_surface_integrand_extractor( new detail::integrand_extractor<InterfaceType>(viennafvm::cell_boundary) );
    viennamath::rt_expr<InterfaceType> integrand(ex.get()->recursive_manipulation( wrapped_surface_integrand_extractor ));

    return prepare_for_evaluation<CellType>(quantities, integrand, u);
  }


  /** @brief Discards all terms other than volume integrals from an expression */
  template <typename CellType, typename QuantityContainerT, typename InterfaceType>
  viennamath::rt_expr<InterfaceType> extract_volume_integrand(QuantityContainerT & quantities, viennamath::rt_expr<InterfaceType> const & ex, viennamath::rt_function_symbol<InterfaceType> const & u)
  {
    viennamath::rt_manipulation_wrapper<InterfaceType> wrapped_volume_integrand_extractor( new detail::integrand_extractor<InterfaceType>(viennafvm::cell_volume) );
    viennamath::rt_expr<InterfaceType> integrand(ex.get()->recursive_manipulation( wrapped_volume_integrand_extractor ));

    return prepare_for_evaluation<CellType>(quantities, integrand, u);
  }


} // end namespace viennafvm

#endif

