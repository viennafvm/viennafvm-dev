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

#ifndef VIENNAFVM_WEAKFORM_HPP
#define VIENNAFVM_WEAKFORM_HPP

// *** local includes
//
#include "viennafvm/integral.hpp"

// *** vienna includes
//
#include "viennamath/expression.hpp"
#include "viennamath/manipulation/substitute.hpp"
#include "viennamath/manipulation/integral.hpp"

namespace viennafvm {

namespace detail {

template <typename InterfaceType>
struct weak_form_creator : public viennamath::manipulation_interface<InterfaceType>
{
   // this is the general case: the recursive term traversal always starts from 
   // the top and progresses forward. this functor is always called at the beginning
   //
   InterfaceType * operator()(InterfaceType const * e) const 
   {
      //std::cout << "weak-form-creator::general: " << viennamath::expr<InterfaceType>(e->clone()) << std::endl;
   
      // if trivial expressions such as 'u' (L^2-projection)
      //
      if( !viennamath::callback_if_castable< viennamath::unary_expr<InterfaceType> >::apply(e, *this) &&
          !viennamath::callback_if_castable< viennamath::binary_expr<InterfaceType> >::apply(e, *this))
      {
         // integrate
         //
         viennamath::expr<InterfaceType> temp(e->clone());
         integrated_expr = viennafvm::integral<viennafvm::Omega>(temp, viennafvm::symbolic_tag());
      }

      // otherwise forward to the specializations: unary/binary
      //
      return integrated_expr.get()->clone();
   }

   void operator()(viennamath::unary_expr<InterfaceType> const & expr) const
   {
      //std::cout << "weak-form-creator::unary-expr: " << expr << std::endl;
      
      typedef typename InterfaceType::numeric_type   NumericType;
      typedef viennamath::op_unary<viennamath::op_divergence<NumericType>, InterfaceType>  DivergenceOperatorType;

      // if expression is of the form div(expression) 
      //
      if (dynamic_cast<const DivergenceOperatorType *>(expr.op()) != NULL) 
      {
         // replace by a dOmega integral (Gauss theorem ..)
         //
         viennamath::expr<InterfaceType>  lhs(expr.lhs()->clone());
         integrated_expr = viennafvm::integral<viennafvm::PartialOmega>(lhs, viennafvm::symbolic_tag());
      }
      else throw "# ViennaFVM::WeakFormCreator: unary expression not implemented .. ";
   }
   
   void operator()(viennamath::binary_expr<InterfaceType> const & expr) const
   {
      throw "# ViennaFVM::WeakFormCreator: binary expression handling not implemented .. ";
   }
   
   // [JW] this function is important, as we need to tell the recursive caller
   // that we are modifiable. otherwise, the bottom level object is called, instead 
   // of the required top-level one
   //
   bool modifies(InterfaceType const * e) const { return true; }
      
   
private:
   mutable viennamath::expr<InterfaceType>   integrated_expr;   
};

} // end namespace detail




template <typename InterfaceType>
viennamath::equation<InterfaceType> 
make_weak_form(viennamath::equation<InterfaceType> const & strong_formulation)
{
   // setup a functor which analyses and transforms specific terms of an equation
   //
   viennamath::manipulation_wrapper<InterfaceType>    weak_former( new viennafvm::detail::weak_form_creator<InterfaceType>() );

   // recursively traverse the terms of the LHS and apply the weak formulation functor on each of the terms
   // 
   viennamath::expr<InterfaceType> weak_lhs(strong_formulation.lhs().get()->recursive_manipulation(weak_former));

   // recursively traverse the terms of the RHS and apply the weak formulation functor on each of the terms
   //
   viennamath::expr<InterfaceType> weak_rhs(strong_formulation.rhs().get()->recursive_manipulation(weak_former));

   // build the equation from the weak forms of the LHS and RHS and return it
   //
   return viennamath::equation<InterfaceType>( weak_lhs, weak_rhs);
}


} // end namespace viennafvm

#endif

