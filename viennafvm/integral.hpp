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



#ifndef VIENNAFVM_INTEGRAL_HPP
#define VIENNAFVM_INTEGRAL_HPP

// *** system includes
// *** local includes
// *** vienna includes
#include "viennamath/forwards.h"
#include "viennamath/runtime/unary_expression.hpp"
#include "viennamath/compiletime/unary_op_tags.hpp"
// *** boost includes

namespace viennafvm {

struct Omega         {};
struct PartialOmega  {};
struct symbolic_tag  {};


template <typename BoundaryTag, typename InterfaceType>
viennamath::expr<InterfaceType> 
integral(viennamath::expr<InterfaceType> const & integrand, viennafvm::symbolic_tag)
{
   return viennamath::expr<InterfaceType>(
      new viennamath::unary_expr<InterfaceType>(
         integrand.get()->clone(),
         new viennamath::op_unary<viennamath::op_symbolic_integration<typename InterfaceType::numeric_type, BoundaryTag>, InterfaceType>()
      )
   );
}

template <typename BoundaryTag, typename InterfaceType>
viennamath::expr<InterfaceType> 
integral(viennamath::binary_expr<InterfaceType> const & integrand, viennafvm::symbolic_tag)
{
   return viennamath::expr<InterfaceType>(
      new viennamath::unary_expr<InterfaceType>(
         integrand.clone(),
         new viennamath::op_unary<viennamath::op_symbolic_integration<typename InterfaceType::numeric_type, BoundaryTag>, InterfaceType>()
      )
   );
}

} // end namespace viennafvm

#endif
