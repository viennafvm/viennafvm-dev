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

#ifndef VIENNAFVM_EXTRACTINTEGRALS_HPP
#define VIENNAFVM_EXTRACTINTEGRALS_HPP

// *** system includes
//
#include <vector>

// *** local includes
//

// *** vienna includes
//
#include "viennamath/expression.hpp"
#include "viennamath/manipulation/substitute.hpp"
#include "viennamath/manipulation/integral.hpp"

namespace viennafvm {

namespace detail {

template <typename InterfaceType, typename IntegratorTag>
struct integral_visitor : public viennamath::traversal_interface<InterfaceType>
{
   void operator()(InterfaceType const * e) const 
   {
      std::cout << "  * visitor: " << viennamath::expr<InterfaceType>(e->clone()) << std::endl;

//      if (viennamath::callback_if_castable< IntegratorTag >::apply(e, *this)) 
//      {
//      }

   }
};

} // end namespace detail


struct extract_integrals
{
   typedef viennamath::expr<>          expr_type;
   typedef std::vector<expr_type>      result_type;
   
   
   template <typename InterfaceType, typename IntegratorTag>
   static result_type
   eval(viennamath::expr<InterfaceType> const & expr, IntegratorTag const&)
   {
      std::cout << "Extracting partial omega integral terms from: " << std::endl;
      std::cout << expr << std::endl;
   
      viennamath::traversal_wrapper<InterfaceType>    visitor( new viennafvm::detail::integral_visitor<InterfaceType, IntegratorTag>() );
   
      expr.get()->recursive_traversal(visitor);   
   
      result_type result;
      return result;
   }         
   
//   template <typename InterfaceType, typename BoundaryTag>
//   static result_type
//   eval(viennamath::expr<InterfaceType> const & expr, 
//        viennamath::op_symbolic_integration<typename InterfaceType::numeric_type, BoundaryTag> const&)
//   {
//      result_type result;
//      return result;
//   }
//   
//   template <typename InterfaceType>
//   static result_type
//   eval(viennamath::expr<InterfaceType> const & expr, 
//        viennamath::op_symbolic_integration<typename InterfaceType::numeric_type, viennafvm::Omega> const&)
//   {
//      result_type result;
//      return result;
//   }
//   
//   template <typename InterfaceType>
//   static result_type
//   eval(viennamath::expr<InterfaceType> const & expr, 
//        viennamath::op_symbolic_integration<typename InterfaceType::numeric_type, viennafvm::PartialOmega> const&)
//   {
//      std::cout << "Extracting partial omega integral terms from: " << std::endl;
//      std::cout << expr << std::endl;
//   
//      viennamath::traversal_wrapper<InterfaceType>    visitor( new viennafvm::detail::integral_visitor<InterfaceType>() );
//   
//      expr.get()->recursive_traversal(visitor);   
//   
//      result_type result;
//      return result;
//   }      
};




} // end namespace viennafvm

#endif

