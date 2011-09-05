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

#ifndef VIENNAFVM_RHSZERO_HPP
#define VIENNAFVM_RHSZERO_HPP

// *** vienna includes
//
#include "viennamath/expression.hpp"

namespace viennafvm {


template <typename InterfaceType>
viennamath::equation<InterfaceType> 
make_rhs_zero(viennamath::equation<InterfaceType> const & equ)
{
   return viennamath::equation<InterfaceType>( equ.lhs()-equ.rhs(), 0); 
}


} // end namespace viennafvm

#endif

