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

#ifndef VIENNAFVM_PDESTACK_HPP
#define VIENNAFVM_PDESTACK_HPP

// *** system includes
//
#include <vector>

// *** vienna includes
//
#include "viennamath/expression.hpp"

// *** boost includes
//
#include "boost/fusion/include/map.hpp"

namespace viennafvm {

namespace tag {
struct pde     {};
struct unknowns {};
}

struct pde_entry :  boost::fusion::map<
                           boost::fusion::pair<viennafvm::tag::pde,           viennamath::equation<>   >,
                           boost::fusion::pair<viennafvm::tag::unknowns,      std::vector< viennamath::function_symbol<> > > 
                        >  {};

struct pde_set : std::vector< pde_entry > {};

} // end namespace viennafvm

#endif

