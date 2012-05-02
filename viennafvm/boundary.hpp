#ifndef VIENNAFVM_BOUNDARY_HPP
#define VIENNAFVM_BOUNDARY_HPP

/* =========================================================================
   Copyright (c) 2012, Institute for Microelectronics,
                       Institute for Analysis and Scientific Computing,
                       TU Wien.
                             -----------------
               ViennaFVM - The Vienna Finite Volume Method Library
                             -----------------

   Author:     Karl Rupp                          rupp@iue.tuwien.ac.at
               (add your name here)

   License:    MIT (X11), see file LICENSE in the ViennaFEM base directory
============================================================================ */

#include "viennafvm/forwards.h"

#include "viennamath/forwards.h"
#include "viennamath/manipulation/substitute.hpp"
#include "viennamath/expression.hpp"
#include "viennadata/api.hpp"

/** @file  boundary.hpp
    @brief Provide convenience routines for setting boundary conditions
*/

namespace viennafvm
{
  template <typename VertexType>
  void set_dirichlet_boundary(VertexType const & v,
                              numeric_type const & value,
                              std::size_t id = 0)
  {
    typedef viennafvm::boundary_key      BoundaryKey;
    
    //set flag:
    viennadata::access<BoundaryKey, bool >(BoundaryKey(id))(v) = true; 
    
    //set data:
    viennadata::access<BoundaryKey, numeric_type >(BoundaryKey(id))(v) = value; 
  }

  template <typename VertexType, typename NumericT>
  void set_dirichlet_boundary(VertexType const & v,
                              std::vector<NumericT> const & value,
                              std::size_t id = 0)
  {
    typedef viennafvm::boundary_key      BoundaryKey;;
    
    //set flag:
    viennadata::access<BoundaryKey, bool >(BoundaryKey(id))(v) = true; 
    
    //set data:
    viennadata::access<BoundaryKey, std::vector<NumericT> >(BoundaryKey(id))(v) = value; 
  }

}
#endif
