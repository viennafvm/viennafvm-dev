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
  template <typename CellType>
  void set_dirichlet_boundary(CellType const & c,
                              numeric_type const & value,
                              std::size_t id = 0)
  {
    typedef viennafvm::boundary_key      BoundaryKey;

    //set flag:
    viennadata::access<BoundaryKey, bool >(BoundaryKey(id))(c) = true;

    //set data:
    viennadata::access<BoundaryKey, numeric_type >(BoundaryKey(id))(c) = value;
  }

  template <typename CellType, typename NumericT>
  void set_dirichlet_boundary(CellType const & c,
                              std::vector<NumericT> const & value,
                              std::size_t id = 0)
  {
    typedef viennafvm::boundary_key      BoundaryKey;;

    //set flag:
    viennadata::access<BoundaryKey, bool >(BoundaryKey(id))(c) = true;

    //set data:
    viennadata::access<BoundaryKey, std::vector<NumericT> >(BoundaryKey(id))(c) = value;
  }

}
#endif