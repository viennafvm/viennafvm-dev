#ifndef VIENNAFVM_COMMON_HPP
#define VIENNAFVM_COMMON_HPP

/* =======================================================================
   Copyright (c) 2011, Institute for Microelectronics, TU Wien
   http://www.iue.tuwien.ac.at
                             -----------------
           ViennaFVM - The Vienna Finite Volume Method Library
                             -----------------

   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at
               (add your name here)

   license:    To be discussed, see file LICENSE in the ViennaFVM base directory
======================================================================= */

#include <iostream>
#include <fstream>
#include <vector>

#include <math.h>

#include "viennagrid/domain.hpp"
#include "viennagrid/iterators.hpp"
#include "viennagrid/point.hpp"

#include "viennamath/runtime/function_symbol.hpp"

#include "viennafvm/forwards.h"

namespace viennafvm
{

  template <typename ConfigType, typename InterfaceType>
  void disable_quantity(viennagrid::segment_t<ConfigType> const & seg,
                        viennamath::rt_function_symbol<InterfaceType> const & fs)
  {
    typedef viennafvm::disable_quantity_key   DisablerKey;
    typedef viennagrid::segment_t<ConfigType> SegmentType;
    typedef typename ConfigType::cell_tag     CellTag;

    typedef typename viennagrid::result_of::ncell<ConfigType, CellTag::dim>::type               CellType;
    typedef typename viennagrid::result_of::const_ncell_range<SegmentType, CellTag::dim>::type  CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

    DisablerKey key(fs.id());

    CellContainer cells = viennagrid::ncells(seg);
    for (CellIterator cit  = cells.begin();
                      cit != cells.end();
                    ++cit)
    {
      //set flag:
      viennadata::access<DisablerKey, bool>(key)(*cit) = true;
    }
  }

  template <typename ConfigType, typename ElementTag>
  bool is_quantity_disabled(viennagrid::element_t<ConfigType, ElementTag> const & cell,
                            long id)
  {
    typedef viennafvm::disable_quantity_key   DisablerKey;

    DisablerKey key(id);

    return viennadata::access<DisablerKey, bool>(key)(cell);
  }

  template <typename ConfigType, typename ElementTag>
  bool is_quantity_enabled(viennagrid::element_t<ConfigType, ElementTag> const & cell,
                           long id)
  {
    return !is_quantity_disabled(cell, id);
  }

} //namespace viennashe

#endif
