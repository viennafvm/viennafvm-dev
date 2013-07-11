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

#include "viennagrid/domain/domain.hpp"
#include "viennagrid/point.hpp"

#include "viennamath/runtime/function_symbol.hpp"

#include "viennafvm/forwards.h"

namespace viennafvm
{

  template <typename SegmentationType, typename AccessorType>
  void disable_quantity(viennagrid::segment_t<SegmentationType> const & seg,
                        AccessorType accessor)
  {
    typedef viennagrid::segment_t<SegmentationType> SegmentType;
    typedef typename viennagrid::result_of::cell_tag<SegmentType>::type CellTag;

    typedef typename viennagrid::result_of::element<SegmentType, CellTag>::type               CellType;
    typedef typename viennagrid::result_of::const_element_range<SegmentType, CellTag>::type  CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

    CellContainer cells = viennagrid::elements(seg);
    for (CellIterator cit  = cells.begin();
                      cit != cells.end();
                    ++cit)
    {
      //set flag:
      accessor(*cit) = true;
    }
  }

  template <typename ElementTag, typename WrappedConfigType, typename AccessorType>
  bool is_quantity_disabled(viennagrid::element_t<ElementTag, WrappedConfigType> const & cell,
                            AccessorType const accessor)
  {
    return accessor(cell);
  }

  template <typename ElementTag, typename WrappedConfigType, typename AccessorType>
  bool is_quantity_enabled(viennagrid::element_t<ElementTag, WrappedConfigType> const & cell,
                            AccessorType const accessor)
  {
    return !is_quantity_disabled(cell, accessor);
  }


  template <typename CellType, typename AccessorType>
  numeric_type get_current_iterate(CellType const & cell, AccessorType const accessor)
  {
    return accessor(cell);
  }

  template <typename CellType, typename AccessorType>
  void set_current_iterate(CellType const & cell, AccessorType accessor, numeric_type value)
  {
    accessor(cell) = value;
  }

} //namespace viennashe

#endif
