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

#include <viennagrid/domain.hpp>
#include <viennagrid/iterators.hpp>
#include <viennagrid/point.hpp>

namespace viennafvm
{

    template <typename DeviceType, typename VectorType>
    void updateQuantity(DeviceType const & domain, VectorType & quantity, VectorType const & update)
    {
      for (unsigned long index_x = 0; index_x < domain.size_x(); ++index_x)
      {
        for (unsigned long index_y = 0; index_y < domain.size_y(); ++index_y)
        {
          if (domain.isDirichletBoundary(index_x, index_y))
              continue;

          quantity(domain.getDOF(index_x, index_y)) += update(domain.getDOF(index_x, index_y));
        }
      }
    }

} //namespace viennashe

#endif
