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

#include "viennagrid/forwards.hpp"
#include "viennagrid/mesh/segmentation.hpp"

#include <map>

/** @file  boundary.hpp
    @brief Provide convenience routines for setting boundary conditions
*/

namespace viennafvm
{

  /** @brief Assigns a value to the quantity of all cells of the segment or mesh */
  template <typename QuantityType, typename DomainSegmentType>
  void set_dirichlet_boundary(QuantityType & quan,
                              DomainSegmentType const & seg,
                              double value)
  {
    typedef typename viennagrid::result_of::cell_tag<DomainSegmentType>::type CellTag;

    typedef typename viennagrid::result_of::const_element_range<DomainSegmentType, CellTag>::type  CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

    CellContainer cells(seg);
    for (CellIterator cit  = cells.begin();
                      cit != cells.end();
                    ++cit)
    {
      quan.set_boundary_value(*cit, value);
      quan.set_boundary_type(*cit, viennafvm::BOUNDARY_DIRICHLET);
    }
  }

  /** @brief Adds a value to an already assigned quantity value on all cells of the segment or mesh */
  template <typename QuantityType, typename DomainSegmentType>
  void addto_dirichlet_boundary(QuantityType & quan,
                              DomainSegmentType const & seg,
                              double value)
  {
    typedef typename viennagrid::result_of::cell_tag<DomainSegmentType>::type CellTag;

    typedef typename viennagrid::result_of::const_element_range<DomainSegmentType, CellTag>::type  CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

    CellContainer cells(seg);
    for (CellIterator cit  = cells.begin();
                      cit != cells.end();
                    ++cit)
    {
      // if this cell has already been assigned a dirichlet value,
      // add the new value to the old one
      if(quan.get_boundary_type(*cit) == viennafvm::BOUNDARY_DIRICHLET)
      {
        quan.set_boundary_value(*cit, quan.get_boundary_value(*cit) + value);
      }
      // if this cell has infact not been assigend a dirichlet value,
      // simply assign the new value
      //
      else
      {
        quan.set_boundary_value(*cit, value);
        quan.set_boundary_type(*cit, viennafvm::BOUNDARY_DIRICHLET);
      }
    }
  }

  /** @brief Set the values of a quantity on all cells with a constant value */
  template <typename QuantityType, typename DomainSegmentType>
  void set_initial_value(QuantityType & quan,
                         DomainSegmentType const & seg,
                         double value)
  {
    typedef typename viennagrid::result_of::cell_tag<DomainSegmentType>::type CellTag;

    typedef typename viennagrid::result_of::const_element_range<DomainSegmentType, CellTag>::type  CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

    CellContainer cells(seg);
    for (CellIterator cit  = cells.begin();
                      cit != cells.end();
                    ++cit)
    {
      quan.set_value(*cit, value);
    }
  }

  /** @brief Set the values of a quantity on all cells according to an externally provided values container. 
             The keys of the 'values' container must correspond to the cell IDs. */
  template <typename QuantityType, typename DomainSegmentType>
  void set_initial_value(QuantityType                         & quan,
                         DomainSegmentType              const & seg,
                         std::map<std::size_t, double>  const & values)
  {
    typedef typename viennagrid::result_of::cell_tag<DomainSegmentType>::type CellTag;

    typedef typename viennagrid::result_of::const_element_range<DomainSegmentType, CellTag>::type  CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

    CellContainer cells(seg);
    for (CellIterator cit  = cells.begin();
                      cit != cells.end();
                    ++cit)
    {
      quan.set_value(*cit, values.at(cit->id().get()));
    }
  }

  /** @brief Transfer the quantity values from another quantity */
  template <typename QuantityType, typename DomainSegmentType>
  void set_initial_value(QuantityType & quan,
                         DomainSegmentType const & seg,
                         QuantityType & source)
  {
    typedef typename viennagrid::result_of::cell_tag<DomainSegmentType>::type CellTag;

    typedef typename viennagrid::result_of::const_element_range<DomainSegmentType, CellTag>::type  CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

    CellContainer cells(seg);
    for (CellIterator cit  = cells.begin();
                      cit != cells.end();
                    ++cit)
    {
      quan.set_value(*cit, source.get_value(*cit));
    }
  }

  /** @brief Assign the quantity based on an unary functor taking the cell id as input and returning the corresponding cell-based values */
  template <typename QuantityType, typename DomainSegmentType, typename FunctorT>
  void set_initial_value(QuantityType            & quan,
                         DomainSegmentType const & seg,
                         FunctorT                  functor)
  {
    set_initial_value(quan, seg, &functor);
  }

  /** @brief Assign the quantity based on a pointer to a unary functor taking the cell id as input and returning the corresponding cell-based values */
  template <typename QuantityType, typename DomainSegmentType, typename FunctorT>
  void set_initial_value(QuantityType            & quan,
                         DomainSegmentType const & seg,
                         FunctorT                * functor)
  {
    typedef typename viennagrid::result_of::cell_tag<DomainSegmentType>::type CellTag;

    typedef typename viennagrid::result_of::const_element_range<DomainSegmentType, CellTag>::type  CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

    CellContainer cells(seg);
    for (CellIterator cit  = cells.begin();
                      cit != cells.end();
                    ++cit)
    {
      quan.set_value(*cit, (*functor)(cit->id().get()));
    }
  }

  /** @brief Tag the quantity as a solution quantity */
  template <typename QuantityType, typename DomainSegmentType>
  void set_unknown(QuantityType & quan,
                   DomainSegmentType const & seg,
                   bool set_to = true)
  {
    typedef typename viennagrid::result_of::cell_tag<DomainSegmentType>::type CellTag;

    typedef typename viennagrid::result_of::const_element_range<DomainSegmentType, CellTag>::type  CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

    CellContainer cells(seg);
    for (CellIterator cit  = cells.begin();
                      cit != cells.end();
                    ++cit)
    {
      quan.set_unknown_mask(*cit, set_to);
    }
  }


}
#endif
