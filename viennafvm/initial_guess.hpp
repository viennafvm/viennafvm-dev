#ifndef VIENNAFVM_INITIAL_GUESS_HPP
#define VIENNAFVM_INITIAL_GUESS_HPP

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

#include <cmath>
#include <assert.h>

#include "viennafvm/forwards.h"
#include "viennafvm/common.hpp"
#include "viennafvm/util.hpp"
#include "viennafvm/boundary.hpp"

#include "viennagrid/mesh/mesh.hpp"

#include "viennadata/api.hpp"


/** @file   initial_guess.hpp
    @brief  Sets the initial guesses for quantities
*/

namespace viennafvm
{

  class arithmetic_mean_smoother
  {
    public:
      template <typename NumericT>
      NumericT operator()(std::vector<NumericT> const & values) const
      {
        NumericT result = 0;
        for (std::size_t i=0; i<values.size(); ++i)
          result += values[i];
        return result / values.size();
      }
  };

  class geometric_mean_smoother
  {
    public:
      template <typename NumericT>
      NumericT operator()(std::vector<NumericT> const & values) const
      {
        NumericT result = 1;
        for (std::size_t i=0; i<values.size(); ++i)
        {
          assert(values[i] >= 0 && bool("Quantity a in geometric smoother negative!"));
          result *= std::pow(values[i], 1.0 / values.size());
        }

        return result;
      }
  };

  template <typename DomainSegmentType, typename QuantityDisabledAccessorType, typename SourceAccessorType, typename DestinationAccessorType>
  void set_initial_guess(DomainSegmentType const & domain,
                         QuantityDisabledAccessorType const quantity_disabled_accessor,
                         SourceAccessorType const source_accessor,
                         DestinationAccessorType destination_accessor)
  {
    typedef typename viennagrid::result_of::cell_tag<DomainSegmentType>::type     CellTag;


    typedef typename viennagrid::result_of::element<DomainSegmentType, CellTag>::type               CellType;
    typedef typename viennagrid::result_of::const_element_range<DomainSegmentType, CellTag>::type   CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

    CellContainer cells(domain);
    for (CellIterator cit  = cells.begin();
                      cit != cells.end();
                    ++cit)
    {
      if (!quantity_disabled_accessor(*cit))
        destination_accessor(*cit) = source_accessor(*cit);
    }
  }


  template <typename DomainSegmentType, typename StorageType, typename QuantityDisabledKeyType, typename SourceKeyType, typename DestinationKeyType>
  void set_initial_guess(DomainSegmentType const & domain,
                         StorageType & storage,
                         QuantityDisabledKeyType const & quantity_disabled_key,
                         SourceKeyType const & source_key,
                         DestinationKeyType const & destination_key)
  {
    typedef typename viennagrid::result_of::cell<DomainSegmentType>::type     CellType;

    set_initial_guess(domain,
                      viennadata::make_accessor<QuantityDisabledKeyType, bool, CellType>(storage, quantity_disabled_key),
                      viennadata::make_accessor<SourceKeyType, numeric_type, CellType>(storage, source_key),
                      viennadata::make_accessor<DestinationKeyType, numeric_type, CellType>(storage, destination_key));
  }


  template <typename DomainSegmentType, typename StorageType, typename InterfaceType, typename SourceKeyType>
  void set_initial_guess(DomainSegmentType const & domain,
                         StorageType & storage,
                         viennamath::rt_function_symbol<InterfaceType> const & func_symbol,
                         SourceKeyType const & source_key)
  {
    set_initial_guess(domain, storage,
                      viennafvm::disable_quantity_key(func_symbol.id()),
                      source_key,
                      viennafvm::current_iterate_key(func_symbol.id()));
  }




  template <typename DomainSegmentType, typename QuantityType, typename SmootherType>
  void smooth_initial_guess(DomainSegmentType const & domseg,
                            QuantityType & quan,
                            SmootherType const & smoother)
  {
    typedef viennafvm::current_iterate_key              IterateKey;

    typedef typename viennagrid::result_of::cell_tag<DomainSegmentType>::type     CellTag;
    typedef typename viennagrid::result_of::facet_tag<CellTag>::type     FacetTag;

    typedef typename viennagrid::result_of::element<DomainSegmentType, FacetTag>::type                FacetType;
    typedef typename viennagrid::result_of::element<DomainSegmentType, CellTag >::type                CellType;

    typedef typename viennagrid::result_of::const_element_range<DomainSegmentType, CellTag>::type  CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type                      CellIterator;

    typedef typename viennagrid::result_of::const_element_range<CellType, FacetTag>::type  FacetOnCellContainer;
    typedef typename viennagrid::result_of::iterator<FacetOnCellContainer>::type               FacetOnCellIterator;

    std::map<CellType const *, std::vector<numeric_type> >  cell_neighbor_values;

    //
    // Phase 1: Gather neighboring values:
    //
    CellContainer cells(domseg);
    for (CellIterator cit = cells.begin(); cit != cells.end(); ++cit)
    {
      cell_neighbor_values[&(*cit)].push_back(quan.get_value(*cit));

      if (quan.get_boundary_type(*cit) == viennafvm::BOUNDARY_DIRICHLET) // Dirichlet boundaries should not be smoothed
        continue;

      FacetOnCellContainer facets_on_cell(*cit);
      for (FacetOnCellIterator focit  = facets_on_cell.begin();
                               focit != facets_on_cell.end();
                             ++focit)
      {
        CellType const * other_cell = util::other_cell_of_facet(*focit, *cit, domseg);

        if (other_cell)
          cell_neighbor_values[&(*cit)].push_back(quan.get_value(*other_cell));
      }
    }

    //
    // Phase 2: Run averaging
    //
    for (CellIterator cit = cells.begin(); cit != cells.end(); ++cit)
    {
      quan.set_value(*cit, smoother(cell_neighbor_values[&(*cit)]));
    }
  }


}

#endif
