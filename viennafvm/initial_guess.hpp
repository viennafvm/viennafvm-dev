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

#include "viennagrid/iterators.hpp"

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

  template <typename DomainType, typename FunctionSymbol, typename KeyType>
  void set_initial_guess(DomainType const & domain, FunctionSymbol const & u, KeyType const & key)
  {
    typedef viennafvm::current_iterate_key    IterateKey;
    typedef typename DomainType::config_type  ConfigType;
    typedef typename ConfigType::cell_tag     CellTag;

    typedef typename viennagrid::result_of::ncell<ConfigType, CellTag::dim>::type               CellType;
    typedef typename viennagrid::result_of::const_ncell_range<DomainType, CellTag::dim>::type   CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

    IterateKey iter_key(u.id());

    CellContainer cells = viennagrid::ncells(domain);
    for (CellIterator cit  = cells.begin();
                      cit != cells.end();
                    ++cit)
    {
      if (is_quantity_enabled(*cit, u.id()))
        viennadata::access<IterateKey, numeric_type>(iter_key)(*cit) = viennadata::access<KeyType, numeric_type>(key)(*cit);
    }
  }


  template <typename DomainType, typename FunctionSymbol, typename SmootherType>
  void smooth_initial_guess(DomainType const & domain, FunctionSymbol const & u, SmootherType const & smoother)
  {
    typedef viennafvm::current_iterate_key              IterateKey;
    typedef typename DomainType::config_type            Config;
    typedef typename Config::cell_tag                   CellTag;
    typedef typename viennagrid::result_of::ncell<Config, CellTag::dim-1>::type                FacetType;
    typedef typename viennagrid::result_of::ncell<Config, CellTag::dim  >::type                CellType;

    typedef typename viennagrid::result_of::const_ncell_range<DomainType, CellTag::dim>::type  CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type                      CellIterator;

    typedef typename viennagrid::result_of::const_ncell_range<CellType, CellTag::dim-1>::type  FacetOnCellContainer;
    typedef typename viennagrid::result_of::iterator<FacetOnCellContainer>::type               FacetOnCellIterator;

    std::map<CellType const *, std::vector<numeric_type> >  cell_neighbor_values;
    IterateKey iter_key(u.id());

    //
    // Phase 1: Gather neighboring values:
    //
    CellContainer cells = viennagrid::ncells(domain);
    for (CellIterator cit = cells.begin(); cit != cells.end(); ++cit)
    {
      if (is_quantity_disabled(*cit, u.id()))
        continue;

      cell_neighbor_values[&(*cit)].push_back(viennadata::access<IterateKey, numeric_type>(iter_key)(*cit));

      if (is_dirichlet_boundary(*cit, u)) // Dirichlet boundaries should not be smoothed
        continue;

      FacetOnCellContainer facets_on_cell = viennagrid::ncells(*cit);
      for (FacetOnCellIterator focit  = facets_on_cell.begin();
                               focit != facets_on_cell.end();
                             ++focit)
      {
        CellType const * other_cell = util::other_cell_of_facet(*focit, *cit, domain);

        if (other_cell)
        {
          if (is_quantity_disabled(*other_cell, u.id()))
            continue;

          numeric_type other_value = is_dirichlet_boundary(*other_cell, u)
              ? get_dirichlet_boundary(*other_cell, u)
              : viennadata::access<IterateKey, numeric_type>(iter_key)(*other_cell);

          cell_neighbor_values[&(*cit)].push_back(other_value);
        }
      }
    }

    //
    // Phase 2: Run averaging
    //
    for (CellIterator cit = cells.begin(); cit != cells.end(); ++cit)
    {
      if (is_quantity_disabled(*cit, u.id()))
        continue;

      viennadata::access<IterateKey, numeric_type>(iter_key)(*cit) = smoother(cell_neighbor_values[&(*cit)]);
    }
  }

}

#endif
