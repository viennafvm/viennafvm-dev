#ifndef VIENNAFVM_MAPPING_HPP
#define VIENNAFVM_MAPPING_HPP

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


#include <vector>
#include "viennafvm/forwards.h"
#include "viennafvm/common.hpp"

#include "viennagrid/forwards.hpp"

namespace viennafvm
{

/** @brief If EntityType is a ViennaGrid segment, returns the domain. If EntityType is already the domain, no changes.
*/
template <typename EntityType>
struct extract_domain
{
   typedef EntityType  type;
   static EntityType       & apply(EntityType       & mesh) { return mesh; }
   static EntityType const & apply(EntityType const & mesh) { return mesh; }
};

template <typename ConfigType>
struct extract_domain<viennagrid::segmentation<ConfigType> >
{
   typedef typename viennagrid::result_of::mesh<ConfigType>::type    type;
   static type       & apply(viennagrid::segment_handle<ConfigType>       & seg) { return seg.mesh(); }
   static type const & apply(viennagrid::segment_handle<ConfigType> const & seg) { return seg.mesh(); }
};



/** @brief Distributes mapping indices over domain or segment
*
*/
template <typename LinPdeSysT, typename DomainType, typename QuantityContainerT>
long create_mapping(LinPdeSysT & pde_system,
                    std::size_t  pde_index,
                    DomainType const & domain,
                    QuantityContainerT & quantities,
                    long start_index = 0)
{
  typedef typename viennagrid::result_of::cell<DomainType>::type CellTag;

  typedef typename viennagrid::result_of::point<DomainType>::type               PointType;
  typedef typename viennagrid::result_of::element<DomainType, CellTag>::type    CellType;

  typedef typename viennagrid::result_of::const_element_range<DomainType, CellTag>::type   CellContainer;
  typedef typename viennagrid::result_of::iterator<CellContainer>::type                    CellIterator;

  typedef typename QuantityContainerT::value_type       QuantityType;

  long map_index = start_index;

  long unknown_id = pde_system.unknown(pde_index)[0].id();

  QuantityType & quan = quantities.at(unknown_id);

  CellContainer cells(domain);
  for (CellIterator cit = cells.begin(); cit != cells.end(); ++cit)
  {
    if (quan.get_unknown_mask(*cit))   // quantity is set to unknown here, so it gets an unknown index assigned
    {
      quan.set_unknown_index(*cit, map_index);
      map_index += pde_system.unknown(pde_index).size();
    }
    else
      quan.set_unknown_index(*cit, -1);
  }

  return map_index;
}


/** @brief Distributes mapping indices over domain or segment
*
*/
template <typename LinPdeSysT, typename DomainType, typename QuantityContainerT>
long create_mapping(LinPdeSysT & pde_system,
                    DomainType const & domain,
                    QuantityContainerT & quantities)
{
  typedef typename viennagrid::result_of::cell<DomainType>::type CellTag;

  typedef typename viennagrid::result_of::point<DomainType>::type                  PointType;
  typedef typename viennagrid::result_of::element<DomainType, CellTag>::type    CellType;

  typedef typename viennagrid::result_of::const_element_range<DomainType, CellTag>::type   CellContainer;
  typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

  typedef typename LinPdeSysT::mapping_key_type   MappingKeyType;
  typedef typename LinPdeSysT::boundary_key_type  BoundaryKeyType;

  long next_index = 0;

  for (std::size_t pde_index = 0; pde_index < pde_system.size(); ++pde_index)
  {
    next_index = create_mapping(pde_system, pde_index, domain, quantities, next_index);
  }

  return next_index;
}


} // end namespace viennafvm

#endif

