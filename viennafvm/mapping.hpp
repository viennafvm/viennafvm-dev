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

#include "viennagrid/forwards.h"
#include "viennagrid/iterators.hpp"

namespace viennafvm
{

/** @brief If EntityType is a ViennaGrid segment, returns the domain. If EntityType is already the domain, no changes.
*/
template <typename EntityType>
struct extract_domain
{
   typedef EntityType  type;
   static EntityType       & apply(EntityType       & domain) { return domain; }
   static EntityType const & apply(EntityType const & domain) { return domain; }
};

template <typename ConfigType>
struct extract_domain<viennagrid::segment_t<ConfigType> >
{
   typedef typename viennagrid::result_of::domain<ConfigType>::type    type;
   static type       & apply(viennagrid::segment_t<ConfigType>       & seg) { return seg.domain(); }
   static type const & apply(viennagrid::segment_t<ConfigType> const & seg) { return seg.domain(); }
};


/** @brief Distributes mapping indices over domain or segment
*
*/
template <typename LinPdeSysT, typename DomainType>
long create_mapping(LinPdeSysT & pde_system,
                    DomainType const & domain)
{
  typedef typename DomainType::config_type              Config;
  typedef typename Config::cell_tag                     CellTag;

  typedef typename viennagrid::result_of::point<Config>::type                  PointType;
  typedef typename viennagrid::result_of::ncell<Config, CellTag::dim>::type    CellType;

  typedef typename viennagrid::result_of::const_ncell_range<DomainType, CellTag::dim>::type   CellContainer;
  typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

  typedef typename LinPdeSysT::mapping_key_type   MappingKeyType;
  typedef typename LinPdeSysT::boundary_key_type  BoundaryKeyType;

  long map_index = 0;
  bool init_done = false;

  for (std::size_t pde_index = 0; pde_index < pde_system.size(); ++pde_index)
  {
    long unknown_id = pde_system.option(pde_index).data_id();

    BoundaryKeyType bnd_key(unknown_id);
    MappingKeyType  map_key(unknown_id);

    //eventually, map indices need to be set to invalid first:
    if (!init_done)
    {
      typedef typename extract_domain<DomainType>::type   TrueDomainType;
      typedef typename viennagrid::result_of::const_ncell_range<TrueDomainType, CellTag::dim>::type  DomainCellContainer;
      typedef typename viennagrid::result_of::iterator<DomainCellContainer>::type                    DomainCellIterator;

      DomainCellContainer cells = viennagrid::ncells<CellTag::dim>(extract_domain<DomainType>::apply(domain));
      for (DomainCellIterator cit  = cells.begin();
                              cit != cells.end();
                            ++cit)
      {
        viennadata::access<MappingKeyType, long>(map_key)(*cit) = viennafvm::QUANTITY_DISABLED; // some negative value
      }
      viennadata::access<MappingKeyType, bool>(map_key)(extract_domain<DomainType>::apply(domain)) = true;
    }

    CellContainer cells = viennagrid::ncells(domain);
    for (CellIterator cit = cells.begin(); cit != cells.end(); ++cit)
    {
      if (viennadata::access<BoundaryKeyType, bool>(bnd_key)(*cit))  // boundary cell
      {
        //std::cout << "boundary cell" << std::endl;
        viennadata::access<MappingKeyType, long>(map_key)(*cit) = viennafvm::DIRICHLET_BOUNDARY; // some negative value
      }
      else if (viennafvm::is_quantity_enabled(*cit, unknown_id))
      {
        //std::cout << "interior cell" << std::endl;
        if (viennadata::access<MappingKeyType, long>(map_key)(*cit) < 0) //only assign if no dof assigned yet
        {
          viennadata::access<MappingKeyType, long>(map_key)(*cit) = map_index;
          map_index += pde_system.unknown(pde_index).size();
        }
      }
    }

    viennadata::access<MappingKeyType, long>(map_key)(extract_domain<DomainType>::apply(domain)) = map_index;

  }

  return map_index;
}


} // end namespace viennafvm

#endif

