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
#include "viennagrid/domain/segmentation.hpp"

/** @file  boundary.hpp
    @brief Provide convenience routines for setting boundary conditions
*/

namespace viennafvm
{
//   template <typename CellType, typename BoundaryAccessorType, typename BoundaryValueAccessorType>
//   void set_dirichlet_boundary(CellType const & c,
//                               numeric_type const & value,
//                               BoundaryAccessorType boundary_accessor,
//                               BoundaryValueAccessorType boundary_value_accessor)
//   {
// //     typedef viennafvm::boundary_key      BoundaryKey;
// 
//     //set flag:
//     boundary_accessor(c) = true;
// //     viennadata::access<BoundaryKey, bool >(BoundaryKey(id))(c) = true;
// 
//     //set data:
//     boundary_value_accessor(c) = value;
// //     viennadata::access<BoundaryKey, numeric_type >(BoundaryKey(id))(c) = value;
//   }

  template <typename DomainSegmentType, typename BoundaryAccessorType, typename BoundaryValueAccessorType>
  void set_dirichlet_boundary(DomainSegmentType const & seg,
                              BoundaryAccessorType boundary_accessor,
                              BoundaryValueAccessorType boundary_value_accessor,
                              numeric_type const & value)
  {
    typedef typename viennagrid::result_of::cell_tag<DomainSegmentType>::type CellTag;

    typedef typename viennagrid::result_of::element<DomainSegmentType, CellTag>::type               CellType;
    typedef typename viennagrid::result_of::const_element_range<DomainSegmentType, CellTag>::type  CellContainer;
    typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

    CellContainer cells = viennagrid::elements(seg);
    for (CellIterator cit  = cells.begin();
                      cit != cells.end();
                    ++cit)
    {
      //set flag:
      boundary_accessor(*cit) = true;

      //set data:
      boundary_value_accessor(*cit) = value;
    }
  }

  
  template <typename DomainSegmentType, typename A, typename B, typename KeyType>
  void set_dirichlet_boundary(DomainSegmentType const & domseg,
                              viennadata::storage<A,B> & storage,
                              KeyType const & key,
                              numeric_type const & value)
  {
    typedef typename viennagrid::result_of::cell< DomainSegmentType >::type CellType;
    set_dirichlet_boundary(domseg,
                           viennadata::accessor<KeyType, bool, CellType>(storage, key),
                           viennadata::accessor<KeyType, numeric_type, CellType>(storage, key),
                           value);
  }

  
//   template <typename SegmentationType, typename KeyType, typename StorageType>
//   void set_dirichlet_boundary(viennagrid::segment_t<SegmentationType> const & seg,
//                               StorageType & storage,
//                               KeyType const & key,
//                               numeric_type const & value,
//                               long id)
//   {
//     typedef typename viennagrid::result_of::cell< viennagrid::segment_t<SegmentationType> >::type CellType;
// 
//     typename viennadata::result_of::accessor<StorageType, viennafvm::boundary_key, bool, CellType>::type boundary_accessor =
//         viennadata::accessor<viennafvm::boundary_key, bool, CellType>(storage, viennafvm::boundary_key(id));
// 
//     typename viennadata::result_of::accessor<StorageType, viennafvm::boundary_key, numeric_type, CellType>::type boundary_value_accessor =
//         viennadata::accessor<viennafvm::boundary_key, numeric_type, CellType>(storage, viennafvm::boundary_key(id));
// 
//     set_dirichlet_boundary(seg, value, boundary_accessor, boundary_value_accessor);
//   }
  
  template <typename DomainSegmentType, typename A, typename B, typename InterfaceType>
  void set_dirichlet_boundary(DomainSegmentType const & domseg,
                              viennadata::storage<A,B> & storage,
                              viennamath::rt_function_symbol<InterfaceType> const & func_symbol,
                              numeric_type const & value)
  {
    set_dirichlet_boundary( domseg, storage, viennafvm::boundary_key(func_symbol.id()), value );
  }

  

//   template <typename CellType, typename NumericT, typename BoundaryAccessorType, typename BoundaryValueAccessorType>
//   void set_dirichlet_boundary(CellType const & c,
//                               std::vector<NumericT> const & value,
//                               BoundaryAccessorType boundary_accessor,
//                               BoundaryValueAccessorType boundary_value_accessor)
//   {
// //     typedef viennafvm::boundary_key      BoundaryKey;;
// 
//     //set flag:
//     boundary_accessor(c) = true;
// //     viennadata::access<BoundaryKey, bool >(BoundaryKey(id))(c) = true;
// 
//     //set data:
//     boundary_value_accessor(c) = value;
// //     viennadata::access<BoundaryKey, std::vector<NumericT> >(BoundaryKey(id))(c) = value;
//   }

  //
  // allow the use of function symbols (more intuitive)
  //

//   template <typename StorageType, typename CellType, typename InterfaceType>
//   void set_dirichlet_boundary(StorageType & storage,
//                               CellType const & c,
//                               numeric_type const & value,
//                               viennamath::rt_function_symbol<InterfaceType> const & func_symbol)
//   {
//     typename viennadata::result_of::accessor<StorageType, viennafvm::boundary_key, bool, CellType>::type boundary_accessor =
//         viennadata::accessor<viennafvm::boundary_key, bool, CellType>(storage, viennafvm::boundary_key(func_symbol.id()));
// 
//     typename viennadata::result_of::accessor<StorageType, viennafvm::boundary_key, numeric_type, CellType>::type boundary_value_accessor =
//         viennadata::accessor<viennafvm::boundary_key, numeric_type, CellType>(storage, viennafvm::boundary_key(func_symbol.id()));
// 
//     
//     set_dirichlet_boundary(c, value, boundary_accessor, boundary_value_accessor);
//   }
// 
//   template <typename StorageType, typename CellType, typename NumericT, typename InterfaceType>
//   void set_dirichlet_boundary(StorageType & storage,
//                               CellType const & c,
//                               std::vector<NumericT> const & value,
//                               viennamath::rt_function_symbol<InterfaceType> const & func_symbol)
//   {
//     typename viennadata::result_of::accessor<StorageType, viennafvm::boundary_key, long, CellType>::type boundary_accessor =
//         viennadata::accessor<viennafvm::boundary_key, long, CellType>(storage, viennafvm::boundary_key(func_symbol.id()));
// 
//     typename viennadata::result_of::accessor<StorageType, viennafvm::boundary_key, std::vector<NumericT>, CellType>::type boundary_value_accessor =
//         viennadata::accessor<viennafvm::boundary_key, std::vector<NumericT>, CellType>(storage, viennafvm::boundary_key(func_symbol.id()));
// 
//     set_dirichlet_boundary(c, value, boundary_accessor, boundary_value_accessor);
//   }





  
//   template <typename StorageType, typename CellType, typename InterfaceType>
//   bool is_dirichlet_boundary(CellType const & c, viennamath::rt_function_symbol<InterfaceType> const & func_symbol)
//   {
//     typedef viennafvm::boundary_key           BoundaryKey;
// 
//     return viennadata::access<BoundaryKey, bool >(BoundaryKey(func_symbol.id()))(c);
//   }
// 
//   template <typename StorageType, typename CellType, typename InterfaceType>
//   numeric_type get_dirichlet_boundary(CellType const & c, viennamath::rt_function_symbol<InterfaceType> const & func_symbol)
//   {
//     typedef viennafvm::boundary_key           BoundaryKey;
// 
//     return viennadata::access<BoundaryKey, numeric_type >(BoundaryKey(func_symbol.id()))(c);
//   }

}
#endif
