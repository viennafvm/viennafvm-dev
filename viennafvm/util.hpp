#ifndef VIENNAFVM_UTIL_HPP
#define VIENNAFVM_UTIL_HPP

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

#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <cmath>

#include "viennadata/api.hpp"
#include "viennafvm/forwards.h"
#include "viennagrid/algorithm/voronoi.hpp"



/** @file viennafvm/util.hpp
    @brief Miscellaneous utilities
*/

namespace viennafvm
{

  namespace util
  {

    ///////////////////// outer unit vectors on facet ///////////////////
    template <typename FacetType, typename CellType>
    typename viennagrid::result_of::point<CellType>::type
    unit_outer_normal(FacetType const & facet, CellType const & cell);

    // 1d
    template <typename FacetType, typename ConfigType>
    typename viennagrid::result_of::point<FacetType>::type
    unit_outer_normal(FacetType const & facet, viennagrid::element_t<ConfigType, viennagrid::line_tag> const & cell)
    {
      typedef typename viennagrid::result_of::point<FacetType>::type   PointType;
      typedef viennagrid::element_t<ConfigType, viennagrid::line_tag> CellType;
      typedef typename ConfigType::cell_tag                           CellTag;

      typedef typename viennagrid::result_of::const_ncell_range<CellType, CellTag::dim-1> ::type FacetOnCellContainer;
      typedef typename viennagrid::result_of::iterator<FacetOnCellContainer>::type               FacetOnCellIterator;

      FacetType const * other_facet = &(viennagrid::ncells<CellTag::dim-1>(cell)[0]);
      if (other_facet->point()[0] == facet.point()[0])
        other_facet = &(viennagrid::ncells<CellTag::dim-1>(cell)[1]);

      PointType n;

      n[0] = (facet.point()[0] < other_facet->point()[0]) ? -1: 1;  //left facet has normal -1, right facet has normal +1

      return n;
    }

    // 2d
    namespace detail
    {
      template <typename FacetType, typename CellType>
      typename viennagrid::result_of::point<FacetType>::type
      unit_outer_normal_2d(FacetType const & facet, CellType const & cell)
      {
        typedef typename CellType::config_type                             ConfigType;
        typedef typename viennagrid::result_of::point<CellType>::type      PointType;

        typedef typename viennagrid::result_of::ncell<ConfigType, 0>::type VertexType;

        //
        // Find a vertex which is not part of the face
        //
        VertexType const * non_facet_vertex = NULL;

        for (std::size_t i=0; i<viennagrid::ncells<0>(cell).size(); ++i)
        {
          non_facet_vertex = &(viennagrid::ncells<0>(cell)[i]);
          bool is_part_of_facet = false;

          for (std::size_t j=0; j<viennagrid::ncells<0>(facet).size(); ++j)
          {
            if (non_facet_vertex == &(viennagrid::ncells<0>(facet)[j]))
            {
              is_part_of_facet = true;
              break;
            }
          }
          if (!is_part_of_facet)
            break;
        }

        //
        // Now compute direction vector
        //
        PointType edge_vec  = viennagrid::ncells<0>(facet)[1].point() - viennagrid::ncells<0>(facet)[0].point();
        edge_vec /= viennagrid::norm(edge_vec);

        PointType other_vec = non_facet_vertex->point() - viennagrid::ncells<0>(facet)[0].point();
        other_vec -= viennagrid::inner_prod(edge_vec, other_vec) * edge_vec; //orthogonalize (one step of Gram-Schmidt)
        other_vec /= -1.0 * viennagrid::norm(other_vec); //make it unit length and flip direction

        return other_vec;
      }
    }

    // interface for triangles
    template <typename FacetType, typename ConfigType>
    typename viennagrid::result_of::point<FacetType>::type
    unit_outer_normal(FacetType const & facet, viennagrid::element_t<ConfigType, viennagrid::triangle_tag> const & cell)
    {
      return detail::unit_outer_normal_2d(facet, cell);
    }

    // interface for quadrilaterals
    template <typename FacetType, typename ConfigType>
    typename viennagrid::result_of::point<FacetType>::type
    unit_outer_normal(FacetType const & facet, viennagrid::element_t<ConfigType, viennagrid::quadrilateral_tag> const & cell)
    {
      return detail::unit_outer_normal_2d(facet, cell);
    }


    // 3d
    namespace detail
    {
      template <typename FacetType, typename CellType>
      typename viennagrid::result_of::point<FacetType>::type
      unit_outer_normal_3d(FacetType const & facet, CellType const & cell)
      {
        typedef typename CellType::config_type                             ConfigType;
        typedef typename viennagrid::result_of::point<CellType>::type      PointType;

        typedef typename viennagrid::result_of::ncell<ConfigType, 0>::type VertexType;

        //
        // Find a vertex which is not part of the face
        //
        VertexType const * non_facet_vertex = NULL;

        for (std::size_t i=0; i<viennagrid::ncells<0>(cell).size(); ++i)
        {
          non_facet_vertex = &(viennagrid::ncells<0>(cell)[i]);
          bool is_part_of_facet = false;

          for (std::size_t j=0; j<viennagrid::ncells<0>(facet).size(); ++j)
          {
            if (non_facet_vertex == &(viennagrid::ncells<0>(facet)[j]))
            {
              is_part_of_facet = true;
              break;
            }
          }
          if (!is_part_of_facet)
            break;
        }

        //
        // Now compute direction vector via cross-prod
        //
        PointType edge_vec1  = viennagrid::ncells<0>(facet)[1].point() - viennagrid::ncells<0>(facet)[0].point();
        PointType edge_vec2  = viennagrid::ncells<0>(facet)[2].point() - viennagrid::ncells<0>(facet)[0].point();
        PointType normal_vec = viennagrid::cross_prod(edge_vec1, edge_vec2);
        normal_vec /= viennagrid::norm(normal_vec);

        // ensure normal vector is pointing outwards:
        PointType other_vec = non_facet_vertex->point() - viennagrid::ncells<0>(facet)[0].point();
        if (viennagrid::inner_prod(normal_vec, other_vec) > 0)
          normal_vec *= -1.0;

        return normal_vec;
      }
    }

    // interface for tetrahedra
    template <typename FacetType, typename ConfigType>
    typename viennagrid::result_of::point<FacetType>::type
    unit_outer_normal(FacetType const & facet, viennagrid::element_t<ConfigType, viennagrid::tetrahedron_tag> const & cell)
    {
      return detail::unit_outer_normal_3d(facet, cell);
    }

    // interface for hexahedra
    template <typename FacetType, typename ConfigType>
    typename viennagrid::result_of::point<FacetType>::type
    unit_outer_normal(FacetType const & facet, viennagrid::element_t<ConfigType, viennagrid::hexahedron_tag> const & cell)
    {
      return detail::unit_outer_normal_3d(facet, cell);
    }





    //
    // Other cell of facet
    //


    template <typename FacetType, typename CellType, typename DomainType>
    CellType const * other_cell_of_facet(FacetType const & facet, CellType const & cell, DomainType const & domain)
    {
      typedef typename CellType::tag      CellTag;

      typedef typename viennagrid::result_of::const_ncell_range<FacetType, CellTag::dim>::type  CellOnFacetRange;
      typedef typename viennagrid::result_of::iterator<CellOnFacetRange>::type                  CellOnFacetIterator;

      CellOnFacetRange    cells = viennagrid::ncells<CellTag::dim>(facet, domain);
      CellOnFacetIterator cofit = cells.begin();

      if (&(*cofit) == &cell) // we know the first cell pointed to by the iterator already, so we pick the 'other'
        ++cofit;

      if (cofit != cells.end())
        return &(*cofit);

      return NULL;  // facet is part of one cell only, so there is no 'other' cell
    }

  } //namespace util

} //namespace viennafvm
#endif
