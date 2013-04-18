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

#ifndef VIENNAFVM_PDE_ASSEMBLER_HPP
#define VIENNAFVM_PDE_ASSEMBLER_HPP

// *** local includes
//
#include "viennafvm/integral_form.hpp"
// #include "viennafvm/pde_stack.hpp" obsolete
// #include "viennafvm/acc.hpp" obsolete
#include "viennafvm/extract_integrals.hpp" 
#include "viennafvm/rhs_zero.hpp"
#include "viennafvm/linear_pde_system.hpp"
#include "viennafvm/mapping.hpp"

#include "viennagrid/forwards.h"
#include "viennagrid/algorithm/voronoi.hpp"
#include "viennagrid/algorithm/volume.hpp"
#include "viennagrid/algorithm/centroid.hpp"

#include "viennamath/manipulation/eval.hpp"

//#define VIENNAFVMDEBUG

namespace viennafvm
{
  namespace detail
  {
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
  }

   struct linear_assembler
   {
      template <typename LinPdeSysT, 
                typename SegmentT, 
                typename MatrixT,
                typename VectorT>  
      void operator()(LinPdeSysT pde_system,
                      SegmentT   & segment,
                      MatrixT    & system_matrix, 
                      VectorT    & load_vector)
      {
        typedef typename SegmentT::config_type                config_type;
        typedef viennamath::equation                          equ_type;
        typedef viennamath::expr                              expr_type;
        typedef typename expr_type::numeric_type              numeric_type;
         
        typedef typename SegmentT::config_type              Config;
        typedef typename Config::cell_tag                   CellTag;
        typedef typename viennagrid::result_of::ncell<Config, CellTag::dim-1>::type             FacetType;
        typedef typename viennagrid::result_of::ncell<Config, CellTag::dim  >::type             CellType;
      
        typedef typename viennagrid::result_of::ncell_range<SegmentT, CellTag::dim>::type        CellContainer;
        typedef typename viennagrid::result_of::iterator<CellContainer>::type                    CellIterator;

        typedef typename viennagrid::result_of::ncell_range<CellType, CellTag::dim-1>::type      FacetOnCellContainer;
        typedef typename viennagrid::result_of::iterator<FacetOnCellContainer>::type             FacetOnCellIterator;


         typedef typename LinPdeSysT::mapping_key_type   MappingKeyType;
         MappingKeyType  map_key(pde_system.option(0).data_id());      
      
         typedef typename LinPdeSysT::boundary_key_type  BoundaryKeyType;
         BoundaryKeyType bnd_key(pde_system.option(0).data_id());


         std::size_t map_index = viennafvm::create_mapping(pde_system, segment);

         std::cout << "map index: " << map_index << std::endl;
         system_matrix.clear();
         system_matrix.resize(map_index, map_index, false);
         load_vector.clear();
         load_vector.resize(map_index);

      #ifdef VIENNAFVM_DEBUG
         std::cout << "strong form: " << std::endl;
         std::cout << pde_system.pde(0) << std::endl;
      #endif

         equ_type integral_form = viennafvm::make_integral_form( pde_system.pde(0) );
         
      #ifdef VIENNAFVM_DEBUG
         std::cout << "integral form: " << std::endl;
         std::cout << integral_form << std::endl;
      #endif 

         /*
        equ_type weak_form_rhs_zero = viennafvm::make_rhs_zero( weak_form );

      #ifdef VIENNAFVMDEBUG
        std::cout << "rhs zero: " << std::endl;
        std::cout << weak_form_rhs_zero << std::endl;
      #endif
      */

        //
        // Preprocess symbolic representation:
        //

        expr_type  partial_omega_integrand(dynamic_cast<viennamath::unary_expr const * >(integral_form.lhs().get())->lhs()->clone()); //TODO: Unhack!
        expr_type          omega_integrand(dynamic_cast<viennamath::unary_expr const * >(integral_form.rhs().get())->lhs()->clone()); //TODO: Unhack!

\
        //
        // Actual assembly:
        //
        CellContainer cells = viennagrid::ncells(segment);
        for (CellIterator cit = cells.begin(); cit != cells.end(); ++cit)
        {
          long row_index = viennadata::access<MappingKeyType, long>(map_key)(*cit);

          if (row_index < 0)
            continue;

          //
          // Boundary integral terms:
          //
          FacetOnCellContainer facets_on_cell = viennagrid::ncells(*cit);
          for (FacetOnCellIterator focit  = facets_on_cell.begin();
                                   focit != facets_on_cell.end();
                                 ++focit)
          {
            CellType const * other_cell = detail::other_cell_of_facet(*focit, *cit, segment);

            if (other_cell)
            {
              long col_index = viennadata::access<MappingKeyType, long>(map_key)(*other_cell);

              double flux = 1.0 / viennagrid::norm(viennagrid::centroid(*cit) - viennagrid::centroid(*other_cell));
              double area = viennagrid::volume(*focit);

              system_matrix(row_index, row_index) -= flux * area;

              if (col_index < 0)
              {
                double boundary_value = viennadata::access<BoundaryKeyType, double>(bnd_key)(*other_cell);

                load_vector(row_index) -= flux * area * boundary_value;
              }
              system_matrix(row_index, col_index) += flux * area;

            }
          }

          //
          // Volume terms
          //
          double cell_volume      = viennagrid::volume(*cit);
          std::vector<double> p(3);
          load_vector(row_index) += viennamath::eval(omega_integrand, p) * cell_volume;
          //std::cout << "Writing " << viennamath::eval(omega_integrand, p) << " * " << cell_volume << " to rhs at " << row_index << std::endl;
        }

         
      } // functor
   };
   
   template <typename InterfaceType, typename SegmentT, typename MatrixT, typename VectorT>  
   void assemble_pde( viennafvm::linear_pde_system<InterfaceType>  & pde_system,
                      SegmentT & segment,
                      MatrixT  & matrix,
                      VectorT  & rhs
                     ) 
   {
      viennafvm::linear_assembler()(pde_system, segment, matrix, rhs);
   }
}

#endif 







