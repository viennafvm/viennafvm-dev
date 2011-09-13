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
#include "viennafvm/weak_form.hpp"
// #include "viennafvm/pde_stack.hpp" obsolete
// #include "viennafvm/acc.hpp" obsolete
#include "viennafvm/extract_integrals.hpp" 
#include "viennafvm/rhs_zero.hpp"
#include "viennafvm/linear_pde_system.hpp"
#include "viennafvm/mapping.hpp"

#define VIENNAFVMDEBUG

namespace viennafvm
{

   struct pde_assembler
   {
      template <typename LinPdeSysT, 
                typename SegmentT, 
                typename MatrixT, 
                typename VectorT>  
      void operator()(LinPdeSysT & pde_system,
                      SegmentT   & segment,
                      MatrixT    & matrix,
                      VectorT    & rhs) 
      {
         typedef typename SegmentT::config_type                config_type;      
         typedef viennamath::equation<>                        equ_type;
         typedef viennamath::expr<>                            expr_type;
         typedef typename expr_type::numeric_type              numeric_type;
         
         typedef viennamath::op_symbolic_integration<numeric_type, viennafvm::PartialOmega>     partial_omega_tag;
         typedef viennamath::op_symbolic_integration<numeric_type, viennafvm::Omega>            omega_tag;         
      
         typedef typename SegmentT::config_type              Config;
         typedef typename viennagrid::result_of::ncell<Config, 0>::type                         VertexType;
         typedef typename viennagrid::result_of::ncell<Config, 1>::type                         EdgeType;      
      
         typedef typename viennagrid::result_of::ncell_range<SegmentT, 0>::type        VertexContainer;
         typedef typename viennagrid::result_of::iterator<VertexContainer>::type       VertexIterator;

         typedef typename viennagrid::result_of::ncell_range<VertexType, 1>::type      EdgeOnVertexContainer;
         typedef typename viennagrid::result_of::iterator<EdgeOnVertexContainer>::type EdgeOnVertexIterator;

         typedef typename viennagrid::result_of::ncell_range<EdgeType, 0>::type        VertexOnEdgeContainer;
         typedef typename viennagrid::result_of::iterator<VertexOnEdgeContainer>::type VertexOnEdgeIterator;      
      
         typedef typename LinPdeSysT::mapping_key_type   MappingKeyType;      
         MappingKeyType  map_key(pde_system.option(0).data_id());      
      
         //static const int DIM = config_type::dimension_tag::value;
      
         // compute the voronoi information and store it on the segment with 
         // the provided keys
         //
         viennagrid::apply_voronoi(segment,  viennagrid::voronoi_interface_area_key(), 
                                             viennagrid::voronoi_box_volume_key());
                                        
         size_t map_index = viennafvm::create_mapping(pde_system, segment);

         std::cout << "map index: " << map_index << std::endl;
     
      #ifdef VIENNAFVMDEBUG
         std::cout << "strong form: " << std::endl;
         std::cout << pde_system.pde(0) << std::endl;
      #endif

         equ_type weak_form = viennafvm::make_weak_form( pde_system.pde(0) );
         
      #ifdef VIENNAFVMDEBUG
         std::cout << "weak form: " << std::endl;
         std::cout << weak_form << std::endl;
      #endif 

//            equ_type weak_form_rhs_zero = viennafvm::make_rhs_zero( weak_form );

//         #ifdef VIENNAFVMDEBUG
//            std::cout << "rhs zero: " << std::endl;
//            std::cout << weak_form_rhs_zero << std::endl;
//         #endif

//            viennafvm::extract_integrals::result_type  partial_omega_terms = viennafvm::extract_integrals::eval(weak_form_rhs_zero.lhs(), partial_omega_tag());
//            viennafvm::extract_integrals::result_type  omega_terms = viennafvm::extract_integrals::eval(weak_form_rhs_zero, omega_tag());

         


         long row_index, col_index;
         double matrix_entry;

         // common part
         //
         VertexContainer vertices = viennagrid::ncells<0>(segment);
         for (VertexIterator vit = vertices.begin(); vit != vertices.end(); ++vit)
         {
            row_index = viennadata::access<MappingKeyType, long>(map_key)(*vit);

            // if vertex is an interior vertex ..
            //   note: we don't assemble the boundary vertices ..
            //
            if (row_index >= 0)
            {
               // dOmega part
               //  
               EdgeOnVertexContainer edges = viennagrid::ncells<1>(*vit, segment);
               for (EdgeOnVertexIterator eovit = edges.begin(); eovit != edges.end(); ++eovit)
               {
                  VertexOnEdgeContainer vertices_on_edge = viennagrid::ncells<0>(*eovit);
                  VertexOnEdgeIterator voeit = vertices_on_edge.begin();
                  
                  if ( &(*voeit) == &(*vit))  //one of the two vertices of the edge is different from *vit
                     ++voeit;                  
                  
                  col_index        = viennadata::access<MappingKeyType, long>(map_key)(*voeit);                  
                  
                  //std::cout << "row: " << row_index << " - col: " << col_index << std::endl; 
                  
                  
                  matrix_entry =  
                     viennadata::access<viennafvm::edge_interface_area_key, double>()(*eovit) 
                     / 
                     viennadata::access<viennafvm::edge_len_key, double>()(*eovit);
                  
                  // if the neighbour vertex is an interior vertex ...
                  // this part works on the off-diagonal part of the matrix
                  //
                  if(col_index >= 0)
                  {
                     matrix(row_index, col_index) += matrix_entry;
                  }
                  // if neighbour vertex is a boundary vertex ... 
                  //
                  else
                  {
                     // multiply with dirichlet value and subtract from rhs
                     //
                     rhs(row_index) -= matrix_entry * viennadata::access<viennamos::tag::potential, double>()(*voeit);
                  }

                  //std::cout << "row: " << row_index << " - col: " << col_index << std::endl;

                  // accumulate main diagonal entries
                  //
                  matrix(row_index, row_index) -= matrix_entry;
                  
               } // end edge-on-vertex traversal
               
               //rhs(row_index) += viennadata::access<box_volume_key, double>()(*vit);

            } // end interior vertex

         } // end vertex traversal
         
      } // functor
   };
   
   template <typename InterfaceType, typename SegmentT, typename MatrixT, typename VectorT>  
   void assemble_pde( viennafvm::linear_pde_system<InterfaceType>  & pde_system,
                      SegmentT & segment,
                      MatrixT  & matrix,
                      VectorT  & rhs
                     ) 
   {
      viennafvm::pde_assembler()(pde_system, segment, matrix, rhs);
   }
}

#endif 







