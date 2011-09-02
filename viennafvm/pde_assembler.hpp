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
#include "viennafvm/pde_stack.hpp"
#include "viennafvm/acc.hpp"

#define VIENNAFVMDEBUG

namespace viennafvm
{

   struct pde_assembler
   {
      template <typename PDESetT, typename SegmentT, typename MatrixT, typename VectorT>  
      void operator()(PDESetT  & pdeset,
                      SegmentT & segment,
                      MatrixT  & matrix,
                      VectorT  & rhs
                     ) 
      {
         typedef typename SegmentT::config_type                config_type;      
         typedef viennamath::equation<>                        equ_type;
      
         static const int DIM = config_type::dimension_tag::value;
      
         // compute the voronoi information and store it on the segment with 
         // the provided keys
         //
         viennagrid::write_voronoi_info<viennafvm::edge_len_key,
                                        viennafvm::edge_interface_area_key,
                                        viennafvm::box_volume_key>(segment);
                                        
         // traverse the pdes which should be assembled
         //
         for(typename viennafvm::pde_set::iterator pdeit = pdeset.begin(); 
             pdeit != pdeset.end(); pdeit++)
         {
         #ifdef VIENNAFVMDEBUG
            std::cout << "strong form: " << std::endl;
            std::cout << viennamos::acc<viennafvm::tag::pde>(*pdeit) << std::endl;
         #endif

            equ_type weak_form = viennafvm::make_weak_form( viennamos::acc<viennafvm::tag::pde>(*pdeit) );
            
         #ifdef VIENNAFVMDEBUG
            std::cout << "weak form: " << std::endl;
            std::cout << weak_form << std::endl;
         #endif

/*
   1. extract all dOmegas and Omegas, and seperate them 
   2. process all dOmegas and Omegas in singular 'loops':
      dOmega integrals: edge_on_vertex/vertex_on_edge 
      Omega integrals:  work on main vertex
*/

         }

      }
   };
   
   template <typename PDESetT, typename SegmentT, typename MatrixT, typename VectorT>  
   void assemble_pde( PDESetT  & pdeset,
                      SegmentT & segment,
                      MatrixT  & matrix,
                      VectorT  & rhs
                     ) 
   {
      viennafvm::pde_assembler()(pdeset, segment, matrix, rhs);
   }
}

#endif 







