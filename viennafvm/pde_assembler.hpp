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


namespace viennafvm
{

   class pde_assembler
   {
      public:
      
      template <typename SystemType, typename DomainType, typename MatrixType, typename VectorType>  //template for operator()
      void operator()(SystemType pde_system,
                      DomainType & domain,
                      MatrixType & system_matrix,
                      VectorType & load_vector
                     ) const
      {
      }
   };
   
}

#endif 



