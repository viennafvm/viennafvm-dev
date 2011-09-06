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

#ifndef VIENNAFVM_LINEARPDESYSTEM_HPP
#define VIENNAFVM_LINEARPDESYSTEM_HPP

#include <vector>

// *** local includes:
#include "viennafvm/forwards.h"
#include "viennafvm/linear_pde_options.hpp"

// *** vienna includes:
#include "viennamath/forwards.h"
#include "viennadata/api.hpp"
#include "viennagrid/domain.hpp"

namespace viennafvm
{
  template <typename InterfaceType = viennamath::expression_interface<viennamath::default_numeric_type>, 
            typename MappingKeyType  = viennafvm::mapping_key, 
            typename BoundaryKeyType = viennafvm::boundary_key >
  struct linear_pde_system
  {
      typedef InterfaceType                                 interface_type;
      typedef MappingKeyType                                mapping_key_type;
      typedef BoundaryKeyType                               boundary_key_type;
      typedef viennamath::equation<InterfaceType>           equation_type;
      typedef viennamath::function_symbol<InterfaceType>    unknown_type;
      typedef std::vector< unknown_type >                   unknown_cont_type;
      typedef std::vector< std::string >                    key_cont_type;
      typedef viennafvm::linear_pde_options                 option_type;      
      
      void add_pde(equation_type      const & pde,
                   unknown_cont_type  const & unknowns, 
                   key_cont_type      const & keys, 
                   option_type        const & option)
      {
        pdes_.push_back(pde); 
        unknowns_.push_back(unknowns);
        keys_.push_back(keys);
        options_.push_back(option);
      }
      
      equation_type              pde(size_t index)       const { return pdes_[index]; }
      unknown_cont_type const &  unknown(size_t index)   const { return unknowns_[index]; }
      key_cont_type const &      keys(size_t index)      const { return keys_[index]; }      
      option_type                option(size_t index)    const { return options_[index]; }      
      
      size_t size() const { return pdes_.size(); }
      
    private:
      std::vector< equation_type >        pdes_;
      std::vector< unknown_cont_type >    unknowns_;
      std::vector< key_cont_type >        keys_;
      std::vector< option_type >          options_;
  };
  
  
   template <typename InterfaceType>
   linear_pde_system<InterfaceType> make_linear_pde_system(viennamath::equation<InterfaceType>            equ_1,
                                                           viennamath::function_symbol<InterfaceType>     unknown_1, 
                                                           std::string                                    key_1, 
                                                           viennafvm::linear_pde_options                  option_1)
   {
      typedef viennafvm::linear_pde_system<InterfaceType>   linear_pde_sys_type;

      linear_pde_sys_type ret;
      
      typename linear_pde_sys_type::unknown_cont_type unknown_vec_1(1);
      unknown_vec_1[0] = unknown_1;
   
      typename linear_pde_sys_type::key_cont_type     key_vec_1(1);
      key_vec_1[0] = key_1;

      ret.add_pde(equ_1, unknown_vec_1, key_vec_1, option_1);
      return ret;
   }
  
////  template <typename InterfaceType>
////  linear_pde_system<InterfaceType> make_linear_pde_system(viennamath::equation<InterfaceType> equ_1,
////                                                          std::vector<viennamath::function_symbol<InterfaceType> > unknowns_1)
////  {
////    linear_pde_system<InterfaceType> ret;
////    ret.add_pde(equ_1, unknowns_1);
////    return ret;
////  }
  
}
#endif

