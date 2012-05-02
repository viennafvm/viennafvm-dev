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

#ifndef VIENNAFVM_LINEARPDEOPTIONS_HPP
#define VIENNAFVM_LINEARPDEOPTIONS_HPP

// *** system includes
//
#include <vector>

// *** local includes:
//
//#include "viennafem/forwards.h"

// *** vienna includes:
//
#include "viennadata/api.hpp"
#include "viennagrid/domain.hpp"

namespace viennafvm
{
  
  class linear_pde_options
  {
    public: 
      explicit linear_pde_options() : data_id_(0), check_mapping_(false) {}
      
      long data_id() const { return data_id_; }
      void data_id(long new_id) { data_id_ = new_id; }
      
      bool check_existing_mapping() const { return check_mapping_; }
      void check_existing_mapping(bool b) { check_mapping_ = b; }
      
    private:
      long data_id_;
      bool check_mapping_;
  };
  
  linear_pde_options make_linear_pde_options(long data_id, bool existing_mapping = false)
  {
    linear_pde_options options;
    options.data_id(data_id);
    options.check_existing_mapping(existing_mapping);
    return options;
  }
  
}
#endif

