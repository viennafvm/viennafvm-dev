#ifndef VIENNAFVM_FORWARDS_H
#define VIENNAFVM_FORWARDS_H

/* =======================================================================
   Copyright (c) 2011, Institute for Microelectronics, TU Wien
   http://www.iue.tuwien.ac.at
                             -----------------
           ViennaFVM - The Vienna Finite Volume Method Library
                             -----------------

   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at
               (add your name here)

   license:    To be discussed, see file LICENSE in the ViennaFVM base directory
======================================================================= */

#include <ostream>
#include "viennadata/api.hpp"

#define NUM_PI 3.1415926535897932384626433832795

namespace viennafvm
{

  /** @brief The default floating point type to be used in ViennaFVM.
   *
   *  Feel free to change this typedef to a high-precision type if required.
   *  Keep in mind that only float and double types can be used for GPU acceleration.
   */
  typedef double             numeric_type;

  enum
  {
    cell_volume = 0,
    cell_boundary = 1
  };


  // define keys and configure ViennaData to use a type-based dispatch:
  class boundary_key
  {
    public:
      boundary_key(long id) : id_(id) {}

      bool operator<(boundary_key const & other) const { return id_ < other.id_; }

    private:
      long id_;
  };

  class mapping_key
  {
    public:
      mapping_key(long id) : id_(id) {}

      bool operator<(mapping_key const & other) const { return id_ < other.id_; }

    private:
      long id_;
  };

  class current_iterate_key
  {
    public:
      current_iterate_key(long id) : id_(id) {}

      bool operator<(current_iterate_key const & other) const { return id_ < other.id_; }

    private:
      long id_;
  };

  //box-integration related:
  class facet_distance_key {};
  class facet_area_key {};

}



namespace viennadata
{
  namespace config
  {
    //box-integration related:
    template <>
    struct key_dispatch<viennafvm::facet_distance_key>
    {
      typedef type_key_dispatch_tag    tag;
    };

    template <>
    struct key_dispatch<viennafvm::facet_area_key>
    {
      typedef type_key_dispatch_tag    tag;
    };
  }
}

#endif
