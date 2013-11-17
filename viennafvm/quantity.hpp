#ifndef VIENNAFVM_QUANTITY_HPP
#define VIENNAFVM_QUANTITY_HPP

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

#include "viennafvm/forwards.h"

#include "viennamath/forwards.h"
#include "viennamath/manipulation/substitute.hpp"
#include "viennamath/expression.hpp"

#include "viennagrid/forwards.hpp"
#include "viennagrid/mesh/segmentation.hpp"

/** @file  quantity.hpp
    @brief Defines the basic quantity which may be entirely known or unknown (with suitable boundary conditions) over a mesh
*/

namespace viennafvm
{

  template<typename AssociatedT, typename ValueT = double>
  class quantity
  {
    public:
      typedef ValueT         value_type;
      typedef AssociatedT    associated_type;

      quantity() {}  // to fulfill default constructible concept!

      quantity(std::string const & quan_name,
               std::size_t num_values,
               value_type default_value = value_type())
        : name_(quan_name),
          values_                  (num_values, default_value),
          boundary_types_          (num_values, BOUNDARY_NONE),
          boundary_values_         (num_values, default_value),
          defined_but_unknown_mask_(num_values, false),
          unknowns_indices_        (num_values, -1)
      {}

      std::string get_name() const { return name_; }

      ValueT get_value(associated_type const & elem) const         { return values_.at(elem.id());         }
      void   set_value(associated_type const & elem, ValueT value) {        values_.at(elem.id()) = value; }

      // Dirichlet and Neumann
      ValueT get_boundary_value(associated_type const & elem) const         { return boundary_values_.at(elem.id());         }
      void   set_boundary_value(associated_type const & elem, ValueT value) {        boundary_values_.at(elem.id()) = value; }

      boundary_type_id get_boundary_type(associated_type const & elem) const                   { return boundary_types_.at(elem.id());         }
      void             set_boundary_type(associated_type const & elem, boundary_type_id value) {        boundary_types_.at(elem.id()) = value; }

      bool   get_unknown_mask(associated_type const & elem) const       { return defined_but_unknown_mask_.at(elem.id());         }
      void   set_unknown_mask(associated_type const & elem, bool value) {        defined_but_unknown_mask_.at(elem.id()) = value; }

      long   get_unknown_index(associated_type const & elem) const       { return unknowns_indices_.at(elem.id());         }
      void   set_unknown_index(associated_type const & elem, long value) {        unknowns_indices_.at(elem.id()) = value; }

      std::size_t get_unknown_num() const
      {
        std::size_t num = 0;
        for (std::size_t i=0; i<unknowns_indices_.size(); ++i)
        {
          if (unknowns_indices_[i] >= 0)
            ++num;
        }
        return num;
      }

      // possible design flaws:
      std::vector<ValueT> const & values() const { return values_; }

    private:
      std::string                    name_;
      std::vector<ValueT>            values_;
      std::vector<boundary_type_id>  boundary_types_;
      std::vector<ValueT>            boundary_values_;
      std::vector<bool>              defined_but_unknown_mask_;
      std::vector<long>              unknowns_indices_;
  };
}


#endif
