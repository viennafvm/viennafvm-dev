#ifndef VIENNAFVM_PROBLEM_DESCRIPTION_HPP
#define VIENNAFVM_PROBLEM_DESCRIPTION_HPP

/* =======================================================================
   Copyright (c) 2011, Institute for Microelectronics, TU Wien
   http://www.iue.tuwien.ac.at
                             -----------------
           ViennaFVM - The Vienna Finite Volume Method Library
                             -----------------

   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at
               Josef Weinbub                   weinbub@iue.tuwien.ac.at
               (add your name here)

   license:    see file LICENSE in the ViennaFVM base directory
======================================================================= */

#include <stdexcept>

#include "viennagrid/mesh/mesh.hpp"

#ifdef VIENNAFVM_VERBOSE
#include "viennafvm/timer.hpp"
#endif
#include "viennafvm/forwards.h"
#include "viennafvm/quantity.hpp"
#include "viennafvm/linear_assembler.hpp"
#include "viennafvm/linear_solvers/viennacl.hpp"

namespace viennafvm
{

  /** @brief Exception for the case that an invalid quantity is accessed */
  class quantity_not_found_exception : public std::runtime_error {
  public:
    quantity_not_found_exception(std::string const & str) : std::runtime_error(str) {}
  };


  template<typename MeshT>
  class problem_description
  {
      typedef typename viennagrid::result_of::cell_tag<MeshT>::type                CellTag;
      typedef typename viennagrid::result_of::element<MeshT, CellTag>::type        CellType;

    public:
      typedef viennafvm::numeric_type numeric_type;

      typedef quantity<CellType, numeric_type>   QuantityType;
      typedef QuantityType                       quantity_type;

      typedef std::vector<quantity_type>   quantity_container_type;
      typedef quantity_container_type      QuantityContainerType;

      typedef MeshT        MeshType;
      typedef MeshType     mesh_type;

      problem_description() : mesh_(NULL) {}  // set mesh via operator= or copy-CTOR

      problem_description(MeshT const & mesh) : mesh_(&mesh) {}

      quantity_container_type const & quantities() const { return quantities_; }
      quantity_container_type       & quantities()       { return quantities_; }

      QuantityType & add_quantity(std::string name, numeric_type default_value = numeric_type())
      {
        quantities_.push_back(QuantityType(name, viennagrid::cells(*mesh_).size(), default_value));
        return quantities_.back();
      }

      QuantityType & get_quantity(std::string name)
      {
        for (std::size_t i=0; i<quantities_.size(); ++i)
          if (quantities_[i].get_name() == name)
            return quantities_[i];

        throw quantity_not_found_exception(name);
      }

      QuantityType const & get_quantity(std::string name) const
      {
        for (std::size_t i=0; i<quantities_.size(); ++i)
          if (quantities_[i].get_name() == name)
            return quantities_[i];

        throw quantity_not_found_exception(name);
      }

      MeshType const & mesh() const { return *mesh_; }

    private:
      MeshType const *         mesh_;
      quantity_container_type  quantities_;
  };

}

#endif // VIENNAFVM_PDE_SOLVER_HPP

