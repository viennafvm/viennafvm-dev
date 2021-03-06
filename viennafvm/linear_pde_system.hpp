#ifndef VIENNAFVM_LINEARPDESYSTEM_HPP
#define VIENNAFVM_LINEARPDESYSTEM_HPP

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

#include <vector>

// *** local includes:
#include "viennafvm/forwards.h"

// *** vienna includes:
#include "viennamath/forwards.h"
#include "viennagrid/mesh/mesh.hpp"

namespace viennafvm
{


  class linear_pde_options
  {
    public:
      explicit linear_pde_options(long id = 0) : data_id_(id), check_mapping_(false), geometric_update_(false), damping_term_(viennamath::rt_constant<numeric_type>(0)) {}

      long data_id() const { return data_id_; }
      void data_id(long new_id) { data_id_ = new_id; }

      bool check_existing_mapping() const { return check_mapping_; }
      void check_existing_mapping(bool b) { check_mapping_ = b; }

      bool geometric_update() const { return geometric_update_; }
      void geometric_update(bool b) { geometric_update_ = b; }

      viennamath::expr damping_term() const { return damping_term_; }
      void damping_term(viennamath::expr const & e) { damping_term_ = e; }

    private:
      long data_id_;
      bool check_mapping_;
      bool geometric_update_;
      viennamath::expr damping_term_;
  };

  inline linear_pde_options make_linear_pde_options(long data_id, bool existing_mapping = false)
  {
    linear_pde_options options;
    options.data_id(data_id);
    options.check_existing_mapping(existing_mapping);
    return options;
  }


  template <typename InterfaceType = viennamath::default_interface_type>
  struct linear_pde_system
  {
      typedef InterfaceType                                 interface_type;
      typedef viennamath::equation                          equation_type;
      typedef viennamath::function_symbol                   unknown_type;
      typedef std::vector< unknown_type >                   unknown_cont_type;
      typedef viennafvm::linear_pde_options                 option_type;

      linear_pde_system() : is_linear_(true) {}

      void add_pde(equation_type      const & pde,
                   unknown_cont_type  const & unknowns,
                   option_type        const & option)
      {
        pdes_.push_back(pde);
        unknowns_.push_back(unknowns);
        options_.push_back(option);
      }

      void add_pde(equation_type      const & pde,
                   unknown_cont_type  const & unknowns)
      {
        std::size_t current_id = size();

        pdes_.push_back(pde);
        unknowns_.push_back(unknowns);
        options_.push_back(option_type(current_id));
      }

      void add_pde(equation_type      const & pde,
                   unknown_type       const & unknown)
      {
        std::size_t current_id = size();

        pdes_.push_back(pde);

        unknown_cont_type new_unknown(1); new_unknown[0] = unknown;
        unknowns_.push_back(new_unknown);

        options_.push_back(option_type(current_id));
      }

      equation_type const &      pde(size_t index)       const { return pdes_[index]; }
      unknown_cont_type const &  unknown(size_t index)   const { return unknowns_[index]; }
      option_type const &        option(size_t index)    const { return options_[index]; }
      option_type       &        option(size_t index)          { return options_[index]; }

      bool is_linear() const { return is_linear_; }
      void is_linear(bool b) { is_linear_ = b; }

      size_t size() const { return pdes_.size(); }

    private:
      std::vector< equation_type >        pdes_;
      std::vector< unknown_cont_type >    unknowns_;
      std::vector< option_type >          options_;
      bool is_linear_;
  };


   template <typename InterfaceType>
   linear_pde_system<InterfaceType> make_linear_pde_system(viennamath::equation             equ_1,
                                                           viennamath::function_symbol      unknown_1,
                                                           viennafvm::linear_pde_options    option_1)
   {
      typedef viennafvm::linear_pde_system<InterfaceType>   linear_pde_sys_type;

      linear_pde_sys_type ret;

      typename linear_pde_sys_type::unknown_cont_type unknown_vec_1(1);
      unknown_vec_1[0] = unknown_1;

      ret.add_pde(equ_1, unknown_vec_1, option_1);
      return ret;
   }

   template <typename InterfaceType>
   linear_pde_system<InterfaceType> make_linear_pde_system(viennamath::rt_equation<InterfaceType>         equ_1,
                                                           viennamath::rt_function_symbol<InterfaceType>  unknown_1)
   {
      typedef viennafvm::linear_pde_system<InterfaceType>   linear_pde_sys_type;

      linear_pde_sys_type ret;

      linear_pde_options option_1;

      typename linear_pde_sys_type::unknown_cont_type unknown_vec_1(1);
      unknown_vec_1[0] = unknown_1;

      ret.add_pde(equ_1, unknown_vec_1, option_1);
      return ret;
   }

   template <typename InterfaceType>
   linear_pde_system<InterfaceType> make_linear_pde_system(viennamath::rt_equation<InterfaceType>         equ_1,
                                                           viennamath::rt_function_symbol<InterfaceType>  unknown_1,
                                                           viennafvm::linear_pde_options    option_1
                                                          )
   {
      typedef viennafvm::linear_pde_system<InterfaceType>   linear_pde_sys_type;

      linear_pde_sys_type ret;

      typename linear_pde_sys_type::unknown_cont_type unknown_vec_1(1);
      unknown_vec_1[0] = unknown_1;

      ret.add_pde(equ_1, unknown_vec_1, option_1);
      return ret;
   }


}
#endif

