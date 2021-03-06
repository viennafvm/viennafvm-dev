#ifndef VIENNAFVM_CELL_QUAN_HPP
#define VIENNAFVM_CELL_QUAN_HPP

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
#include "viennadata/api.hpp"

#include "viennagrid/forwards.hpp"
#include "viennagrid/mesh/segmentation.hpp"
#include "viennagrid/algorithm/centroid.hpp"

/** @file  ncell_quantity.hpp
    @brief Defines ViennaMath extensions: Piecewise constants (constants on each cell) and flux evaluators on interfaces
*/

namespace viennafvm
{
  namespace detail
  {
    /** @brief The runtime interface for cell quantities.
    *
    * @param CellType    Type of the ViennaGrid cell
    * @param NumericT    Floating point type of the value to be returned (typically 'double')
    */
    template <typename CellType, typename NumericT = viennafvm::numeric_type>
    class ncell_quantity_interface
    {
      protected:
        typedef NumericT          numeric_type;

      public:
        virtual numeric_type eval(CellType const & cell, CellType const * outer_cell, numeric_type v) const = 0;
        virtual numeric_type eval(CellType const & cell, CellType const * outer_cell, std::vector<numeric_type> const & v) const = 0;

        virtual ncell_quantity_interface<CellType, NumericT> * clone() const = 0;
    };

    /** @brief Implementation of a function which is piecewise constant on each cell. Function values are retrieved from ViennaData.
    *
    * @param CellType    Type of the ViennaGrid cell
    * @param KeyType     The key type to be used with ViennaData
    * @param DataType    The data type to be used with ViennaData
    */
    template <typename CellType, typename QuantityT>
    class ncell_quantity_constant : public ncell_quantity_interface<CellType>
    {
        typedef ncell_quantity_constant<CellType, QuantityT>        self_type;
        typedef typename ncell_quantity_interface<CellType>::numeric_type    numeric_type;

      public:
        ncell_quantity_constant(QuantityT const & quan, bool is_gradient) : quan_(quan), is_gradient_(is_gradient) {}

        numeric_type eval(CellType const & cell, CellType const * outer_cell, numeric_type /*v*/) const
        {
          if (is_gradient_)
            return (quan_.get_value(*outer_cell) - quan_.get_value(cell)) / viennagrid::norm(viennagrid::centroid(*outer_cell) - viennagrid::centroid(cell));
          return quan_.get_value(cell);
        }

        numeric_type eval(CellType const & cell, CellType const * outer_cell, std::vector<numeric_type> const & /*v*/) const
        {
          if (is_gradient_)
            return (quan_.get_value(*outer_cell) - quan_.get_value(cell)) / viennagrid::norm(viennagrid::centroid(*outer_cell) - viennagrid::centroid(cell));
          return quan_.get_value(cell);
        }

        ncell_quantity_interface<CellType> * clone() const { return new self_type(quan_, is_gradient_); }

      private:
        QuantityT const & quan_;
        bool is_gradient_;
    };



    /** @brief A type erasure class which enables to store cell_quan_constants and cell_quan_exprs with different template arguments in a single array.
    *
    * @param CellType    Type of the ViennaGrid cell
    * @param NumericT    Floating point type of the value to be returned (typically 'double')
    */
    template <typename CellType, typename NumericT = viennafvm::numeric_type>
    class ncell_quantity_wrapper
    {
      public:
        template <typename T>
        ncell_quantity_wrapper(T const * t) : functor_(t) {}

        ncell_quantity_wrapper() {}

        ncell_quantity_wrapper & operator=(ncell_quantity_wrapper & other)
        {
          functor_ = other.functor_;
          return *this;
        }

        NumericT eval(CellType const & cell,
                      CellType const * outer_cell,
                      numeric_type v) const
        {
          return functor_->eval(cell, outer_cell, v);
        }

        NumericT eval(CellType const & cell,
                      CellType const * outer_cell,
                      std::vector<numeric_type> const & v) const
        {
          return functor_->eval(cell, outer_cell, v);
        }

        ncell_quantity_interface<CellType> * clone() const { return functor_->clone(); }

      private:
        std::auto_ptr< const ncell_quantity_interface<CellType> > functor_;
    };

  } //namespace detail


  /** @brief The main cell quantity class for using piecewise constant or piecewise expressions (in local coordinates) with ViennaMath.
   *
    * @param CellType       Type of the ViennaGrid cell
    * @param InterfaceType  The runtime interface class of ViennaMath.
   */
  template <typename CellType, typename InterfaceType>
  class ncell_quantity : public InterfaceType
  {
      typedef ncell_quantity<CellType, InterfaceType>     self_type;
    public:
      typedef typename InterfaceType::numeric_type            numeric_type;

      explicit ncell_quantity(CellType const * cell,
                              detail::ncell_quantity_wrapper<CellType, numeric_type> const & wrapper) : current_cell(cell), outer_cell(NULL), is_gradient_(false), accessor(wrapper.clone()) {}
      explicit ncell_quantity(CellType const * cell,
                              CellType const * cell_outer,
                              bool is_gradient,
                              detail::ncell_quantity_wrapper<CellType, numeric_type> const & wrapper) : current_cell(cell), outer_cell(cell_outer), is_gradient_(is_gradient), accessor(wrapper.clone()) {}

      //template <typename T>
      //explicit cell_quan(T const & t) : current_cell(NULL), accessor( new quan_accessor<CellType, T, numeric_type>() ) {}

      explicit ncell_quantity() : current_cell(NULL), outer_cell(NULL) {}

      //interface requirements:
      InterfaceType * clone() const { return new self_type(current_cell, outer_cell, is_gradient_, accessor); }
      numeric_type eval(std::vector<numeric_type> const & v) const
      {
        return accessor.eval(*current_cell, outer_cell, v);
      }
      numeric_type eval(numeric_type v) const
      {
        return accessor.eval(*current_cell, outer_cell, v);
      }

      std::string deep_str() const
      {
        std::stringstream ss;
        ss << "cell_quan<" << CellType::tag::dim << ">(" << name_ << ")";
        return ss.str();
      }
      numeric_type unwrap() const { throw "Cannot evaluate unknown_func!"; }

      InterfaceType * substitute(const InterfaceType * e,
                                 const InterfaceType * repl) const
      {
        if (deep_equal(e))
          return repl->clone();

        return clone();
      }

      InterfaceType * substitute(std::vector<const InterfaceType *> const &  e,
                                 std::vector<const InterfaceType *> const &  repl) const
      {
        //std::cout << "Comparing variable<" << id << "> with " << e->str() << ", result: ";
        for (std::size_t i=0; i<e.size(); ++i)
          if (deep_equal(e[i]))
            return repl[i]->clone();

        //std::cout << "FALSE" << std::endl;
        return clone();
      }

      bool deep_equal(const InterfaceType * other) const
      {
        //TODO: Include comparison of accessor
        return dynamic_cast< const self_type *>(other) != NULL;
      }

      bool shallow_equal(const InterfaceType * other) const
      {
        return dynamic_cast< const self_type *>(other) != NULL;
      }

      InterfaceType * diff(const InterfaceType * /*diff_var*/) const
      {
        return new viennamath::rt_constant<double, InterfaceType>(0);
      }


      //additional members:
      void update(CellType const & cell) const
      {
        current_cell = &cell;
      }

      void update(CellType const & cell, CellType const & cell_outer) const
      {
        if (is_gradient_)
        {
          current_cell = &cell;
          outer_cell   = &cell_outer;
        }
      }

      template <typename QuantityT>
      void wrap_constant(QuantityT const & quan, bool is_gradient)
      {
        detail::ncell_quantity_wrapper<CellType, numeric_type> temp(new detail::ncell_quantity_constant<CellType, QuantityT>(quan, is_gradient));
        accessor = temp;
        is_gradient_ = is_gradient;
        name_ = quan.get_name();
      }

      detail::ncell_quantity_wrapper<CellType, numeric_type> const & wrapper() const { return accessor; }

    private:
      mutable const CellType * current_cell;
      mutable const CellType * outer_cell;
      bool is_gradient_;
      std::string name_;
      detail::ncell_quantity_wrapper<CellType, numeric_type> accessor;
  };

  //TODO: Check whether cell_quan can be injected directly into existing ViennaMath overloads

  /** @brief Operator overload for the multiplication of a cell quantity with a ViennaMath variable */
  template <typename CellType, typename InterfaceType>
  viennamath::rt_expr<InterfaceType> operator*(viennamath::rt_variable<InterfaceType> const & lhs,
                                               ncell_quantity<CellType, InterfaceType> const & rhs)
  {
    return viennamath::rt_expr<InterfaceType>(new viennamath::rt_binary_expr<InterfaceType>(lhs.clone(),
                                                            new viennamath::op_binary<viennamath::op_mult<viennamath::default_numeric_type>, InterfaceType >(),
                                                            rhs.clone()));
  }


  /** @brief Operator overload for the multiplication of a cell quantity with a ViennaMath expression wrapper */
  template <typename CellType, typename InterfaceType>
  viennamath::rt_expr<InterfaceType> operator*(viennamath::rt_expr<InterfaceType> const & lhs,
                                               ncell_quantity<CellType, InterfaceType> const & rhs)
  {
    return viennamath::rt_expr<InterfaceType>(new viennamath::rt_binary_expr<InterfaceType>(lhs.get()->clone(),
                                                            new viennamath::op_binary<viennamath::op_mult<viennamath::default_numeric_type>, InterfaceType >(),
                                                            rhs.clone()));
  }

  /** @brief Operator overload for the multiplication of a cell quantity with a ViennaMath unary expression */
  template <typename CellType, typename InterfaceType>
  viennamath::rt_expr<InterfaceType> operator*(ncell_quantity<CellType, InterfaceType> const & lhs,
                                               viennamath::rt_unary_expr<InterfaceType> const & rhs
                               )
  {
    return viennamath::rt_expr<InterfaceType>(new viennamath::rt_binary_expr<InterfaceType>(lhs.clone(),
                                                            new viennamath::op_binary<viennamath::op_mult<viennamath::default_numeric_type>, InterfaceType >(),
                                                            rhs.clone()));
  }

  /** @brief Operator overload for the multiplication of a cell quantity with a ViennaMath binary expression */
  template <typename CellType, typename InterfaceType>
  viennamath::rt_expr<InterfaceType> operator*(ncell_quantity<CellType, InterfaceType> const & lhs,
                                               viennamath::rt_binary_expr<InterfaceType> const & rhs
                               )
  {
    return viennamath::rt_expr<InterfaceType>(new viennamath::rt_binary_expr<InterfaceType>(lhs.clone(),
                                                            new viennamath::op_binary<viennamath::op_mult<viennamath::default_numeric_type>, InterfaceType >(),
                                                            rhs.clone()));
  }

  /** @brief Operator overload for the multiplication of a cell quantity with a ViennaMath expression wrapper */
  template <typename CellType, typename InterfaceType>
  viennamath::rt_expr<InterfaceType> operator/(viennamath::rt_expr<InterfaceType> const & lhs,
                                               ncell_quantity<CellType, InterfaceType> const & rhs)
  {
    return viennamath::rt_expr<InterfaceType>(new viennamath::rt_binary_expr<InterfaceType>(lhs.get()->clone(),
                                                            new viennamath::op_binary<viennamath::op_div<viennamath::default_numeric_type>, InterfaceType >(),
                                                            rhs.clone()));
  }
}

// generic clone:
namespace viennamath
{

  /** @brief Any generic free functions for unifying interfaces are defined here. */
  namespace traits
  {
    template <typename InterfaceType, typename CellType>
    InterfaceType * clone(viennafvm::ncell_quantity<CellType, InterfaceType> const & c) { return c.clone(); }
  }
}


#endif
