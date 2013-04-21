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

#include "viennagrid/forwards.h"
#include "viennagrid/segment.hpp"

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
        virtual numeric_type eval(CellType const & cell, numeric_type v) const = 0;
        virtual numeric_type eval(CellType const & cell, std::vector<numeric_type> const & v) const = 0;

        virtual ncell_quantity_interface<CellType, NumericT> * clone() const = 0;
    };

    /** @brief Implementation of a function which is piecewise constant on each cell. Function values are retrieved from ViennaData.
    *
    * @param CellType    Type of the ViennaGrid cell
    * @param KeyType     The key type to be used with ViennaData
    * @param DataType    The data type to be used with ViennaData
    */
    template <typename CellType, typename KeyType, typename DataType>
    class ncell_quantity_constant : public ncell_quantity_interface<CellType>
    {
        typedef ncell_quantity_constant<CellType, KeyType, DataType>        self_type;
        typedef typename ncell_quantity_interface<CellType>::numeric_type    numeric_type;

      public:
        ncell_quantity_constant(KeyType const & key) : key_(key) {}

        numeric_type eval(CellType const & cell, numeric_type v) const
        {
          return viennadata::access<KeyType, DataType>(key_)(cell);
        }

        numeric_type eval(CellType const & cell, std::vector<numeric_type> const & v) const
        {
          return viennadata::access<KeyType, DataType>(key_)(cell);
        }

        ncell_quantity_interface<CellType> * clone() const { return new self_type(key_); }

      private:
        KeyType key_;
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
                      numeric_type v) const
        {
          return functor_->eval(cell, v);
        }

        NumericT eval(CellType const & cell,
                      std::vector<numeric_type> const & v) const
        {
          return functor_->eval(cell, v);
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
      typedef typename InterfaceType::numeric_type            numeric_type;
    public:

      explicit ncell_quantity(CellType const * cell, detail::ncell_quantity_wrapper<CellType, numeric_type> const & wrapper) : current_cell(cell), accessor(wrapper.clone()) {}

      //template <typename T>
      //explicit cell_quan(T const & t) : current_cell(NULL), accessor( new quan_accessor<CellType, T, numeric_type>() ) {}

      explicit ncell_quantity() : current_cell(NULL) {}

      //interface requirements:
      InterfaceType * clone() const { return new self_type(current_cell, accessor); }
      numeric_type eval(std::vector<numeric_type> const & v) const
      {
        return accessor.eval(*current_cell, v);
      }
      numeric_type eval(numeric_type v) const
      {
        return accessor.eval(*current_cell, v);
      }

      std::string deep_str() const
      {
        std::stringstream ss;
        ss << "cell_quan<" << CellType::tag::dim << ">(" << current_cell << ")";
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
      };

      bool deep_equal(const InterfaceType * other) const
      {
        //TODO: Include comparison of accessor
        return dynamic_cast< const self_type *>(other) != NULL;
      }

      bool shallow_equal(const InterfaceType * other) const
      {
        return dynamic_cast< const self_type *>(other) != NULL;
      }

      InterfaceType * diff(const InterfaceType * diff_var) const
      {
        throw "Cannot differentiate cell_quan!";
        return NULL;
      }


      //additional members:
      void update(CellType const & cell) const
      {
        current_cell = &cell;
      }

      template <typename KeyType>
      void wrap_constant(KeyType const & k)
      {
        detail::ncell_quantity_wrapper<CellType, numeric_type> temp( new detail::ncell_quantity_constant<CellType, KeyType, numeric_type>(k) );
        accessor = temp;
      }

      detail::ncell_quantity_wrapper<CellType, numeric_type> const & wrapper() const { return accessor; }

    private:
      mutable const CellType * current_cell;
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

  //
  // Convenience setter
  //

  //
  // quantity region:
  //

  template <typename KeyType, typename CellType>
  void set_quantity_region(KeyType const & key, CellType const & c, bool b)
  {
    viennadata::access<KeyType, bool >(key)(c) = b;
  }

  namespace detail
  {
    template <typename KeyType, typename DomSegType>
    void set_quantity_region_impl(KeyType const & key,
                                  DomSegType const & seg,
                                  bool b)
    {
      typedef typename DomSegType::config_type  ConfigType;
      typedef typename ConfigType::cell_tag     CellTag;

      typedef typename viennagrid::result_of::ncell<ConfigType, CellTag::dim>::type               CellType;
      typedef typename viennagrid::result_of::const_ncell_range<DomSegType, CellTag::dim>::type  CellContainer;
      typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

      CellContainer cells = viennagrid::ncells(seg);
      for (CellIterator cit  = cells.begin();
                        cit != cells.end();
                      ++cit)
      {
        viennadata::access<KeyType, bool >(key)(*cit) = b;
      }
    }

  }

  template <typename KeyType, typename ConfigType>
  void set_quantity_region(KeyType const & key,
                           viennagrid::segment_t<ConfigType> const & seg,
                           bool b)
  {
    detail::set_quantity_region_impl(key, seg, b);
  }

  template <typename KeyType, typename ConfigType>
  void set_quantity_region(KeyType const & key,
                           viennagrid::domain_t<ConfigType> const & dom,
                           bool b)
  {
    detail::set_quantity_region_impl(key, dom, b);
  }


  //
  // quantity value:
  //

  template <typename KeyType, typename CellType>
  void set_quantity_value(KeyType const & key, CellType const & c, numeric_type val)
  {
    viennadata::access<KeyType, numeric_type >(key)(c) = val;
  }

  namespace detail
  {
    template <typename KeyType, typename DomSegType>
    void set_quantity_value_impl(KeyType const & key,
                                 DomSegType const & seg,
                                 numeric_type val)
    {
      typedef typename DomSegType::config_type  ConfigType;
      typedef typename ConfigType::cell_tag     CellTag;

      typedef typename viennagrid::result_of::ncell<ConfigType, CellTag::dim>::type               CellType;
      typedef typename viennagrid::result_of::const_ncell_range<DomSegType, CellTag::dim>::type  CellContainer;
      typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

      CellContainer cells = viennagrid::ncells(seg);
      for (CellIterator cit  = cells.begin();
                        cit != cells.end();
                      ++cit)
      {
        viennadata::access<KeyType, numeric_type >(key)(*cit) = val;
      }
    }

  }

  template <typename KeyType, typename ConfigType>
  void set_quantity_value(KeyType const & key,
                          viennagrid::segment_t<ConfigType> const & seg,
                          numeric_type val)
  {
    detail::set_quantity_value_impl(key, seg, val);
  }

  template <typename KeyType, typename ConfigType>
  void set_quantity_value(KeyType const & key,
                          viennagrid::domain_t<ConfigType> const & dom,
                          numeric_type val)
  {
    detail::set_quantity_value_impl(key, dom, val);
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
