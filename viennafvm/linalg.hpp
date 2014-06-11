#ifndef VIENNAFVM_LINALG_HPP
#define VIENNAFVM_LINALG_HPP

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
#include <map>
#include <iostream>
#include <cmath>
#include <cassert>


/** @file viennafvm/linalg.hpp
    @brief Some basic linear algebra functionality
*/

namespace viennafvm
{


  /** @brief Sparse matrix class based on a vector of binary trees for holding the entries.
    *
    * Introduced to overcome the scalability issues of ublas::compressed_matrix
    */
  template <typename NumericT>
  class sparse_matrix
  {
    typedef sparse_matrix<NumericT>  self_type;

  public:
    typedef std::size_t    size_type;
    typedef size_type      SizeType;

  private:
    typedef std::map<size_type, NumericT>     MapType;
    typedef std::vector<MapType>              DataContainerType;

  public:
    typedef typename MapType::iterator                      iterator2;
    typedef typename MapType::const_iterator          const_iterator2;

    typedef MapType     row_type;
    typedef MapType     RowType;

    sparse_matrix() : zero_(0) {}
    sparse_matrix(size_type num_rows, size_type num_cols) : data_(num_rows), zero_(0)
    {
      assert(num_cols == num_rows && bool("Only square matrices supported."));
      (void)num_cols;
    }

    size_type size1() const { return data_.size(); }
    size_type size2() const { return data_.size(); }

    void resize(size_type num_rows, size_type num_cols, bool preserve = false)
    {
      (void)num_cols;
      if (!preserve)
        data_.clear();
      data_.resize(num_rows);

    }

    void clear() { data_.clear(); }

    /** @brief Non-const access to the entries of the matrix. Use this for assembly. */
    NumericT & operator()(size_type i, size_type j)
    {
      assert(i < data_.size() && j < data_.size() && bool("Access to matrix out of range!"));
      return data_.at(i)[j];
    }

    NumericT const & operator()(size_type i, size_type j) const
    {
      assert(i < data_.size() && j < data_.size() && bool("Access to matrix out of range!"));
      const_iterator2 it = data_.at(i).find(j);
      if (it == data_.at(i).end())
        return zero_;
      else
        return it->second;
    }

    MapType const & row(size_t i) const { return data_.at(i); }
    MapType       & row(size_t i)       { return data_.at(i);   }

    size_type nnz() const
    {
      size_type result = 0;
      for (size_type i=0; i<data_.size(); ++i)
        result += data_[i].size();
      return result;
    }

    self_type trans() const
    {
      self_type A_trans(data_.size(), data_.size());

      for (size_type i=0; i<data_.size(); ++i)
      {
        for (const_iterator2 it = data_[i].begin(); it != data_[i].end(); ++it)
          A_trans(it->first, i) = it->second;
      }
      return A_trans;
    }

    self_type operator*(NumericT factor) const
    {
      self_type result = *this;
      for (size_type i=0; i<data_.size(); ++i)
      {
        for (const_iterator2 it = data_[i].begin(); it != data_[i].end(); ++it)
          result(it->first, i) *= factor;
      }
      return result;
    }

    self_type & operator+=(self_type const & B)
    {
      for (size_type i=0; i<B.size1(); ++i)
      {
        row_type row_i = B.row(i);
        for (const_iterator2 it = row_i.begin(); it != row_i.end(); ++it)
          data_[i][it->first] += it->second;
      }
      return *this;
    }

    DataContainerType const & get() const { return data_; }

  private:
    DataContainerType  data_;
    NumericT           zero_;   //helper member for implementing operator() const
  };


  /** @brief Normalizes an equation system such that all diagonal entries are non-negative, and such that all 2-norms of the rows are unity. */
  template <typename NumericT, typename VectorT>
  VectorT row_normalize_system(sparse_matrix<NumericT> & A, VectorT & b)
  {
    typedef typename sparse_matrix<NumericT>::row_type     RowType;
    typedef typename sparse_matrix<NumericT>::iterator2    AlongRowIterator;

    VectorT scale_factors(b.size());

    for (std::size_t i=0; i<A.size1(); ++i)
    {
      RowType & row_i = A.row(i);

      // obtain current norm of row:
      double row_norm = 0.0;
      for (AlongRowIterator iter  = row_i.begin();
                            iter != row_i.end();
                          ++iter)
      {
        row_norm += iter->second * iter->second;
      }

      row_norm = sqrt(row_norm);

      // normalize such that diagonal entry becomes positive:
      if (A(i, i) < 0.0)
        row_norm *= -1.0;

      // scale row accordingly:
      for (AlongRowIterator iter  = row_i.begin();
                            iter != row_i.end();
                          ++iter)
      {
        iter->second /= row_norm;
      }
      b[i] /= row_norm;

      // remember scaling factor:
      scale_factors[i] = 1.0 / row_norm;
    }

    return scale_factors;
  }

  /** @brief Converts a sparse matrix to a string. Handy debugging facility. */
  template <typename NumericT>
  std::ostream & operator<<(std::ostream & os, sparse_matrix<NumericT> const & A)
  {
    typedef typename sparse_matrix<NumericT>::const_iterator2     AlongRowIterator;

    for (std::size_t i=0; i<A.size1(); ++i)
    {
      os << std::endl << "Row " << i << ": ";
      for (AlongRowIterator iter  = A.row(i).begin();
                            iter != A.row(i).end();
                          ++iter)
      {
        os << "(" << iter->first << ", " << iter->second << "),  ";
      }
    }

    return os;
  }


  /** @brief Computes A * x  for a sparse A and a vector x.
    *
    * Does not apply expression template fancyness, but this is not a performance-critical routine anyway...
   */
  template <typename NumericT, typename VectorType>
  VectorType prod(sparse_matrix<NumericT> const & A,
                  VectorType const & x)
  {
    typedef typename sparse_matrix<NumericT>::const_iterator2    Iterator2;
    typedef typename sparse_matrix<NumericT>::row_type           RowType;

    VectorType result(A.size1());
    for (std::size_t i=0; i<A.size1(); ++i)
    {
      RowType const & row_i = A.row(i);
      double val = 0.0;
      for (Iterator2 iter2  = row_i.begin();
                     iter2 != row_i.end();
                   ++iter2)
      {
        val += iter2->second * x[iter2->first];
      }
      result[i] = val;
    }

    return result;
  }

  /////////////////////////////////////////////////////////////////////////

  template <typename NumericT>
  class dense_matrix
  {
  public:
    dense_matrix(std::size_t rows, std::size_t cols) : num_rows_(rows), num_cols_(cols), data_(rows * cols) {}

    NumericT const & operator()(std::size_t i, std::size_t j) const { return data_[i*num_cols_ + j]; }
    NumericT       & operator()(std::size_t i, std::size_t j)       { return data_[i*num_cols_ + j]; }

    std::size_t rows() const { return num_rows_; }
    std::size_t cols() const { return num_cols_; }

  private:
    std::size_t num_rows_;
    std::size_t num_cols_;
    std::vector<NumericT> data_;
  };

  template <typename NumericT>
  std::vector<NumericT> solve(dense_matrix<NumericT> const & A, std::vector<NumericT> const & b)
  {
    assert(A.rows() == A.cols() && bool("Square input matrix required!"));

    std::vector<NumericT> result(A.rows());

    if (A.rows() == 1)
    {
      result[0] = b[0] / A(0,0);
    }
    else if (A.rows() == 2)
    {
      //
      // A = [a, b; c, d] results in A^{-1} = [d, -b; -c, a] / det(A)
      //

      NumericT det_A = A(0,0) * A(1,1) - A(1,0)*A(0,1);

      assert(det_A != 0.0 && bool("Singular matrix detected!"));

      result[0] = (b[0] * A(1,1) - b[1] * A(0,1)) / det_A;
      result[1] = (b[0] * A(1,0) - b[1] * A(0,0)) / det_A;
    }
    else if (A.rows() == 3)
    {
      //
      // Use Cramer's rule
      //
      double det_A =   A(0,0)*A(1,1)*A(2,2) + A(0,1)*A(1,2)*A(3,0) + A(0,2)*A(1,0)*A(2,1)
                     - A(2,0)*A(1,1)*A(0,2) - A(2,1)*A(1,2)*A(0,0) - A(2,2)*A(1,0)*A(0,1);

      // compute entries of inverse matrix from cofactors
      double inv_A_00 =  A(1,1)*A(2,2) - A(2,1)*A(2,1);
      double inv_A_01 = -A(0,1)*A(2,2) + A(2,1)*A(0,2);
      double inv_A_02 =  A(0,1)*A(1,2) - A(1,1)*A(0,2);

      double inv_A_10 = -A(1,0)*A(2,2) + A(2,1)*A(1,2);
      double inv_A_11 =  A(0,0)*A(2,2) - A(2,0)*A(0,2);
      double inv_A_12 = -A(0,0)*A(1,2) + A(1,0)*A(0,2);

      double inv_A_20 =  A(1,0)*A(2,1) - A(2,0)*A(1,1);
      double inv_A_21 = -A(0,0)*A(2,1) + A(2,0)*A(0,1);
      double inv_A_22 =  A(0,0)*A(1,1) - A(1,0)*A(0,1);

      result[0] = (inv_A_00 * b[0] + inv_A_01 * b[1] + inv_A_02 * b[2]) / det_A;
      result[1] = (inv_A_10 * b[0] + inv_A_11 * b[1] + inv_A_12 * b[2]) / det_A;
      result[2] = (inv_A_20 * b[0] + inv_A_21 * b[1] + inv_A_22 * b[2]) / det_A;
    }
    else
    {
      throw "Implementation of systems larger than 3-by-3 not implemented!";
    }

    return result;
  }

} //namespace viennafvm
#endif
