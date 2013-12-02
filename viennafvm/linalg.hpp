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
#include <cmath>
#include <cassert>


/** @file viennafvm/linalg.hpp
    @brief Some basic linear algebra functionality
*/

namespace viennafvm
{

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
