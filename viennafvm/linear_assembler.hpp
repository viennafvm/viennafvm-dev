#ifndef VIENNAFVM_PDE_ASSEMBLER_HPP
#define VIENNAFVM_PDE_ASSEMBLER_HPP

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


// *** local includes
//
#include "viennafvm/integral_form.hpp"
#include "viennafvm/extract_integrals.hpp"
#include "viennafvm/linear_pde_system.hpp"
#include "viennafvm/mapping.hpp"
#include "viennafvm/util.hpp"
#include "viennafvm/flux.hpp"
#include "viennafvm/ncell_quantity.hpp"

#include "viennagrid/forwards.hpp"
#include "viennagrid/algorithm/voronoi.hpp"
#include "viennagrid/algorithm/volume.hpp"
#include "viennagrid/algorithm/centroid.hpp"
#include "viennagrid/mesh/coboundary_iteration.hpp"

#include "viennamath/manipulation/eval.hpp"
#include "viennamath/manipulation/diff.hpp"

#include "viennadata/api.hpp"

//#define VIENNAFVMDEBUG

namespace viennafvm
{
  template <typename CellType, typename FacetType, typename CellValueAccessorType, typename FacetValueAccessorType, typename FacetDistanceAccessorType>
  void compute_gradients_for_cell(CellType const & inner_cell, FacetType const & facet, CellType const & outer_cell,
                                  CellValueAccessorType const cell_value_accessor, FacetValueAccessorType facet_value_accessor, FacetDistanceAccessorType const facet_distance_accessor)
  {
    double value_inner = cell_value_accessor(inner_cell);
    double value_outer = cell_value_accessor(outer_cell);
    double distance    = facet_distance_accessor(facet);

    facet_value_accessor(facet) = (value_outer - value_inner) / distance;
  }




  class linear_assembler
  {
    public:

      /** @brief  Assembles the full PDE system into the same matrix */
      template <typename LinPdeSysT,
                typename SegmentT,
                typename QuantityContainer,
                typename MatrixT,
                typename VectorT>
      void operator()(LinPdeSysT const  & pde_system,
                      SegmentT   const  & segment,
                      QuantityContainer & quantities,
                      MatrixT           & system_matrix,
                      VectorT           & load_vector)
      {
        std::size_t map_index = viennafvm::create_mapping(pde_system, segment, quantities);

        system_matrix.clear();
        system_matrix.resize(map_index, map_index, false);
        load_vector.clear();
        load_vector.resize(map_index);


        for (std::size_t pde_index = 0; pde_index < pde_system.size(); ++pde_index)
        {
#ifdef VIENNAFVM_DEBUG
          std::cout << std::endl;
          std::cout << "//" << std::endl;
          std::cout << "//   Equation " << pde_index << std::endl;
          std::cout << "//" << std::endl;
#endif
          assemble(pde_system, pde_index,
                   segment, quantities,
                   system_matrix, load_vector);

        } // for pde_index
      } // functor


      /** @brief  Assembles one PDE out of the PDE system into the matrix */
      template <typename LinPdeSysT,
                typename SegmentT,
                typename QuantityContainer,
                typename MatrixT,
                typename VectorT>
      void operator()(LinPdeSysT const  & pde_system,
                      std::size_t         pde_index,
                      SegmentT   const  & segment,
                      QuantityContainer & quantities,
                      MatrixT           & system_matrix,
                      VectorT           & load_vector)
      {
        std::size_t map_index = viennafvm::create_mapping(pde_system, pde_index, segment, quantities);

        system_matrix.clear();
        system_matrix.resize(map_index, map_index, false);
        load_vector.clear();
        load_vector.resize(map_index);


#ifdef VIENNAFVM_DEBUG
        std::cout << std::endl;
        std::cout << "//" << std::endl;
        std::cout << "//   Equation " << pde_index << std::endl;
        std::cout << "//" << std::endl;
#endif
        assemble(pde_system, pde_index,
                 segment, quantities,
                 system_matrix, load_vector);

      } // functor


    private:

      template <typename PDESystemType,
                typename SegmentT,
                typename QuantityContainer,
                typename MatrixT,
                typename VectorT>
      void assemble(PDESystemType const & pde_system,
                    std::size_t           pde_index,
                    SegmentT      const & segment,
                    QuantityContainer   & quantities,
                    MatrixT             & system_matrix,
                    VectorT             & load_vector)
      {
        typedef viennamath::equation                          equ_type;
        typedef viennamath::expr                              expr_type;
        typedef typename expr_type::interface_type            interface_type;

        typedef typename viennagrid::result_of::cell_tag<SegmentT>::type CellTag;
        typedef typename viennagrid::result_of::facet_tag<CellTag>::type FacetTag;

        typedef typename viennagrid::result_of::element<SegmentT, FacetTag>::type                FacetType;
        typedef typename viennagrid::result_of::element<SegmentT, CellTag  >::type               CellType;

        typedef typename viennagrid::result_of::const_element_range<SegmentT, CellTag>::type     CellContainer;
        typedef typename viennagrid::result_of::iterator<CellContainer>::type                    CellIterator;

        typedef typename viennagrid::result_of::const_element_range<CellType, FacetTag>::type    FacetOnCellContainer;
        typedef typename viennagrid::result_of::iterator<FacetOnCellContainer>::type             FacetOnCellIterator;

        typedef typename QuantityContainer::value_type     QuantityType;

        viennamath::equation          const & pde         = pde_system.pde(pde_index);
        viennamath::function_symbol   const & u           = pde_system.unknown(pde_index)[0];
        viennafvm::linear_pde_options const & pde_options = pde_system.option(pde_index);

        QuantityType const & quan = quantities.at(u.id());


#ifdef VIENNAFVM_DEBUG
        std::cout << " - Strong form: " << pde << std::endl;
#endif

        equ_type integral_form = viennafvm::make_integral_form( pde );

#ifdef VIENNAFVM_DEBUG
        std::cout << " - Integral form: " << integral_form << std::endl;
#endif

        //
        // Preprocess symbolic representation:
        //

        //Note: Assuming that LHS holds all matrix terms, while RHS holds all load vector terms
        expr_type  partial_omega_integrand = extract_surface_integrand<CellType>(quantities, integral_form.lhs(), u);
        expr_type   matrix_omega_integrand = extract_volume_integrand<CellType>(quantities, integral_form.lhs(), u);
        expr_type      rhs_omega_integrand = extract_volume_integrand<CellType>(quantities, integral_form.rhs(), viennamath::function_symbol(1337)); //workaround for substituting *all* function symbols
        expr_type  stabilization_integrand = prepare_for_evaluation<CellType>(quantities, pde_options.damping_term(), u);

        expr_type substituted_matrix_omega_integrand_no_diff  = viennamath::diff(matrix_omega_integrand, u);
        viennamath::rt_manipulation_wrapper<interface_type> wrapped_function_symbol_replacer( new detail::function_symbol_replacer<QuantityContainer, CellType, interface_type>(quantities, viennamath::function_symbol(1337)) );
        viennamath::rt_expr<interface_type> replaced_integrand(substituted_matrix_omega_integrand_no_diff.get()->recursive_manipulation( wrapped_function_symbol_replacer ));

        expr_type substituted_matrix_omega_integrand = viennamath::simplify(replaced_integrand);

#ifdef VIENNAFVM_DEBUG
        std::cout << " - Surface integrand for matrix: " << partial_omega_integrand << std::endl;
        std::cout << " - Volume integrand for matrix:  " <<  matrix_omega_integrand << std::endl;
        std::cout << " - Stabilization for matrix:     " << stabilization_integrand << std::endl;
        std::cout << " - Volume integrand for rhs:     " <<     rhs_omega_integrand << std::endl;
        std::cout << " - Substituted volume integrand for matrix: " << substituted_matrix_omega_integrand << std::endl;
#endif

        viennafvm::flux_handler<QuantityContainer, CellType, FacetType, interface_type>  flux(quantities, partial_omega_integrand, u);


        std::vector<double> p(3); //dummy vector for evaluation

        typedef viennadata::storage<>   StorageType;
        StorageType storage;

        //
        // Preprocess domain
        //
        setup(segment, storage);

        typename viennadata::result_of::accessor<StorageType, viennafvm::facet_area_key, double, FacetType>::type facet_area_accessor =
            viennadata::make_accessor(storage, viennafvm::facet_area_key());

        typename viennadata::result_of::accessor<StorageType, facet_distance_key, double, FacetType>::type facet_distance_accessor =
            viennadata::make_accessor(storage, facet_distance_key());

        //
        // Actual assembly:
        //
        CellContainer cells(segment);
        for (CellIterator cit = cells.begin(); cit != cells.end(); ++cit)
        {
          long row_index = quan.get_unknown_index(*cit);

          if (row_index < 0)
            continue;

          //
          // Boundary integral terms:
          //
          FacetOnCellContainer facets_on_cell(*cit);
          for (FacetOnCellIterator focit  = facets_on_cell.begin();
                                   focit != facets_on_cell.end();
                                 ++focit)
          {
            CellType const * other_cell = util::other_cell_of_facet(*focit, *cit, segment);

            if (other_cell)
            {
              long col_index = quan.get_unknown_index(*other_cell);

              double effective_facet_area = facet_area_accessor(*focit);
              double distance             = facet_distance_accessor(*focit);

              if (quan.get_boundary_type(*other_cell) == viennafvm::BOUNDARY_DIRICHLET)
              {
                double boundary_value = quan.get_boundary_value(*other_cell);
                double current_value  = quan.get_value(*cit);

                // updates are homogeneous, hence no direct contribution to RHS here. Might change later when boundary values are slowly increased.
                load_vector(row_index)              -= flux.out(*cit, *focit, *other_cell, distance) * effective_facet_area * (boundary_value - current_value);
                system_matrix(row_index, row_index) -= flux.in(*cit, *focit, *other_cell, distance) * effective_facet_area;

                load_vector(row_index) -= flux.out(*cit, *focit, *other_cell, distance) * effective_facet_area * current_value;
                load_vector(row_index) += flux.in(*cit, *focit, *other_cell, distance) * effective_facet_area * quan.get_value(*cit);
              }
              else if (col_index >= 0)
              {
                system_matrix(row_index, col_index) += flux.out(*cit, *focit, *other_cell, distance) * effective_facet_area;
                system_matrix(row_index, row_index) -= flux.in(*cit, *focit, *other_cell, distance) * effective_facet_area;

                load_vector(row_index) -= flux.out(*cit, *focit, *other_cell, distance) * effective_facet_area * quan.get_value(*other_cell);
                load_vector(row_index) += flux.in(*cit, *focit, *other_cell, distance) * effective_facet_area * quan.get_value(*cit);
              }
              // else: nothing to do because other cell is not considered for this quantity

            }
          }

          //
          // Volume terms
          //
          double cell_volume      = viennagrid::volume(*cit);

          // Matrix (including residual contributions)
          viennamath::rt_traversal_wrapper<interface_type> cell_updater(new detail::ncell_updater<CellType, interface_type>(*cit));
          substituted_matrix_omega_integrand.get()->recursive_traversal(cell_updater);
          system_matrix(row_index, row_index) += viennamath::eval(substituted_matrix_omega_integrand, p) * cell_volume;
          load_vector(row_index) -= viennamath::eval(substituted_matrix_omega_integrand, p) * cell_volume * quan.get_value(*cit);

          stabilization_integrand.get()->recursive_traversal(cell_updater);
          system_matrix(row_index, row_index) += viennamath::eval(stabilization_integrand, p) * cell_volume;

          // RHS
          rhs_omega_integrand.get()->recursive_traversal(cell_updater);
          load_vector(row_index) += viennamath::eval(rhs_omega_integrand, p) * cell_volume;
          //std::cout << "Writing " << viennamath::eval(rhs_omega_integrand, p) << " * " << cell_volume << " to rhs at " << row_index << std::endl;

        } // for cells

      } // assemble

      template <typename SegmentT, typename StorageType>
      void setup(SegmentT const & segment, StorageType & storage)
      {
        typedef typename viennagrid::result_of::cell_tag<SegmentT>::type CellTag;
        typedef typename viennagrid::result_of::facet_tag<CellTag>::type FacetTag;

        typedef typename viennagrid::result_of::point<SegmentT>::type                  PointType;
        typedef typename viennagrid::result_of::element<SegmentT, FacetTag>::type  FacetType;

        typedef typename viennagrid::result_of::const_element_range<SegmentT, FacetTag>::type  FacetContainer;
        typedef typename viennagrid::result_of::iterator<FacetContainer>::type                     FacetIterator;

        typedef typename viennagrid::result_of::const_coboundary_range<SegmentT, FacetType, CellTag>::type CellOnFacetRange;
//         typedef typename viennagrid::result_of::const_element_range<FacetType, CellTag>::type   CellOnFacetRange;
        typedef typename viennagrid::result_of::iterator<CellOnFacetRange>::type                   CellOnFacetIterator;



        typename viennadata::result_of::accessor<StorageType, viennafvm::facet_area_key, double, FacetType>::type facet_area_accessor =
            viennadata::make_accessor(storage, viennafvm::facet_area_key());

        typename viennadata::result_of::accessor<StorageType, viennafvm::facet_distance_key, double, FacetType>::type facet_distance_accessor =
            viennadata::make_accessor(storage, viennafvm::facet_distance_key());


        FacetContainer facets(segment);
        for (FacetIterator fit  = facets.begin();
                           fit != facets.end();
                         ++fit)
        {

          CellOnFacetRange    cells = viennagrid::coboundary_elements<FacetType, CellTag>(segment, fit.handle());

          if (cells.size() == 2)
          {
            CellOnFacetIterator cofit    = cells.begin();

            PointType centroid_1         = viennagrid::centroid(*cofit); ++cofit;
            PointType centroid_2         = viennagrid::centroid(*cofit);
            PointType center_connection  = centroid_1 - centroid_2;
            PointType outer_normal       = util::unit_outer_normal(*fit, *cofit, viennagrid::default_point_accessor(segment)); //note: consistent orientation of center_connection and outer_normal is important here!

            double center_connection_len = viennagrid::norm(center_connection);
            double effective_facet_ratio = viennagrid::inner_prod(center_connection, outer_normal) / center_connection_len;  // inner product of unit vectors
            double effective_facet_area  = viennagrid::volume(*fit) * effective_facet_ratio;

            facet_area_accessor(*fit) = effective_facet_area;
            facet_distance_accessor(*fit) = center_connection_len;
          }
        }
    }
  };

  template <typename InterfaceType, typename SegmentT, typename MatrixT, typename VectorT>
  void assemble_pde( viennafvm::linear_pde_system<InterfaceType>  & pde_system,
                     SegmentT & segment,
                     MatrixT  & matrix,
                     VectorT  & rhs
                    )
  {
     viennafvm::linear_assembler()(pde_system, segment, matrix, rhs);
  }
}

#endif







