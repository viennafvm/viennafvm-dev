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

#include "viennagrid/forwards.h"
#include "viennagrid/algorithm/voronoi.hpp"
#include "viennagrid/algorithm/volume.hpp"
#include "viennagrid/algorithm/centroid.hpp"

#include "viennamath/manipulation/eval.hpp"
#include "viennamath/manipulation/diff.hpp"

#include "viennadata/api.hpp"

//#define VIENNAFVMDEBUG

namespace viennafvm
{
  template <typename CellType, typename FacetType>
  void compute_gradients_for_cell(CellType const & inner_cell, FacetType const & facet, CellType const & outer_cell, long id)
  {
    current_iterate_key key(id);

    double value_inner = viennadata::access<current_iterate_key, double>(key)(inner_cell);
    double value_outer = viennadata::access<current_iterate_key, double>(key)(outer_cell);
    double distance    = viennadata::access<facet_distance_key, double>(facet_distance_key())(facet);

    viennadata::access<current_iterate_key, double>(key)(facet) = (value_outer - value_inner) / distance;
  }




  class linear_assembler
  {
    public:
      template <typename LinPdeSysT,
                typename SegmentT,
                typename MatrixT,
                typename VectorT>
      void operator()(LinPdeSysT const & pde_system,
                      SegmentT   const & segment,
                      MatrixT          & system_matrix,
                      VectorT          & load_vector)
      {
        typedef typename SegmentT::config_type                config_type;
        typedef viennamath::equation                          equ_type;
        typedef viennamath::expr                              expr_type;
        typedef typename expr_type::interface_type            interface_type;
        typedef typename expr_type::numeric_type              numeric_type;

        typedef typename SegmentT::config_type              Config;
        typedef typename Config::cell_tag                   CellTag;
        typedef typename viennagrid::result_of::ncell<Config, CellTag::dim-1>::type                FacetType;
        typedef typename viennagrid::result_of::ncell<Config, CellTag::dim  >::type                CellType;

        typedef typename viennagrid::result_of::const_ncell_range<SegmentT, CellTag::dim>::type    CellContainer;
        typedef typename viennagrid::result_of::iterator<CellContainer>::type                      CellIterator;

        typedef typename viennagrid::result_of::const_ncell_range<CellType, CellTag::dim-1>::type  FacetOnCellContainer;
        typedef typename viennagrid::result_of::iterator<FacetOnCellContainer>::type               FacetOnCellIterator;

         std::size_t map_index = viennafvm::create_mapping(pde_system, segment);

         std::cout << "Total number of unknowns: " << map_index << std::endl;
         system_matrix.clear();
         system_matrix.resize(map_index, map_index, false);
         load_vector.clear();
         load_vector.resize(map_index);


         for (std::size_t pde_index = 0; pde_index < pde_system.size(); ++pde_index)
         {
           std::cout << std::endl;
           std::cout << "//" << std::endl;
           std::cout << "//   Equation " << pde_index << std::endl;
           std::cout << "//" << std::endl;
           typedef typename LinPdeSysT::mapping_key_type   MappingKeyType;
           MappingKeyType  map_key(pde_system.option(pde_index).data_id());

           typedef typename LinPdeSysT::boundary_key_type  BoundaryKeyType;
           BoundaryKeyType bnd_key(pde_system.option(pde_index).data_id());


        #ifdef VIENNAFVM_DEBUG
           std::cout << " - Strong form: " << pde_system.pde(pde_index) << std::endl;
        #endif

           equ_type integral_form = viennafvm::make_integral_form( pde_system.pde(pde_index) );

        #ifdef VIENNAFVM_DEBUG
           std::cout << " - Integral form: " << integral_form << std::endl;
        #endif


            //
            // Preprocess symbolic representation:
            //

            //Note: Assuming that LHS holds all matrix terms, while RHS holds all load vector terms
            expr_type  partial_omega_integrand = extract_surface_integrand<FacetType>(integral_form.lhs(), pde_system.unknown(pde_index)[0]);
            expr_type   matrix_omega_integrand = extract_volume_integrand<CellType>(integral_form.lhs(), pde_system.unknown(pde_index)[0]);
            expr_type      rhs_omega_integrand = extract_volume_integrand<CellType>(integral_form.rhs(), pde_system.unknown(pde_index)[0]);
            expr_type  stabilization_integrand = prepare_for_evaluation<CellType>(pde_system.option(pde_index).damping_term(), pde_system.unknown(pde_index)[0]);

            std::cout << " - Surface integrand for matrix: " << partial_omega_integrand << std::endl;
            std::cout << " - Volume integrand for matrix:  " <<  matrix_omega_integrand << std::endl;
            std::cout << " - Stabilization for matrix:     " << stabilization_integrand << std::endl;
            std::cout << " - Volume integrand for rhs:     " <<     rhs_omega_integrand << std::endl;

            viennafvm::flux_handler<CellType, FacetType, interface_type>  flux(partial_omega_integrand, pde_system.unknown(pde_index)[0]);

            expr_type substituted_matrix_omega_integrand  = viennamath::diff(matrix_omega_integrand, pde_system.unknown(pde_index)[0]);

            std::vector<double> p(3); //dummy vector for evaluation

            //
            // Preprocess domain
            //
            setup(segment);

            //
            // Actual assembly:
            //
            CellContainer cells = viennagrid::ncells(segment);
            for (CellIterator cit = cells.begin(); cit != cells.end(); ++cit)
            {
              long row_index = viennadata::access<MappingKeyType, long>(map_key)(*cit);

              if (row_index < 0)
                continue;

              //flux.set_inner_cell(*cit);

              //
              // Boundary integral terms:
              //
              FacetOnCellContainer facets_on_cell = viennagrid::ncells(*cit);
              for (FacetOnCellIterator focit  = facets_on_cell.begin();
                                       focit != facets_on_cell.end();
                                     ++focit)
              {
                CellType const * other_cell = util::other_cell_of_facet(*focit, *cit, segment);

                if (other_cell)
                {
                  long col_index = viennadata::access<MappingKeyType, long>(map_key)(*other_cell);
                  double effective_facet_area = viennadata::access<viennafvm::facet_area_key, double>()(*focit);

                  for (std::size_t i=0; i<pde_system.size(); ++i)
                    compute_gradients_for_cell(*cit, *focit, *other_cell, pde_system.unknown(i)[0].id());

                  if (col_index == viennafvm::DIRICHLET_BOUNDARY)
                  {
                    double boundary_value = viennadata::access<BoundaryKeyType, double>(bnd_key)(*other_cell);

                    // updates are homogeneous, hence no direct contribution to RHS here. Might change later when boundary values are slowly increased.
                    system_matrix(row_index, row_index) -= flux.in(*cit, *focit, *other_cell) * effective_facet_area;

                    load_vector(row_index) -= flux.out(*cit, *focit, *other_cell) * effective_facet_area * boundary_value;
                    load_vector(row_index) += flux.in(*cit, *focit, *other_cell) * effective_facet_area * get_current_iterate(*cit, pde_system.unknown(pde_index)[0]);
                  }
                  else if (col_index >= 0)
                  {
                    system_matrix(row_index, col_index) += flux.out(*cit, *focit, *other_cell) * effective_facet_area;
                    system_matrix(row_index, row_index) -= flux.in(*cit, *focit, *other_cell) * effective_facet_area;

                    load_vector(row_index) -= flux.out(*cit, *focit, *other_cell) * effective_facet_area * get_current_iterate(*other_cell, pde_system.unknown(pde_index)[0]);
                    load_vector(row_index) += flux.in(*cit, *focit, *other_cell) * effective_facet_area * get_current_iterate(*cit, pde_system.unknown(pde_index)[0]);
                  }
                  // else: nothing to do because other cell is not considered for this quantity

                }
              }

              //
              // Volume terms
              //
              double cell_volume      = viennagrid::volume(*cit);

              // Matrix (including residual contributions)
              system_matrix(row_index, row_index) += viennamath::eval(substituted_matrix_omega_integrand, p) * cell_volume;
              load_vector(row_index) -= viennamath::eval(substituted_matrix_omega_integrand, p) * cell_volume * get_current_iterate(*cit, pde_system.unknown(pde_index)[0]);

              system_matrix(row_index, row_index) += viennamath::eval(stabilization_integrand, p) * cell_volume;

              // RHS
              load_vector(row_index) += viennamath::eval(rhs_omega_integrand, p) * cell_volume;
              //std::cout << "Writing " << viennamath::eval(omega_integrand, p) << " * " << cell_volume << " to rhs at " << row_index << std::endl;

            } // for cells

          } // for pde_index
      } // functor

    private:

      template <typename SegmentT>
      void setup(SegmentT const & segment)
      {
        typedef typename SegmentT::config_type                config_type;
        typedef viennamath::equation                          equ_type;
        typedef viennamath::expr                              expr_type;
        typedef typename expr_type::numeric_type              numeric_type;

        typedef typename SegmentT::config_type              Config;
        typedef typename Config::cell_tag                   CellTag;
        typedef typename viennagrid::result_of::point<Config>::type                  PointType;
        typedef typename viennagrid::result_of::ncell<Config, CellTag::dim-1>::type  FacetType;
        typedef typename viennagrid::result_of::ncell<Config, CellTag::dim  >::type  CellType;

        typedef typename viennagrid::result_of::const_ncell_range<SegmentT, CellTag::dim-1>::type  FacetContainer;
        typedef typename viennagrid::result_of::iterator<FacetContainer>::type                     FacetIterator;

        typedef typename viennagrid::result_of::const_ncell_range<FacetType, CellTag::dim>::type   CellOnFacetRange;
        typedef typename viennagrid::result_of::iterator<CellOnFacetRange>::type                   CellOnFacetIterator;

        FacetContainer facets = viennagrid::ncells(segment);
        for (FacetIterator fit  = facets.begin();
                           fit != facets.end();
                         ++fit)
        {

          CellOnFacetRange    cells = viennagrid::ncells<CellTag::dim>(*fit, segment);

          if (cells.size() == 2)
          {
            CellOnFacetIterator cofit    = cells.begin();

            PointType centroid_1         = viennagrid::centroid(*cofit); ++cofit;
            PointType centroid_2         = viennagrid::centroid(*cofit);
            PointType center_connection  = centroid_1 - centroid_2;
            PointType outer_normal       = util::unit_outer_normal(*fit, *cofit); //note: consistent orientation of center_connection and outer_normal is important here!

            double center_connection_len = viennagrid::norm(center_connection);
            double effective_facet_ratio = viennagrid::inner_prod(center_connection, outer_normal) / center_connection_len;  // inner product of unit vectors
            double effective_facet_area  = viennagrid::volume(*fit) * effective_facet_ratio;

            viennadata::access<viennafvm::facet_area_key,     double>()(*fit) = effective_facet_area;
            viennadata::access<viennafvm::facet_distance_key, double>()(*fit) = center_connection_len;
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







