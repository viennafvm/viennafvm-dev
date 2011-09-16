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

#ifndef VIENNAFVM_POISSON_ASSEMBLER_HPP
#define VIENNAFVM_POISSON_ASSEMBLER_HPP

#include <iostream>
#include <fstream>
#include <vector>

#include "viennafvm/forwards.h"
#include "viennafvm/common.hpp"

#include "viennagrid/domain.hpp"
#include "viennagrid/iterators.hpp"
#include "viennagrid/algorithm/volume.hpp"


//
// [KR]
// General notes: This is in fact still far from a good finite volume implementation - ViennaMath still needs to be integrated here.
//                For the moment, this is therefore nothing but a draft for the internal assembly iterations.
//

namespace viennafvm
{

    /** @brief Computes the potential update in a Gummel iteration for given electron and hole concentrations */
    struct poisson_assembler
    {

      template <typename DeviceType, typename MatrixType, typename VectorType>
      void operator()(DeviceType & domain,  
                      MatrixType & system_matrix,  
                      VectorType & rhs)
      {
        typedef typename DeviceType::config_type           Config;
        typedef typename Config::cell_tag                  CellTag;
        typedef typename viennagrid::result_of::ncell<Config, 0>::type                         VertexType;
        typedef typename viennagrid::result_of::ncell<Config, 1>::type                         EdgeType;
        typedef typename viennagrid::result_of::ncell<Config, CellTag::dim>::type   CellType;
        
        typedef typename viennagrid::result_of::ncell_range<DeviceType, 0>::type     VertexContainer;
        typedef typename viennagrid::result_of::iterator<VertexContainer>::type          VertexIterator;

        typedef typename viennagrid::result_of::ncell_range<VertexType, 1>::type     EdgeOnVertexContainer;
        typedef typename viennagrid::result_of::iterator<EdgeOnVertexContainer>::type    EdgeOnVertexIterator;
        
        typedef typename viennagrid::result_of::ncell_range<EdgeType, 0>::type       VertexOnEdgeContainer;
        typedef typename viennagrid::result_of::iterator<VertexOnEdgeContainer>::type    VertexOnEdgeIterator;
        
        typedef viennafvm::boundary_key                             BoundaryKeyType;
        typedef viennafvm::mapping_key                              MappingKeyType;
        
        
        std::cout << "* poisson_assembler::operator(): Create Mapping:" << std::endl;
      
        size_t map_index = 0;
        BoundaryKeyType bnd_key(0);
        MappingKeyType map_key(0);
        
        VertexContainer vertices = viennagrid::ncells<0>(domain);
        for (VertexIterator vit = vertices.begin();
            vit != vertices.end();
            ++vit)
        {  
          if (viennadata::access<BoundaryKeyType, bool>(bnd_key)(*vit))
          {
            //std::cout << "boundary vertex" << std::endl;
            viennadata::access<MappingKeyType, long>(map_key)(*vit) = -1;
          }
          else
          {
            //std::cout << "interior vertex" << std::endl;
            viennadata::access<MappingKeyType, long>(map_key)(*vit) = map_index;
            map_index += 1;
          }
        }
        std::cout << "---------------------------" << std::endl;
        
        std::cout << "* poisson_assembler::operator(): Assigned degrees of freedom: " << map_index << std::endl;
        
        //resize global system matrix and load vector if needed:
        if (map_index > system_matrix.size1())
        {
          std::cout << "Resizing system matrix..." << std::endl;
          system_matrix.resize(map_index, map_index, false);
          system_matrix.clear();
          system_matrix.resize(map_index, map_index, false);
        }
        
        if (map_index > rhs.size())
        {
          std::cout << "Resizing load vector..." << std::endl;
          rhs.resize(map_index, false);
          rhs.clear();
          rhs.resize(map_index, false);
        }
        
        
        //        
        //Poisson equation:  laplace psi = 1
        //

        //VertexContainer vertices = viennagrid::ncells<0>(domain);
        for (VertexIterator vit = vertices.begin();
             vit != vertices.end();
             ++vit)
        {
          //std::cout << "* poisson_solver::assemble(): Iterating over vertex: " << std::endl;
          //vit->print();
          long row_index = viennadata::access<mapping_key, long>(mapping_key(0))(*vit);
          
          if (row_index < 0)
            continue;
          
          EdgeOnVertexContainer edges = viennagrid::ncells<1>(*vit, domain);
          for (EdgeOnVertexIterator eovit = edges.begin();
               eovit != edges.end();
               ++eovit)
          {
            //std::cout << "* poisson_solver::assemble(): Iterating over edge: " << std::endl;
            //eovit->print_short();
            
            VertexOnEdgeContainer vertices_on_edge = viennagrid::ncells<0>(*eovit);
            VertexOnEdgeIterator voeit = vertices_on_edge.begin();
            
            if ( &(*voeit) == &(*vit))  //one of the two vertices of the edge is different from *vit
              ++voeit;
            
            //std::cout << "* poisson_solver::assemble(): Other vertex: " <<std::endl;
            //voeit->print();
            
            //std::cout << "* poisson_solver::assemble(): Getting data... " <<std::endl;
            long col_index        = viennadata::access<mapping_key, long>(mapping_key(0))(*voeit);
            double edge_len       = viennagrid::volume(*eovit);
            double interface_area = viennadata::access<edge_interface_area_key, double>()(*eovit);
            
            //std::cout << "edge_len: " << edge_len << std::endl;
            //std::cout << "interface_area: " << interface_area << std::endl;

            //std::cout << "* poisson_solver::assemble(): Assembling... " <<std::endl;
            //std::cout << "row_index: " << row_index << std::endl;
            //std::cout << "col_index: " << col_index << std::endl;
            if (col_index < 0)
            {
              //*voeit is a Dirichlet boundary: do nothing for the moment
            }
            else
            { 
              //std::cout << "Writing " << (-interface_area / edge_len) << " at " << row_index << ", " << col_index << std::endl;
              system_matrix(row_index, col_index) = - interface_area / edge_len;
            }
            
            //std::cout << "* poisson_solver::assemble(): Writing matrix diagonal... " <<std::endl;
            //std::cout << "Writing " << (+interface_area / edge_len) << " at " << row_index << ", " << row_index << std::endl;
            system_matrix(row_index, row_index) += interface_area / edge_len;
            

            double box_volume_part = viennadata::access<box_volume_key, double>()(*eovit) / 2.0;  //Note: volume stored on edges consists of volumes of both adjacent boxes.
            
            //write rhs vector entry: 
            rhs(row_index) +=  box_volume_part;  //flux term equals '1' here
            
          } //for edges
           
        } //for vertices   
      } //operator()
      
    }; //poisson_assembler

    
} //namespace viennafvm

#endif
