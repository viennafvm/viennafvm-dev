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


#ifndef VIENNAFVM_IO_VTKWRITER_HPP
#define VIENNAFVM_IO_VTKWRITER_HPP

// include necessary system headers
#include <iostream>

#include "viennafvm/forwards.h"

// ViennaGrid includes:
#include "viennagrid/domain.hpp"
#include "viennagrid/iterators.hpp"
#include "viennagrid/io/vtk_writer.hpp"

// ViennaData includes:
#include "viennadata/api.hpp"

namespace viennafvm
{
  namespace io
  {

    template <typename VectorType,
              typename DomainType>
    void write_solution_to_VTK_file(VectorType const & result,
                                    std::string filename,
                                    DomainType const & domain,
                                    long id)
    {
      typedef typename DomainType::config_type                                              ConfigType;
      typedef typename ConfigType::cell_tag                                                 CellTag;

      typedef typename viennagrid::result_of::ncell<ConfigType, CellTag::dim>::type         CellType;

      typedef typename viennagrid::result_of::const_ncell_range<DomainType, CellTag::dim>::type   CellContainer;
      typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;
      
      typedef viennafvm::mapping_key          MappingKeyType;
      typedef viennafvm::boundary_key         BoundaryKeyType;
      
      MappingKeyType map_key(id);
      BoundaryKeyType bnd_key(id);
      
      
      std::cout << "* write_solution_to_VTK_file(): Writing result on mesh for later export" << std::endl;
      CellContainer cells = viennagrid::ncells<CellTag::dim>(domain);
      for (CellIterator cit = cells.begin();
                        cit != cells.end();
                      ++cit)
      {
        long cur_index = viennadata::access<MappingKeyType, long>(map_key)(*cit);
        if (cur_index > -1)
          viennadata::access<std::string, double>("vtk_data")(*cit) = result[cur_index];
        else //use Dirichlet boundary data:
          viennadata::access<std::string, double>("vtk_data")(*cit) = viennadata::access<BoundaryKeyType, double>(bnd_key)(*cit);
      }

      std::cout << "* write_solution_to_VTK_file(): Writing data to '"
                << filename
                << "' (can be viewed with e.g. Paraview)" << std::endl;

      viennagrid::io::vtk_writer<DomainType> my_vtk_writer;
      viennagrid::io::add_scalar_data_on_cells<std::string, double>(my_vtk_writer, "vtk_data", "fvm_result");
      my_vtk_writer(domain, filename);  
    }


  }
}
#endif

