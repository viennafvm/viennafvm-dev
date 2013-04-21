#ifndef VIENNAFVM_IO_VTKWRITER_HPP
#define VIENNAFVM_IO_VTKWRITER_HPP

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
                                    std::vector<long> id_vector)
    {
      typedef typename DomainType::config_type                                              ConfigType;
      typedef typename ConfigType::cell_tag                                                 CellTag;

      typedef typename viennagrid::result_of::ncell<ConfigType, CellTag::dim>::type         CellType;

      typedef typename viennagrid::result_of::const_ncell_range<DomainType, CellTag::dim>::type   CellContainer;
      typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

      typedef viennafvm::mapping_key          MappingKeyType;
      typedef viennafvm::boundary_key         BoundaryKeyType;

      std::cout << "* write_solution_to_VTK_file(): Writing result on mesh for later export" << std::endl;
      viennagrid::io::vtk_writer<DomainType> my_vtk_writer;

      for (std::size_t i=0; i<id_vector.size(); ++i)
      {
        long id = id_vector[i];
        MappingKeyType map_key(id);
        BoundaryKeyType bnd_key(id);

        std::stringstream ss;
        ss << "fvm_result" << id;
        std::string result_string = ss.str(); // also used for passing staff over to VTK

        CellContainer cells = viennagrid::ncells<CellTag::dim>(domain);
        for (CellIterator cit = cells.begin();
                          cit != cells.end();
                        ++cit)
        {
          long cur_index = viennadata::access<MappingKeyType, long>(map_key)(*cit);
          if (cur_index > -1)
            viennadata::access<std::string, double>(result_string)(*cit) = result[cur_index];
          else //use Dirichlet boundary data:
            viennadata::access<std::string, double>(result_string)(*cit) = viennadata::access<BoundaryKeyType, double>(bnd_key)(*cit);
        }

        viennagrid::io::add_scalar_data_on_cells<std::string, double>(my_vtk_writer, result_string, result_string);
      }

      std::cout << "* write_solution_to_VTK_file(): Writing data to '"
                << filename
                << "' (can be viewed with e.g. Paraview)" << std::endl;

      my_vtk_writer(domain, filename);
    }

    template <typename VectorType,
              typename DomainType>
    void write_solution_to_VTK_file(VectorType const & result,
                                    std::string filename,
                                    DomainType const & domain,
                                    long id)
    {
      std::vector<long> id_vector(1);
      id_vector[0] = id;

      write_solution_to_VTK_file(result, filename, domain, id_vector);
    }
  }
}
#endif

