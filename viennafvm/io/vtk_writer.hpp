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
#include "viennagrid/io/vtk_writer.hpp"

// ViennaData includes:
#include "viennadata/api.hpp"

namespace viennafvm
{
  namespace io
  {

    template <typename StorageType,
              typename VectorType,
              typename DomainType,
              typename SegmentationType>
    void write_solution_to_VTK_file(StorageType const & storage,
                                    VectorType const & result,
                                    std::string filename,
                                    DomainType const & domain,
                                    SegmentationType const & segmentation,
                                    std::vector<long> id_vector)
    {
//       typedef typename DomainType::config_type                                              ConfigType;
//       typedef typename ConfigType::cell_tag                                                 CellTag;

      typedef typename viennagrid::result_of::cell_tag<DomainType>::type CellTag;

      typedef typename viennagrid::result_of::element<DomainType, CellTag>::type         CellType;

      typedef typename viennagrid::result_of::const_element_range<DomainType, CellTag>::type   CellContainer;
      typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

      typedef viennafvm::mapping_key          MappingKeyType;
      typedef viennafvm::boundary_key         BoundaryKeyType;

      std::cout << "* write_solution_to_VTK_file(): Writing result on mesh for later export" << std::endl;
      viennagrid::io::vtk_writer<DomainType> my_vtk_writer;



      std::map< std::string, std::deque<double> > output_values;

      for (std::size_t i=0; i<id_vector.size(); ++i)
      {
        long id = id_vector[i];
        MappingKeyType map_key(id);
        BoundaryKeyType bnd_key(id);

        typename viennadata::result_of::accessor<const StorageType, viennafvm::mapping_key, long, CellType>::type cell_mapping_accessor =
          viennadata::accessor<viennafvm::mapping_key, long, CellType>(storage, map_key);
        
        typename viennadata::result_of::accessor<const StorageType, BoundaryKeyType, double, CellType>::type boudnary_accessor =
          viennadata::accessor<BoundaryKeyType, double, CellType>(storage, bnd_key);

        std::stringstream ss;
        ss << "fvm_result" << id;
        std::string result_string = ss.str(); // also used for passing staff over to VTK

        typename viennagrid::result_of::accessor< std::deque<double>, CellType >::type output_value_accessor( output_values[result_string] );

        CellContainer cells = viennagrid::elements(domain);
        for (CellIterator cit = cells.begin();
                          cit != cells.end();
                        ++cit)
        {
          long cur_index = cell_mapping_accessor(*cit);
          if (cur_index > -1)
            output_value_accessor(*cit) = result[cur_index];
          else //use Dirichlet boundary data:
            output_value_accessor(*cit) = boudnary_accessor(*cit);
        }

        my_vtk_writer.add_scalar_data_on_cells( output_value_accessor, result_string );
      }

      std::cout << "* write_solution_to_VTK_file(): Writing data to '"
                << filename
                << "' (can be viewed with e.g. Paraview)" << std::endl;

      my_vtk_writer(domain, segmentation, filename);
    }

    template <typename StorageType,
              typename VectorType,
              typename DomainType,
              typename SegmentationType>
    void write_solution_to_VTK_file(StorageType const & storage,
                                    VectorType const & result,
                                    std::string filename,
                                    DomainType const & domain,
                                    SegmentationType const & segmentation,
                                    long id)
    {
      std::vector<long> id_vector(1);
      id_vector[0] = id;

      write_solution_to_VTK_file(storage, result, filename, domain, segmentation, id_vector);
    }
  }
}
#endif

