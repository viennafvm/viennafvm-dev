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

    template <typename QuantityContainerT,
              typename MeshT,
              typename SegmentationType>
    void write_solution_to_VTK_file(QuantityContainerT const & quantities,
                                    std::string filename,
                                    MeshT const & mesh,
                                    SegmentationType const & segmentation)
    {
      typedef typename viennagrid::result_of::cell_tag<MeshT>::type CellTag;

      typedef typename viennagrid::result_of::element<MeshT, CellTag>::type         CellType;

      typedef typename viennagrid::result_of::const_element_range<MeshT, CellTag>::type   CellContainer;
      typedef typename viennagrid::result_of::iterator<CellContainer>::type                       CellIterator;

    #ifdef VIENNAFVM_VERBOSE
      std::cout << "* write_solution_to_VTK_file(): Writing result on mesh for later export" << std::endl;
    #endif
      viennagrid::io::vtk_writer<MeshT> my_vtk_writer;

      std::map< std::string, std::deque<double> > output_values;

      for (std::size_t i=0; i<quantities.size(); ++i)
      {
        typename viennagrid::result_of::accessor< std::deque<double>, CellType >::type output_value_accessor( output_values[quantities.at(i).get_name()] );

        CellContainer cells(mesh);
        for (CellIterator cit = cells.begin();
                          cit != cells.end();
                        ++cit)
        {
          output_value_accessor(*cit) = quantities.at(i).get_value(*cit);
        }

        my_vtk_writer.add_scalar_data_on_cells( output_value_accessor, quantities.at(i).get_name() );
      }

    #ifdef VIENNAFVM_VERBOSE
      std::cout << "* write_solution_to_VTK_file(): Writing data to '"
                << filename
                << "' (can be viewed with e.g. Paraview)" << std::endl;
    #endif
      my_vtk_writer(mesh, segmentation, filename);
    }

  }
}
#endif

