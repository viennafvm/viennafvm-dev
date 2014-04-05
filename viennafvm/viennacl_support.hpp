#ifndef VIENNAFVM_VIENNACL_SUPPORT_HPP
#define VIENNAFVM_VIENNACL_SUPPORT_HPP

/* =======================================================================
   Copyright (c) 2011, Institute for Microelectronics, TU Wien
   http://www.iue.tuwien.ac.at
                             -----------------
           ViennaFVM - The Vienna Finite Volume Method Library
                             -----------------

   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at
               Josef Weinbub                   weinbub@iue.tuwien.ac.at
               (add your name here)

   license:    see file LICENSE in the ViennaFVM base directory
======================================================================= */

#ifdef VIENNACL_WITH_OPENCL

#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/platform.hpp"

namespace viennafvm {

inline void print_current_device(std::ostream& stream = std::cout)
{
  stream << "  -----------------------------------------" << std::endl;
  stream << "  Current OpenCL Device:" << std::endl;
  stream << viennacl::ocl::current_device().info();
  stream << "  -----------------------------------------" << std::endl;
}

inline void print_platform_devices(std::ostream& stream = std::cout)
{
   //
   //  retrieve the devices
   //
   typedef std::vector< viennacl::ocl::platform > platforms_type;
   platforms_type platforms = viennacl::ocl::get_platforms();

   bool is_first_element = true;
   for (platforms_type::iterator platform_iter  = platforms.begin();
                                 platform_iter != platforms.end();
                               ++platform_iter)
   {
    typedef std::vector<viennacl::ocl::device> devices_type;
    devices_type devices = platform_iter->devices(CL_DEVICE_TYPE_ALL);

    //
    // print some platform info
    //
    stream << "# =========================================" << std::endl;
    stream << "#         Platform Information             " << std::endl;
    stream << "# =========================================" << std::endl;
    
    stream << "#" << std::endl;
    stream << "# Vendor and version: " << platform_iter->info() << std::endl;
    stream << "#" << std::endl;
    
    if (is_first_element)
    {
      stream << "# ViennaCL uses this OpenCL platform by default." << std::endl;
      is_first_element = false;
    }
    
    
    //
    //  traverse the devices and print the information
    //
    stream << "# " << std::endl;
    stream << "# Available Devices: " << std::endl;
    stream << "# " << std::endl;
    for(devices_type::iterator iter = devices.begin(); iter != devices.end(); iter++)
    {
        stream << std::endl;

        if(*iter == viennacl::ocl::current_device())
          stream << "  (Current) -------------------------------" << std::endl;
        else
          stream << "  -----------------------------------------" << std::endl;
          stream << iter->info();
        stream << "  -----------------------------------------" << std::endl;
    }
    stream << std::endl;
    stream << "###########################################" << std::endl;
    stream << std::endl;
   }
   
}

} // viennafvm

#endif
#endif
