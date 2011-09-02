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



#ifndef VIENNAFVM_ACC_HPP
#define VIENNAFVM_ACC_HPP

// *** system includes
// *** local includes
// *** vienna includes
// *** boost includes
#include <boost/fusion/include/at_key.hpp>
#include <boost/fusion/include/value_at_key.hpp>


namespace viennafvm {

// -----------------------------------------------------------------------------
template<typename AccessTag, typename MapSequenceT>
inline typename boost::fusion::result_of::value_at_key<MapSequenceT, AccessTag>::type const&
acc(MapSequenceT const& map)
{
   return boost::fusion::at_key<AccessTag>(map);
}

template<typename AccessTag, typename MapSequenceT>
inline typename boost::fusion::result_of::value_at_key<MapSequenceT, AccessTag>::type &
acc(MapSequenceT & map)
{
   return boost::fusion::at_key<AccessTag>(map);
}
// -----------------------------------------------------------------------------

} // end namespace viennafvm

#endif
