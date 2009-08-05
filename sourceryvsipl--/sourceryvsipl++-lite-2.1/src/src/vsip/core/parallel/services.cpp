/* Copyright (c) 2005 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/parallel.cpp
    @author  Jules Bergmann
    @date    2005-03-25
    @brief   VSIPL++ Library: Declarations for parallel services.

*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/parallel/services.hpp>
#include <vsip/core/vector.hpp>



/***********************************************************************
  Declarations
***********************************************************************/

vsip::impl::Communicator vsip::impl::Par_service::default_communicator_;



/***********************************************************************
  Definitions
***********************************************************************/

namespace vsip
{

namespace impl
{

namespace pas
{

/// Counter to generate a unique tag for global PAS pbuffer allocations.
long             global_tag = 1;

} // namespace vspi::impl::pas
} // namespace vspi::impl

/// Return the number of processors in the data parallel clique.

length_type
num_processors()
  VSIP_NOTHROW
{
  return impl::default_communicator().size();
}



/// Return the set of processors in the data parallel clique.

const_Vector<processor_type>
processor_set()
{
  static Dense<1, processor_type>* pset_block_ = NULL;

  if (pset_block_ == NULL)
  {
    impl::Communicator::pvec_type const& 
      pvec = impl::default_communicator().pvec(); 

    pset_block_ = new Dense<1, processor_type>(Domain<1>(pvec.size()));
    for (index_type i=0; i<pvec.size(); ++i)
      pset_block_->put(i, pvec[i]);
  }

  return Vector<processor_type>(*pset_block_);
}



/// Return the local processor.

processor_type
local_processor()
  VSIP_NOTHROW
{
  return impl::default_communicator().rank();
}


} // namespace vsip
