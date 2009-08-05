/* Copyright (c) 2006 by CodeSourcery.  All rights reserved. */

/** @file    vsip/opt/choose_dist_block.hpp
    @author  Jules Bergmann
    @date    2006-08-29
    @brief   VSIPL++ Library: Choose distributed block implementation.

*/

#ifndef VSIP_OPT_CHOOSE_DIST_BLOCK_HPP
#define VSIP_OPT_CHOOSE_DIST_BLOCK_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/config.hpp>
#include <vsip/core/metaprogramming.hpp>
#include <vsip/core/block_traits.hpp>
#include <vsip/core/parallel/map_traits.hpp>
#include <vsip/core/dense_fwd.hpp>



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

/// Forward Declaration
template <typename Block,
	  typename Map>
class Distributed_block;


/// Forward Declaration
template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       MapT>
class Pas_block;


template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       MapT>
struct Choose_dist_block
{
#if VSIP_IMPL_PAR_SERVICE == 1
  typedef Distributed_block<Dense<Dim, T, OrderT, Local_map>, MapT> type;
#elif VSIP_IMPL_PAR_SERVICE == 2
  typedef Pas_block<Dim, T, OrderT, MapT> type;
#else
  typedef Distributed_block<Dense<Dim, T, OrderT, Local_map>, MapT> type;
#endif
};
	  

} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_OPT_CHOOSE_DIST_BLOCK_HPP
