/* Copyright (c) 2007 by CodeSourcery.  All rights reserved.
*/
/** @file    examples/extdata.cpp
    @author  Don McCoy
    @date    2007-07-23
    @brief   VSIPL++ Library: Simple example for external data access.
*/

#include <iostream>
#include <vsip/initfin.hpp>
#include <vsip/matrix.hpp>
#include <vsip/opt/extdata.hpp>


int 
main(int argc, char **argv)
{
  vsip::vsipl init(argc, argv);

  const vsip::length_type M = 128;
  const vsip::length_type N = 64;
  
  typedef vsip::scalar_f T;
  typedef vsip::Dense<2, T, vsip::row2_type> block_type;
  typedef vsip::Matrix<T, block_type> view_type;

  view_type grid(M, N, 1.234f);

  std::cout << "first and last values (before): " << grid.get(0,0) 
            << ", " << grid.get(M-1,N-1) << std::endl;

  // Control the scope of the Ext_data object.  Synchronization is forced
  // when the object is destroyed in accordance with the SYNC_OUT parameter.
  {
    vsip::impl::Ext_data<block_type> ext(grid.block(), vsip::impl::SYNC_OUT);

    T* ptr = ext.data();
    for (vsip::length_type i = 0; i < M; ++i)
      for (vsip::length_type j = 0; j < N; ++j)
        *(ptr + i*N + j) = i*N + j;
    std::cout << "ptr = " << ptr << std::endl;
  }
   
  std::cout << "first and last values (after):  " << grid.get(0,0) 
            << ", " << grid.get(M-1,N-1) << std::endl;
  std::cout << "expected values:  0, " << static_cast<float>(M*N) << std::endl;

  return 0;
}
