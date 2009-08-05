/* Copyright (c) 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/

/** @file    apps/utils/lsmat.cpp
    @author  Jules Bergmann
    @date    2007-07-13
    @brief   VSIPL++ Library: List views contained in .mat file.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <iostream>
#include <fstream>

#include <vsip/initfin.hpp>
#include <vsip/vector.hpp>
#include <vsip/matrix.hpp>
#include <vsip/tensor.hpp>
#include <vsip/selgen.hpp>
#include <vsip/map.hpp>

#include <vsip_csl/matlab_file.hpp>
#include <vsip_csl/output.hpp>

using namespace std;
using namespace vsip;
using namespace vsip_csl;



/***********************************************************************
  Definitions
***********************************************************************/

template <dimension_type Dim,
	  typename       T>
void
load_view(
  Matlab_file&           mf,
  Matlab_file::iterator& iter,
  Domain<Dim> const&     dom)
{
  typedef Dense<Dim, T> block_type;
  typedef typename impl::View_of_dim<Dim, T, block_type>::type view_type;

  block_type block(dom);
  view_type  view(block);

  mf.read_view(view, iter);

  cout << view;
}



void
lsmat(char const* fname)
{
  // Create Matlab_file object for 'sample.mat' file.
  Matlab_file mf(fname);
  Matlab_file::iterator cur = mf.begin();
  Matlab_file::iterator end = mf.end();

  Matlab_bin_hdr hdr = mf.header();

  cout << "Matlab file: " << fname << endl;
  cout << "  descr  : " << hdr.description << endl;
  cout << "  version: " << hdr.version << endl;
  cout << "  endian : " << hdr.endian 
       << " (swap: " 
       << (hdr.endian == ('I'<<8 | 'M') ? "yes" :
	   hdr.endian == ('M'<<8 | 'I') ? "no"  : "*unknown*")
       << ")" << endl;


  // Iterate through views in file.
  while (cur != end)
  {
    Matlab_view_header* hdr = *cur;

    cout << "view: " << hdr->array_name << endl;

    // Dump array_name out by byte
    // for (int i=0; hdr->array_name[i] != 0 && i<128; ++i)
    //   cout << "  [" << i << "]: " << (int)hdr->array_name[i] << endl;

    cout << "  dim       : " << hdr->num_dims;
    if (hdr->num_dims > 0)
    {
      char* sep = " (";
      for (index_type i=0; i<hdr->num_dims; ++i)
      {
	cout << sep << hdr->dims[i];
	sep = ", ";
      }
      cout << ")";
    }
    cout << endl;

    cout << "  is_complex: " << (hdr->is_complex ? "true" : "false") << endl;
    cout << "  class_type: " << vsip_csl::matlab::class_type(hdr->class_type);
    cout << " (" << (int)hdr->class_type << ")" << endl;

    if (hdr->class_type == vsip_csl::matlab::mxDOUBLE_CLASS)
    {
      if (hdr->num_dims == 2 && hdr->dims[0] == 1)
	load_view<1, double>(mf, cur, Domain<1>(hdr->dims[1]));
      else if (hdr->num_dims == 2 && hdr->dims[1] == 1)
	load_view<1, double>(mf, cur, Domain<1>(hdr->dims[0]));
      else if (hdr->num_dims == 2)
	load_view<2, double>(mf, cur, Domain<2>(hdr->dims[0], hdr->dims[1]));
//      else if (hdr->num_dims == 3)
//	load_view<3, double>(mf, cur, Domain<3>(hdr->dims[0], hdr->dims[1],
//						hdr->dims[2]));
    }
    else if (hdr->class_type == vsip_csl::matlab::mxSINGLE_CLASS)
    {
      if (hdr->num_dims == 2 && hdr->dims[0] == 1)
	load_view<1, float>(mf, cur, Domain<1>(hdr->dims[1]));
      else if (hdr->num_dims == 2 && hdr->dims[1] == 1)
	load_view<1, float>(mf, cur, Domain<1>(hdr->dims[0]));
      else if (hdr->num_dims == 2)
	load_view<2, float>(mf, cur, Domain<2>(hdr->dims[0], hdr->dims[1]));
//      else if (hdr->num_dims == 3)
//	load_view<3, float>(mf, cur, Domain<3>(hdr->dims[0], hdr->dims[1],
//						hdr->dims[2]));
    }

    ++cur; // Move to next view stored in the file.
  }
}


int
main(int argc, char** argv)
{
  (void)argc;
  char* fname = argv[1];

  lsmat(fname);
} 

