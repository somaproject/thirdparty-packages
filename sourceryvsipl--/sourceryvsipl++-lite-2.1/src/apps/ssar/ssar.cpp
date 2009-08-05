/* Copyright (c) 200, 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    ssar.cpp
    @author  Don McCoy
    @date    2006-10-26
    @brief   VSIPL++ implementation of HPCS Challenge Benchmarks 
               Scalable Synthetic Compact Applications - 
             SSCA #3: Sensor Processing and Knowledge Formation
*/

#include <iostream>
#include <fstream>
#include <cerrno>

#include <vsip/initfin.hpp>
#include <vsip/math.hpp>
#include <vsip/signal.hpp>
#include <vsip_csl/output.hpp>
#include <vsip_csl/memory_pool.hpp>

using namespace vsip_csl;
using namespace vsip;
using namespace std;

// Files required to be in the data directory (must be included 
// before kernel1.hpp):
char const* SAR_DIMENSIONS =                         "dims.txt";
char const* RAW_SAR_DATA =                           "sar.view";
char const* FAST_TIME_FILTER =                       "ftfilt.view";
char const* SLOW_TIME_WAVENUMBER =                   "k.view";
char const* SLOW_TIME_COMPRESSED_APERTURE_POSITION = "uc.view";
char const* SLOW_TIME_APERTURE_POSITION =            "u.view";
char const* SLOW_TIME_SPATIAL_FREQUENCY =            "ku.view";

#include "kernel1.hpp"

struct
ssar_options
{
  scalar_f scale;
  length_type n;
  length_type mc;
  length_type m;
};

void
process_ssar_options(int argc, char** argv, ssar_options& options);

// Options set in process_ssar_options

int loop = 1;  // Number of process_image iterations to perform (default 1)

#if _BIG_ENDIAN
bool swap_bytes = true;   // Whether or not to swap bytes during file I/O
#else
bool swap_bytes = false;
#endif

int
main(int argc, char** argv)
{
  using vsip_csl::load_view_as;
  using vsip_csl::save_view_as;

  vsip::vsipl init(argc, argv);

  Local_map huge_map;

  ssar_options opt;
  process_ssar_options(argc, argv, opt);

  typedef SSAR_BASE_TYPE T;

#if VSIP_IMPL_ENABLE_HUGE_PAGE_POOL && VSIP_IMPL_CBE_SDK
  set_pool(huge_map, new vsip_csl::Huge_page_pool("/huge/benchmark.bin", 20));
#endif

  // Setup for Stage 1, Kernel 1 
  vsip::impl::profile::Acc_timer t0;
  t0.start();
  Kernel1<T> k1(opt.scale, opt.n, opt.mc, opt.m, swap_bytes, huge_map); 
  t0.stop();
  cout << "setup:   " << t0.delta() << " (s)" << endl;

  // Retrieve the raw radar image data from disk.  This Data I/O 
  // component is currently untimed.
  Kernel1<T>::complex_matrix_type s_raw(opt.n, opt.mc, huge_map);
  load_view_as<complex<float>, 
    Kernel1<T>::complex_matrix_type>(RAW_SAR_DATA, s_raw, swap_bytes);

  // Resolve the image.  This Computation component is timed.
  Kernel1<T>::real_matrix_type 
    image(k1.output_size(0), k1.output_size(1), huge_map);

  vsip::impl::profile::Acc_timer t1;
  vsip::Vector<double> process_time(loop);
  for (int l=0; l<loop; ++l)
  {
    t1.start();
    k1.process_image(s_raw, image);
    t1.stop();
    process_time.put(l, t1.delta());
  }

  // Display statistics
  if (loop > 0)
  {
    Index<1> idx;
    double mean = vsip::meanval(process_time);
    cout << "loops:   " << loop << endl;
    cout << "mean:    " << mean << endl;
    cout << "min:     " << vsip::minval(process_time, idx) << endl;
    cout << "max:     " << vsip::maxval(process_time, idx) << endl;
    cout << "std-dev: " << sqrt(vsip::meansqval(process_time - mean)) << endl;
  }

  // Store the image on disk for later processing (not timed).
  save_view_as<float>("image.view", image, swap_bytes); 
}


void
process_ssar_options(int argc, char** argv, ssar_options& options)
{
  char* dir = 0;

  for (int i=1; i<argc; ++i)
  {
    if (!strcmp(argv[i], "-loop")) loop = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-swap")) swap_bytes = true;
    else if (!strcmp(argv[i], "-noswap")) swap_bytes = false;
    else if (dir == 0) dir = argv[i];
    else
    {
      cerr << "Unknown arg: " << argv[i] << endl;
      cerr << "Usage: " << argv[0] << " [-loop] <data dir> [-swap|-noswap]" << endl;
      exit(-1);
    }
  }

  if (dir == 0)
  {
      cerr << "No dir given" << endl;
      cerr << "Usage: " << argv[0] << " [-loop] <data dir> [-swap|-noswap]" << endl;
      exit(-1);
  }

  if (chdir(dir) < 0)
  {
    perror(dir);
    exit(-1);
  }

  ifstream fp_in(SAR_DIMENSIONS);
  if (fp_in.fail())
  {
    perror(SAR_DIMENSIONS);
    exit(-1);
  }

  fp_in >> options.scale;
  fp_in >> options.n;
  fp_in >> options.mc;
  fp_in >> options.m;

  if (fp_in.fail())
  {
    cerr << "Error reading dimension data" << endl;
    exit(-1);
  }

  fp_in.close();
}

