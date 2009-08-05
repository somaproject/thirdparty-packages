/* Copyright (c) 2009 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vadd.cpp
    @author  Mike LeBlanc
    @date    2009-04-07
    @brief   ...
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

void
process_vadd_options(int argc, char** argv);

int loop = 1;
int nn = 1024;

int
main(int argc, char** argv)
{

  vsip::vsipl init(argc, argv);

  process_vadd_options(argc, argv);

  impl::Communicator comm = impl::default_communicator();

  pid_t pid = getpid();

  int r = comm.rank();

#if 0
  cout << "rank: "   << r
       << "  size: " << comm.size()
       << "  pid: "  << pid 
       << endl;

  // Enable this section for easier debugging.
  // Stop each process, allow debugger to be attached.
  if (r == 0) fgetc(stdin);
#endif
  comm.barrier();

  const length_type N = nn;

//cout << "[" << r << "] start, N=" << N << endl;
  //
  // Declare a map type to control the partitioning of vectors
  // among several processors.  We don't know the exact number
  // until run-time.
  //
  typedef Map<Block_dist> vector_map_type;
  //
  // Declare an instance of that map type and initialize it
  // with the number of processors.
  //
  vector_map_type map = vector_map_type(num_processors());
  //
  // Declare a block type based on our map type.
  //
  typedef Dense<1, float, row2_type, vector_map_type> vector_block_type;
  //
  // Declare a vector view based on our block type.
  //
  typedef Vector<float, vector_block_type> vector_type;
  //
  // Finally, declare vectors of our view type, using the
  // run-time knowledge encapsulated in our map instance.
  //
  vector_type A(N,map);
  vector_type B(N,map);
  vector_type C(N,map);
  //
  // Initialize the vectors.  Since our algorithm is A=B+C, the
  // result in A(i) is 3*i.
  //
  A = 0.0f;
  B = vsip::ramp(0.0f,1.0f,N);
  C = vsip::ramp(0.0f,2.0f,N);
  //
  // Create "local" views of the full vectors.  Each one
  // will see only the partition distributed to the
  // processor on which it is created.
  //
  vector_type::local_type A_local = A.local();
  vector_type::local_type B_local = B.local();
  vector_type::local_type C_local = C.local();
  //
  // Note the number of elements in the distributed
  // partitions.
  //
  length_type asize = A_local.size(0);
  
  using impl::profile::Scope;
  using impl::profile::user;
  // 
  // Run the algorithm over and over and collect timings.
  // 
  vsip::impl::profile::Acc_timer t1;
  vsip::Vector<double> process_time(loop);
  comm.barrier();
  for (int l=0; l<loop; ++l)
  {
    t1.start();
    {
      Scope<user> scope("algorithm", asize*sizeof(float));
      //
      // Run our simple algorithm.  Add partitions from vectors
      // B and C and write the result into C's partition.
      //
      A_local = B_local + C_local;
    }
    t1.stop();
    process_time.put(l, t1.delta());
  }
  comm.barrier();
  if (r==0 && loop>0)
    printf("Mean time: %12.10f\n", vsip::meanval(process_time));
  //
  // Declare a globally mapped vector where each processor
  // can store its results.
  //
  Vector<
    float,
    Dense<
      1,
      float,
      row2_type,
      Global_map<1>
      >
    >
    R(N);
  //
  // Copy the distributed vector to the global vector.
  //
  R = A;
  //
  // Wait for all processors to arrive here.
  //
  comm.barrier();
//cout << "[" << r << "] finished" << endl;
  //
  // Processor 0 checks results.  It looks into the non-local
  // view to see all elements of A.
  //
  if (r == 0)
  {
//  cout << "[" << r << "] " << "checking ..." << endl;
    for (index_type i=0; i<N; ++i)
    {
      float x = 3.0f*i;
      if (R(i) != x)
      {
        cout << "[" << r << "] " << "A(" << i << ") should be " << x << " but is " << R(i) << endl;
        break;
      }
    }
  }
  comm.barrier();

#if 0
  // Display statistics
  if (loop > 0)
  {
    Index<1> idx;
    double mean = vsip::meanval(process_time);
    cout << "[" << r << "] " << "loops:   " << loop << endl;
    cout << "[" << r << "] " << "mean:    " << mean << endl;
    cout << "[" << r << "] " << "min:     " << vsip::minval(process_time, idx) << endl;
    cout << "[" << r << "] " << "max:     " << vsip::maxval(process_time, idx) << endl;
    cout << "[" << r << "] " << "std-dev: " << sqrt(vsip::meansqval(process_time - mean)) << endl;
  }
#endif

  return 0;
}


void
process_vadd_options(int argc, char** argv)
{

  for (int i=1; i<argc; ++i)
  {
    if (!strcmp(argv[i], "-loop")) loop = atoi(argv[++i]);
    else
    if (!strcmp(argv[i], "-N")) nn = atoi(argv[++i]);
    else
    {
      cerr << "Unknown arg: " << argv[i] << endl;
      cerr << "Usage: " << argv[0] << " [-loop n] [-N n] " << endl;
      exit(-1);
    }
  }

}

