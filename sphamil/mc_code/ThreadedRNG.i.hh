//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ThreadedRNG.i.hh
 * \author Steven Hamilton
 * \brief  Thread local generation of random numbers
 */
//---------------------------------------------------------------------------//

#ifndef MCREX_ThreadedRNG_i_hh
#define MCREX_ThreadedRNG_i_hh

#include "ThreadedRNG.hh"

#include <vector>
#include <random>
#include <limits>

#include "Teuchos_Assert.hpp"
#include "MCREX_Typedefs.hh"

namespace mcrex
{

// Constructor, initialize generators
ThreadedRNG::ThreadedRNG(const int num_threads)
{
    // Create unique RNG engine and distribution
    // for each thread.  Seed engines separately.
    d_engines.resize(num_threads);
    std::random_device rd;
    for( int i=0; i<num_threads; ++i )
    {
        d_engines[i].seed(rd());
    }
}

// Get random number on calling thread
double ThreadedRNG::getRandom(const int thread) const
{
    TEUCHOS_ASSERT( static_cast<size_t>(thread) < d_engines.size() );

    // We could create a uniform_real_distribution object
    // and create a random number using that.  However,
    // c++11 RNGs natively produce their output on the
    // interval [0,1).  We can instead directly call
    // generate_canonical to avoid an unnecessary transformation.
    // This has a slight performance benefit.
    const int digits = std::numeric_limits<SCALAR>::digits;
    return std::generate_canonical<SCALAR,digits>(d_engines[thread]);
}

} // namespace mcrex
#endif // MCREX_ThreadedRNG_i_hh

