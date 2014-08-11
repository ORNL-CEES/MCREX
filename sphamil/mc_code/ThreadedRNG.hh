//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ThreadedRNG.hh
 * \author Steven Hamilton
 * \brief  Thread local generation of random numbers
 */
//---------------------------------------------------------------------------//

#ifndef ThreadedRNG_hh
#define ThreadedRNG_hh

#include <vector>
#include <random>

namespace mcrex
{

//
//---------------------------------------------------------------------------//
/*!
 * \class ThreadedRNG
 * \brief Thread-safe random number generation.
 *
 *  Random number generation with c++11 involves creation of engine and
 *  distribution objects.  Creating random numbers on multiple threads
 *  from the same distribution/engine objects creates a performance
 *  bottleneck due to locks imposed in the implementation.  We address
 *  this by creating a separate objects for each thread.  Requests for
 *  random numbers are routed to the appropriate objects based on the thread
 *  id provided.
 */
//---------------------------------------------------------------------------//
class ThreadedRNG
{
  public:

    /*!
     * \brief Constructor.
     *
     * \param num_threads Number of threads that will be calling generator
     */
    explicit inline ThreadedRNG(const int num_threads);

    //! \brief Get next random number for specified thread.
    inline double getRandom(const int thread) const;

  private:

    mutable std::vector<std::mt19937_64> d_engines;

};

} // namespace mcrex

#include "ThreadedRNG.i.hh"

#endif // ThreadedRNG_hh

