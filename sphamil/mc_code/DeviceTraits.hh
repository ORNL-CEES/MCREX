//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DeviceTraits.hh
 * \author Steven Hamilton
 * \brief  Templated interface for Kokkos devices.
 */
//---------------------------------------------------------------------------//

#ifndef MCREX_DeviceTraits_hh
#define MCREX_DeviceTraits_hh

#include <omp.h>

#include "Kokkos_hwloc.hpp"

// "New" Kokkos nodes
#include "Kokkos_Serial.hpp"
#include "Kokkos_OpenMP.hpp"
#include "Kokkos_Threads.hpp"

namespace mcrex
{

//---------------------------------------------------------------------------//
/*!
 * \class DeviceTraits
 * \brief Templated interface for initializing/finalizing Kokkos devices.
 *
 * This is a generic implementation for a single thread with trivial
 * initialization/finalization which should be suitable for Kokkos::Serial.
 */
//---------------------------------------------------------------------------//
template <class DeviceType>
class DeviceTraits
{
  public:

    //! \brief Initialize device on specified number of threads.
    static inline void initializeDevice(int threads)
    {
        TEUCHOS_ASSERT(threads==1);
    }

    //! \brief Finalize device.
    static inline void finalizeDevice(){}
};

//---------------------------------------------------------------------------//
/*!
 * \class DeviceTraits<Kokkos::OpenMP>
 * \brief Specialization of DeviceTraits for Kokkos::OpenMP.
 */
//---------------------------------------------------------------------------//
template <>
class DeviceTraits<Kokkos::OpenMP>
{
  public:

    //! \brief Initialize device on specified number of threads.
    static inline void initializeDevice(int threads)
    {
        if( Kokkos::OpenMP::is_initialized() )
            Kokkos::OpenMP::finalize();
        Kokkos::OpenMP::initialize( threads );
    }

    //! \brief Finalize device.
    static inline void finalizeDevice()
    {
        Kokkos::OpenMP::finalize();
    }
};

//---------------------------------------------------------------------------//
/*!
 * \class DeviceTraits<Kokkos::Threads>
 * \brief Specialization of DeviceTraits for Kokkos::Threads.
 */
//---------------------------------------------------------------------------//
template <>
class DeviceTraits<Kokkos::Threads>
{
  public:

    //! \brief Initialize device on specified number of threads.
    static inline void initializeDevice(int threads)
    {
        if( Kokkos::Threads::is_initialized() )
            Kokkos::Threads::finalize();
        Kokkos::Threads::initialize( threads );
        Kokkos::Threads::print_configuration(std::cout,true);
    }

    //! \brief Finalize device.
    static inline void finalizeDevice()
    {
        Kokkos::Threads::finalize();
    }
};


} // namespace mcrex

#endif // MCREX_DeviceTraits_hh

