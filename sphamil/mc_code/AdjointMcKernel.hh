//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   AdjointMcKernel.hh
 * \author Steven Hamilton
 * \brief  Perform adjoint MC histories to solve linear system.
 */
//---------------------------------------------------------------------------//

#ifndef MCREX_AdjointMcKernel_hh
#define MCREX_AdjointMcKernel_hh

#include "MC_Data.hh"
#include "MCREX_Solver.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Kokkos_Parallel.hpp"

#include "ThreadedRNG.hh"
#include "MCREX_Typedefs.hh"

namespace mcrex
{

//---------------------------------------------------------------------------//
/*!
 * \class AdjointMcKernel
 * \brief Perform Monte Carlo random walks on linear system.
 *
 * This class performs random walks using the adjoint Monte Carlo algorithm.
 * The interface of this function conforms to the Kokkos "parallel_reduce"
 * functor API to enable automated shared memory parallelism over MC histories.
 */
//---------------------------------------------------------------------------//

class AdjointMcKernel
{
  public:

    // Required typedefs for Kokkos functor API

    //! Type of device where kernel will be executed
    typedef DEVICE  device_type;
    //! Type of data kernel will operator on
    typedef SCALAR  value_type[];

    // Required public member variable for Kokkos array functor API
    //! Number of entries in a value_type array
    const LO value_count;

    AdjointMcKernel(const Teuchos::ArrayView<const Teuchos::ArrayView<const SCALAR> > H,
                    const Teuchos::ArrayView<const Teuchos::ArrayView<const SCALAR> > P,
                    const Teuchos::ArrayView<const Teuchos::ArrayView<const SCALAR> > W,
                    const Teuchos::ArrayView<const Teuchos::ArrayView<const LO> > inds,
                    const Teuchos::ArrayView<const SCALAR>  coeffs,
                    const LO local_length,
                    const Teuchos::ArrayView<const SCALAR>  start_cdf,
                    const Teuchos::ArrayView<const SCALAR>  start_wt,
                    const ThreadedRNG &rng,
                    const int histories_per_thread,
                    bool use_expected_value,
                    bool print);

    // Kokkos "parallel_reduce" API functions

    //! \brief Initialize a thread
    KOKKOS_INLINE_FUNCTION
    void init( SCALAR *update ) const;

    //! \brief Compute kernel
    KOKKOS_INLINE_FUNCTION
    void operator()(device_type dev, SCALAR *y) const;

    //! \brief Join threads together via a reduce
    KOKKOS_INLINE_FUNCTION
    void join(      volatile SCALAR *update,
              const volatile SCALAR *input) const;

  private:

    inline void getNewRow(const GO state,
                          Teuchos::ArrayView<const SCALAR> &h_vals,
                          Teuchos::ArrayView<const SCALAR> &p_vals,
                          Teuchos::ArrayView<const SCALAR> &w_vals,
                          Teuchos::ArrayView<const LO>     &inds) const;
    inline void tallyContribution(const GO state, const SCALAR wt,
                                  const Teuchos::ArrayView<const SCALAR> h_vals,
                                  const Teuchos::ArrayView<const LO>     inds,
                                  SCALAR *y) const;
    inline GO getNewState(const Teuchos::ArrayView<const SCALAR> cdf,
                          const int thread ) const;


    // Warning, these are non reference-counted views
    // The underlying reference-counted objects from which they were
    // created must stay alive for the entire life of this object
    const Teuchos::ArrayView<const Teuchos::ArrayView<const SCALAR> > d_H;
    const Teuchos::ArrayView<const Teuchos::ArrayView<const SCALAR> > d_P;
    const Teuchos::ArrayView<const Teuchos::ArrayView<const SCALAR> > d_W;
    const Teuchos::ArrayView<const Teuchos::ArrayView<const LO> > d_inds;
    const Teuchos::ArrayView<const SCALAR> d_coeffs;
    const Teuchos::ArrayView<const SCALAR> d_start_cdf;
    const Teuchos::ArrayView<const SCALAR> d_start_wt;

    const ThreadedRNG d_rng;

    const int d_histories_per_thread;
    const bool d_use_expected_value;
    const bool d_print;
    const int d_max_history_length;

};

} // namespace mcrex

#include "AdjointMcKernel.i.hh"

#endif // MCREX_MonteCarloSolver_hh

