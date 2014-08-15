//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   AdjointMcKernel.cc
 * \author Steven Hamilton
 * \brief  Perform single history of adjoint MC
 */
//---------------------------------------------------------------------------//

#ifndef MCREX_AdjointMcKernel_i_hh
#define MCREX_AdjointMcKernel_i_hh

#include <iterator>
#include <random>

#include "AdjointMcKernel.hh"

namespace mcrex
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param P Views into entries of probability matrix
 * \param W Views into entries of weight matrix
 * \param inds Views into nonzeros indices
 * \param coeffs Polynomial coefficients
 * \param local_length Number of local elements in vector
 * \param start_cdf CDF corresponding to random walk starting locations.
 * \param start_wt  Weights corresponding to random walk starting locations.
 * \param rng ThreadedRNG object for generating local random values
 * \param histories_per_thread Number of histories to be computed
 * \param use_expected_value Should expected value estimator be used?
 * \param print Should debug info be printed?
 */
//---------------------------------------------------------------------------//
AdjointMcKernel::AdjointMcKernel(
        const Teuchos::ArrayView<const Teuchos::ArrayView<const SCALAR> > H,
        const Teuchos::ArrayView<const Teuchos::ArrayView<const SCALAR> > P,
        const Teuchos::ArrayView<const Teuchos::ArrayView<const SCALAR> > W,
        const Teuchos::ArrayView<const Teuchos::ArrayView<const LO> > inds,
        const Teuchos::ArrayView<const SCALAR> coeffs,
        const LO local_length,
        const Teuchos::ArrayView<const SCALAR> start_cdf,
        const Teuchos::ArrayView<const SCALAR> start_wt,
        const ThreadedRNG &rng,
        const int histories_per_thread,
        bool use_expected_value,
        bool print)
  : value_count(local_length)
  , d_H(H)
  , d_P(P)
  , d_W(W)
  , d_inds(inds)
  , d_coeffs(coeffs)
  , d_start_cdf(start_cdf)
  , d_start_wt(start_wt)
  , d_rng(rng)
  , d_histories_per_thread(histories_per_thread)
  , d_use_expected_value(use_expected_value)
  , d_print(print)
  , d_max_history_length(d_coeffs.size()-1)
{
}


//---------------------------------------------------------------------------//
// Kokkos init
//---------------------------------------------------------------------------//
void AdjointMcKernel::init( SCALAR *update ) const
{
    for( LO i=0; i<value_count; ++i )
    {
        update[i] = SCALAR_TRAITS::zero();
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform adjoint Monte Carlo process
 */
//---------------------------------------------------------------------------//
void AdjointMcKernel::operator()(device_type dev, SCALAR *y) const
{
    LO     new_ind;
    GO     state;
    Teuchos::ArrayView<const SCALAR> row_h, row_cdf, row_wts;
    Teuchos::ArrayView<const GO>     row_inds;

    for( int ihist=0; ihist<d_histories_per_thread; ++ihist )
    {
        // Get starting position and weight
        state = getNewState(d_start_cdf,dev.team_rank());
        if( state == GO_TRAITS::invalid() )
            continue;

        if( SCALAR_TRAITS::magnitude(d_start_wt[state]) == SCALAR_TRAITS::zero())
            continue;

        SCALAR weight = d_start_wt[state];
        SCALAR initial_weight = weight;

        if( d_print )
        {
            std::cout << "Starting history in state " << state
                << " with initial weight " << initial_weight << std::endl;
        }

        // Collision estimator starts tallying on zeroth order term
        // Expected value estimator gets this added explicitly at the end
        int stage = 0;
        if( d_use_expected_value )
            stage++;

        // Transport particle until done
        while(true)
        {
            // Get data and add to tally
            getNewRow(state,row_h,row_cdf,row_wts,row_inds);
            tallyContribution(state,d_coeffs[stage]*weight,
                              row_h,row_inds,y);

            if( stage >= d_max_history_length )
                break;

            // Get new state index
            new_ind = getNewState(row_cdf,dev.team_rank());
            if( new_ind == GO_TRAITS::invalid() )
                break;

            // Modify weight and update state
            weight *=  row_wts[new_ind];
            state   = row_inds[new_ind];
            stage++;

            if( d_print )
            {
                std::cout << "Transitioning to state " << state
                    << " with new weight " << weight << std::endl;
            }

        } // while
    } // for ihist
}

//---------------------------------------------------------------------------//
// Kokkos join
//---------------------------------------------------------------------------//
void AdjointMcKernel::join(      volatile SCALAR *update,
                           const volatile SCALAR *input) const
{
    for( LO i=0; i<value_count; ++i )
    {
        update[i] = update[i] + input[i];
    }
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Tally contribution into vector
 */
//---------------------------------------------------------------------------//
void AdjointMcKernel::getNewRow( const GO state,
        Teuchos::ArrayView<const SCALAR> &h_vals,
        Teuchos::ArrayView<const SCALAR> &p_vals,
        Teuchos::ArrayView<const SCALAR> &w_vals,
        Teuchos::ArrayView<const LO>     &inds) const
{
    h_vals = d_H[state];
    p_vals = d_P[state];
    w_vals = d_W[state];
    inds   = d_inds[state];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tally contribution into vector
 */
//---------------------------------------------------------------------------//
void AdjointMcKernel::tallyContribution(
        const GO state, const SCALAR wt,
        const Teuchos::ArrayView<const SCALAR> h_vals,
        const Teuchos::ArrayView<const LO>     inds,
        SCALAR *y) const
{
    if( d_use_expected_value )
    {
        LO num_entries = inds.size();

        if( num_entries > 0 )
        {
            y[inds[0]] += wt*h_vals[0];
            for( LO i=1; i<num_entries; ++i )
            {
                // P is cdf, we want value of pdf
                y[inds[i]] += wt*h_vals[i];
            }
        }
    }
    else
    {
        y[state] += wt;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get new state by sampling from cdf
 */
//---------------------------------------------------------------------------//
GO AdjointMcKernel::getNewState(const Teuchos::ArrayView<const SCALAR> cdf,
                                const int thread) const
{
    // Generate random number
    double rand = d_rng.getRandom(thread);

    // Sample cdf to get new state
    Teuchos::ArrayView<const SCALAR>::iterator elem =
        std::lower_bound(cdf.begin(),cdf.end(),rand);

    if( elem == cdf.end() )
        return GO_TRAITS::invalid();

    return elem - cdf.begin();
}

} // namespace mcrex

#endif // MCREX_AdjointMcKernel_i_hh
