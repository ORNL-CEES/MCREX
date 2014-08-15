//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MonteCarloSolver.cc
 * \author Steven Hamilton
 * \brief  Perform Adjoint Monte Carlo on linear system
 */
//---------------------------------------------------------------------------//

#include <omp.h>

#include <iterator>
#include <string>

#include "MonteCarloSolver.hh"
#include "AdjointMcKernel.hh"
#include "ForwardMcKernel.hh"
#include "PolynomialFactory.hh"

namespace mcrex
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param A  Problem matrix
 * \param pl ParameterList
 *
 * Behavior is controlled by following PL entries on the nested "Monte Carlo"
 * sublist:
 *  - mc_type(string)         : "forward" or ("adjoint")
 *  - estimator(string)       : "collision" or ("expected_value")
 *  - num_histories(int)      : >0 (1000)
 *  - weight_cutoff(SCALAR)   : >0.0 (1.0e-6)
 *  - verbosity(string)       : "none", ("low"), "medium", "high"
 */
//---------------------------------------------------------------------------//
MonteCarloSolver::MonteCarloSolver(Teuchos::RCP<const MATRIX> A,
                                   Teuchos::RCP<Teuchos::ParameterList> pl )
  : MCREX_Solver(A,pl)
{
    // Get Monte Carlo sublist
    Teuchos::RCP<Teuchos::ParameterList> mc_pl =
        Teuchos::sublist(pl,"Monte Carlo");

    // Override verbosity if present on sublist
    MCREX_Solver::setParameters(mc_pl);

    // Determine forward or adjoint
    std::string type = mc_pl->get("mc_type","adjoint");
    if( type == "forward" )
        d_type = FORWARD;
    else
        d_type = ADJOINT;

    // Get parameters off of PL
    std::string estimator = mc_pl->get<std::string>("estimator","expected_value");
    TEUCHOS_TEST_FOR_EXCEPT( estimator != "collision" &&
                             estimator != "expected_value" );
    d_use_expected_value = (estimator == "expected_value");

    if( d_type == FORWARD )
        d_use_expected_value = false;

    // Initialize device node with number of threads
    d_num_threads = mc_pl->get<int>("num_threads",1);
    DeviceTraits<DEVICE>::initializeDevice(d_num_threads);
    d_rng = Teuchos::rcp( new ThreadedRNG(d_num_threads) );

    d_num_histories      = mc_pl->get<int>("num_histories",1000);
    d_weight_cutoff      = mc_pl->get<SCALAR>("weight_cutoff",1.0e-6);
    d_start_wt_factor    = mc_pl->get<SCALAR>("start_weight_factor",1.0);
    d_init_count = 0;
    d_initialized = false;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build data necessary for MC solve
 *
 * This function builds polynomial and MC data based on currently defined
 * matrices.  This is separate from the constructor to allow this object to
 * operate on a different matrix than the one that was used at construction.
 */
//---------------------------------------------------------------------------//
void MonteCarloSolver::initialize()
{
    TEUCHOS_ASSERT( b_A != Teuchos::null );

    // Create Polynomial
    Teuchos::RCP<Polynomial> poly = PolynomialFactory::buildPolynomial(b_A,b_pl);
    TEUCHOS_ASSERT( poly != Teuchos::null );

    // Determine basis
    Teuchos::RCP<Teuchos::ParameterList> poly_pl =
        Teuchos::sublist(b_pl,"Polynomial");
    std::string basis_type = poly_pl->get("polynomial_basis","neumann");
    Teuchos::RCP<PolynomialBasis> basis( new PolynomialBasis(basis_type) );
    if( basis_type == "arbitrary" )
    {
        TEUCHOS_ASSERT( poly_pl->isType<SCALAR>("polynomial_basis_alpha") );
        TEUCHOS_ASSERT( poly_pl->isType<SCALAR>("polynomial_basis_beta") );
        SCALAR alpha = poly_pl->get<SCALAR>("polynomial_basis_alpha");
        SCALAR beta  = poly_pl->get<SCALAR>("polynomial_basis_beta");
        basis->setBasisCoefficients(alpha,beta);
    }

    // Get coefficients of polynomial in desired basis
    d_coeffs = poly->getCoeffs(*basis);
    TEUCHOS_ASSERT( !d_coeffs.is_null() );

    // Get Monte Carlo sublist
    Teuchos::RCP<Teuchos::ParameterList> mc_pl =
        Teuchos::sublist(b_pl,"Monte Carlo");

    // Create Monte Carlo data
    d_mc_data = Teuchos::rcp(
        new MC_Data(b_A,basis,b_pl) );
    convertMatrices(d_mc_data->getIterationMatrix(),
                    d_mc_data->getProbabilityMatrix(),
                    d_mc_data->getWeightMatrix());

    if( d_num_histories%d_num_threads!= 0  && b_verbosity>=LOW )
    {
        std::cout << "WARNING: Requested number of histories ("
            << d_num_histories << ") is not divisible by the number "
            << "of threads (" << d_num_threads << "), ";
        d_num_histories = (d_num_histories/d_num_threads+1)*d_num_threads;
        std::cout << d_num_histories << " histories will be performed."
            << std::endl;
    }

    b_label = "MonteCarloSolver";
    d_initialized = true;
    d_init_count++;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform adjoint Monte Carlo process
 */
//---------------------------------------------------------------------------//
void MonteCarloSolver::applyImpl(const MV &x, MV &y) const
{
    TEUCHOS_ASSERT( d_initialized );

    d_apply_count++;

    // For now we only support operating on a single vector
    TEUCHOS_ASSERT( x.getNumVectors() == 1 );

    LO N = x.getLocalLength();

    // Cache data components of x and y, on-the-fly access is SLOW
    const Teuchos::ArrayRCP<const SCALAR> x_data = x.getData(0);
    const Teuchos::ArrayRCP<SCALAR> y_data = y.getDataNonConst(0);

    Kokkos::ParallelWorkRequest request;
    request.team_size = d_num_threads;
    request.league_size = 1;

    if( d_type == FORWARD )
    {
        int histories_per_state = d_num_histories / N;
        if( d_num_histories % N != 0 )
            histories_per_state++;
        TEUCHOS_ASSERT( histories_per_state > 0 );

        // Construct and execute Kokkos kernel
        ForwardMcKernel kernel(d_P(),d_W(),d_inds(),d_coeffs(),N,x_data(),
                               y_data(), *d_rng,histories_per_state,
                               b_verbosity>=HIGH);

        // Create kernel for performing group of MC histories
        Kokkos::parallel_for( request, kernel );

        SCALAR scale_factor = static_cast<SCALAR>(N) /
                              static_cast<SCALAR>(d_num_histories);
        std::transform(y_data.begin(),y_data.end(),y_data.begin(),
                       [scale_factor](SCALAR x){return x*scale_factor;});

    }
    else if( d_type == ADJOINT )
    {
        // Build initial probability and weight distributions
        Teuchos::ArrayRCP<SCALAR> start_cdf(N);
        Teuchos::ArrayRCP<SCALAR> start_wt(N);
        for( LO i=0; i<N; ++i )
        {
            start_cdf[i] =
                SCALAR_TRAITS::pow(SCALAR_TRAITS::magnitude(x_data[i]),
                                   d_start_wt_factor);
        }
        SCALAR pdf_sum = std::accumulate(start_cdf.begin(),start_cdf.end(),0.0);
        TEUCHOS_ASSERT( pdf_sum > 0.0 );
        std::transform(start_cdf.begin(),start_cdf.end(),start_cdf.begin(),
                       [pdf_sum](SCALAR x){return x/pdf_sum;});
        std::transform(x_data.begin(),x_data.end(),start_cdf.begin(),
                       start_wt.begin(),
                       [](SCALAR x, SCALAR y){return y==0.0 ? 0.0 : x/y;});
        std::partial_sum(start_cdf.begin(),start_cdf.end(),start_cdf.begin());

        int histories_per_thread = d_num_histories / d_num_threads;

        // Create temporary storage on device and mirror it on the host
        Kokkos::View<SCALAR*,DEVICE> y_device("result",N);
        Kokkos::View<SCALAR*,DEVICE>::HostMirror y_mirror =
            Kokkos::create_mirror(y_device);

        // Create kernel for performing group of MC histories
        AdjointMcKernel kernel(d_H(),d_P(),d_W(),d_inds(),d_coeffs(),N,
                               start_cdf(),start_wt(),*d_rng,histories_per_thread,
                               d_use_expected_value,b_verbosity>=HIGH);

        // Execute Kokkos kernel on device
        Kokkos::parallel_reduce( request, kernel, y_mirror.ptr_on_device());

        SCALAR scale_factor = 1.0 / static_cast<SCALAR>(d_num_histories);

        for( LO i=0; i<N; ++i )
        {
            y_data[i] = scale_factor*y_mirror(i);
        }

        // For expected value estimator, need to add C_0*x
        if( d_use_expected_value )
        {
            y.update(d_coeffs[0],x,1.0);
        }
    }

    if( b_verbosity >= LOW )
    {
        std::cout << "Performed " << d_num_histories
            << " histories" << std::endl;
    }

    // There isn't a proper iteration count for MC
    // We use the number of histories as a stand-in
    b_num_iters = d_num_histories;
}

//---------------------------------------------------------------------------//
// PRIVATE MEMBER FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Convert MC matrices from CrsMatrix format to lightweight views.
 *
 * Data access in the Tpetra CrsMatrix implementation is not thread safe so
 * we get threadsafe views here to hand to the MC kernel.
 */
//---------------------------------------------------------------------------//
void MonteCarloSolver::convertMatrices(Teuchos::RCP<const MATRIX> H,
                                       Teuchos::RCP<const MATRIX> P,
                                       Teuchos::RCP<const MATRIX> W)
{
    TEUCHOS_ASSERT( H->isLocallyIndexed() );
    TEUCHOS_ASSERT( P->isLocallyIndexed() );
    TEUCHOS_ASSERT( W->isLocallyIndexed() );
    TEUCHOS_ASSERT( H->supportsRowViews() );
    TEUCHOS_ASSERT( P->supportsRowViews() );
    TEUCHOS_ASSERT( W->supportsRowViews() );

    LO numRows = H->getNodeNumRows();
    d_H.resize(numRows);
    d_P.resize(numRows);
    d_W.resize(numRows);
    d_inds.resize(numRows);

    for( LO irow=0; irow<numRows; ++irow )
    {
        H->getLocalRowView(irow,d_inds[irow],d_H[irow]);
        P->getLocalRowView(irow,d_inds[irow],d_P[irow]);
        W->getLocalRowView(irow,d_inds[irow],d_W[irow]);
    }
}

} // namespace mcrex

