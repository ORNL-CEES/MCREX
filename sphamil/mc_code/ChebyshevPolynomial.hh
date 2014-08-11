//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ChebyshevPolynomial.hh
 * \author Steven Hamilton
 * \brief  ChebyshevPolynomial class declaration.
 */
//---------------------------------------------------------------------------//

#ifndef MCREX_ChebyshevPolynomial_hh
#define MCREX_ChebyshevPolynomial_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "MCREX_Typedefs.hh"

#include "PolynomialBasis.hh"
#include "Polynomial.hh"

namespace mcrex
{

//---------------------------------------------------------------------------//
/*!
 * \class ChebyshevPolynomial
 * \brief Chebyshev series polynomial.
 */
//---------------------------------------------------------------------------//
class ChebyshevPolynomial : public Polynomial
{

  public:

    // Constructor
    ChebyshevPolynomial(Teuchos::RCP<const MATRIX> A,
               Teuchos::RCP<Teuchos::ParameterList> pl);

    //! \brief Access minimum eigenvalue.
    SCALAR getLambdaMin() const
    {
        TEUCHOS_ASSERT( !SCALAR_TRAITS::isnaninf(d_lambda_min) );
        return d_lambda_min;
    }

    //! \brief Access maximum eigenvalue.
    SCALAR getLambdaMax() const
    {
        TEUCHOS_ASSERT( !SCALAR_TRAITS::isnaninf(d_lambda_max) );
        return d_lambda_max;
    }

  private:

    void computeEigenvalues();

    // Eigenvalue data
    SCALAR d_lambda_min;
    SCALAR d_lambda_max;

};

}

#endif // MCREX_ChebyshevPolynomial_hh

