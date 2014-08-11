//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PolynomialUtils.hh
 * \author Steven Hamilton
 * \brief  Pure static class for various polynomial-related utilities.
 */
//---------------------------------------------------------------------------//

#ifndef MCREX_PolynomialUtils_hh
#define MCREX_PolynomialUtils_hh

#include <cmath>

#include "Teuchos_Assert.hpp"
#include "Teuchos_ArrayRCP.hpp"

#include "MCREX_Typedefs.hh"

namespace mcrex
{
//---------------------------------------------------------------------------//
/*!
 * \class PolynomialUtils
 * \brief Helper functions for polynomial construction.
 */
//---------------------------------------------------------------------------//
class PolynomialUtils
{
  public:

    //! Binomial coefficient (choose function)
    static unsigned long choose(unsigned long n, unsigned long k);

    //! Compute coefficients for specified order Chebyshev polynomial
    static Teuchos::ArrayRCP<const SCALAR> getChebyshevCoefficients(LO n);
};

}

#endif // MCREX_PolynomialUtils_hh

