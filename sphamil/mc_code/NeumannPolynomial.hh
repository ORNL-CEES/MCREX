//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   NeumannPolynomial.hh
 * \author Steven Hamilton
 * \brief  NeumannPolynomial class declaration.
 */
//---------------------------------------------------------------------------//

#ifndef MCREX_NeumannPolynomial_hh
#define MCREX_NeumannPolynomial_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "MCREX_Typedefs.hh"

#include "PolynomialBasis.hh"
#include "Polynomial.hh"

namespace mcrex
{

//---------------------------------------------------------------------------//
/*!
 * \class NeumannPolynomial
 * \brief Neumann series polynomial.
 */
//---------------------------------------------------------------------------//
class NeumannPolynomial : public Polynomial
{

  public:

    // Constructor
    NeumannPolynomial(Teuchos::RCP<const MATRIX> A,
               Teuchos::RCP<Teuchos::ParameterList> pl);
};

}

#endif // MCREX_NeumannPolynomial_hh

