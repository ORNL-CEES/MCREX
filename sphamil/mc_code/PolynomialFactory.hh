
#ifndef MCREX_PolynomialFactory_hh
#define MCREX_PolynomialFactory_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "AnasaziTypes.hpp"

#include "MCREX_Typedefs.hh"

#include "Polynomial.hh"

namespace mcrex
{

//---------------------------------------------------------------------------//
/*!
 * \class PolynomialFactory
 * \brief Create Polynomial object
 */
//---------------------------------------------------------------------------//
class PolynomialFactory
{

  private:

    // No construction
    PolynomialFactory(){};

  public:

    static Teuchos::RCP<Polynomial>
    buildPolynomial(Teuchos::RCP<const MATRIX> A,
                    Teuchos::RCP<Teuchos::ParameterList> pl );
};

}

#endif // MCREX_PolynomialFactory_hh

