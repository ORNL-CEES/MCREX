//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PolynomialPreconditioner.hh
 * \author Steven Hamilton
 * \brief  Perform Chebyshev iteration.
 */
//---------------------------------------------------------------------------//

#ifndef MCREX_PolynomialPreconditioner_hh
#define MCREX_PolynomialPreconditioner_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "MCREX_Typedefs.hh"

#include "MCREX_Solver.hh"

namespace mcrex
{

//---------------------------------------------------------------------------//
/*!
 * \class PolynomialPreconditioner
 * \brief Apply matrix polynomial.
 *
 * This class generates a matrix polynomial (using a Polynomial
 * constructed by the PolynomialFactory class)
 * and applies it to a vector using Horner's method.  This class is intended
 * to be used as a preconditioner to a linear solver.
 */
//---------------------------------------------------------------------------//
class PolynomialPreconditioner : public MCREX_Solver
{
  public:

    PolynomialPreconditioner(Teuchos::RCP<const MATRIX> A,
                             Teuchos::RCP<Teuchos::ParameterList> pl);

  private:

    // Implementation of apply
    void applyImpl(const MV &x, MV &y) const;

    Teuchos::ArrayRCP<const SCALAR> d_coeffs;
};

}

#endif // MCREX_PolynomialPreconditioner_hh

