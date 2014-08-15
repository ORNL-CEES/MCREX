//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ChebyshevIteration.hh
 * \author Steven Hamilton
 * \brief  Perform Chebyshev iteration.
 */
//---------------------------------------------------------------------------//

#ifndef MCREX_ChebyshevIteration_hh
#define MCREX_ChebyshevIteration_hh

#include "MCREX_Solver.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "MCREX_Typedefs.hh"

namespace mcrex
{

//---------------------------------------------------------------------------//
/*!
 * \class ChebyshevIteration
 * \brief Perform Chebyshev iteration to solve a linear system
 *
 * The iteration in this class is the implementation given by Algorithm 4
 * of Gutknecht & Roellin, "The Chebyshev iteration revisited."
 * Specifically, updates to the solution vector are made using a "delta"
 * formulation and residuals are computed explicitly.
 */
//---------------------------------------------------------------------------//

class ChebyshevIteration : public MCREX_Solver
{
  public:

    ChebyshevIteration(Teuchos::RCP<const MATRIX> A,
                        Teuchos::RCP<Teuchos::ParameterList> pl);

  private:

    // Implementation of apply
    void applyImpl(const MV &x, MV &y) const;

    SCALAR d_c;
    SCALAR d_d;
};

}

#endif // MCREX_ChebyshevIteration_hh

