//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   AdditiveSchwarzWrapper.hh
 * \author Steven Hamilton
 * \brief  Wrap a Belos solver into an operator.
 */
//---------------------------------------------------------------------------//

#ifndef MCREX_AdditiveSchwarzWrapper_hh
#define MCREX_AdditiveSchwarzWrapper_hh

#include "MCREX_Solver.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "MCREX_Typedefs.hh"

#include "Ifpack2_AdditiveSchwarz.hpp"

namespace mcrex
{

//---------------------------------------------------------------------------//
/*!
 * \class AdditiveSchwarzWrapper
 * \brief Wrap an Ifpack2 AdditiveSchwarz object into an MCREX_Solver
 */
//---------------------------------------------------------------------------//
class AdditiveSchwarzWrapper : public MCREX_Solver
{
  public:

    AdditiveSchwarzWrapper(Teuchos::RCP<const MATRIX> A,
        Teuchos::RCP<Ifpack2::Preconditioner<SCALAR,LO,GO,NODE> > prec,
        Teuchos::RCP<Teuchos::ParameterList> pl);

  private:

    // Implementation of apply
    void applyImpl(const MV &x, MV &y) const;

    // Belos solver
    Teuchos::RCP<Ifpack2::AdditiveSchwarz<CRS_MATRIX> > d_schwarz;
};

}

#endif // MCREX_AdditiveSchwarzWrapper_hh

