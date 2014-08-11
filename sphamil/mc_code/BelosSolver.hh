//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   BelosSolver.hh
 * \author Steven Hamilton
 * \brief  Wrap a Belos solver into an operator.
 */
//---------------------------------------------------------------------------//

#ifndef MCREX_BelosSolver_hh
#define MCREX_BelosSolver_hh

#include "MCREX_Solver.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "MCREX_Typedefs.hh"

#include "BelosTpetraAdapter.hpp"
#include "BelosSolverManager.hpp"

namespace mcrex
{

//---------------------------------------------------------------------------//
/*!
 * \class BelosSolver
 * \brief Wrap a Belos solver into an operator.
 *
 * This class provides a thin wrapper around a Belos linear solver
 * (i.e. a Belos::SolverManager).  The "Belos" sublist of the provided
 * ParameterList is provided directly to the Belos::SolverFactory, so any
 * Belos solvers and corresponding options are automatically available.
 * Consult Belos documentation for full list of solver options.
 */
//---------------------------------------------------------------------------//
class BelosSolver : public MCREX_Solver
{
  public:

    BelosSolver(Teuchos::RCP<const MATRIX> A,
                Teuchos::RCP<Teuchos::ParameterList> pl);

  private:

    // Implementation of apply
    void applyImpl(const MV &x, MV &y) const;

    // Preconditioner
    Teuchos::RCP<OP> d_prec;

    // Belos solver
    Teuchos::RCP<Belos::SolverManager<SCALAR,MV,OP> > d_belos_solver;
};

}

#endif // MCREX_BelosSolver_hh

