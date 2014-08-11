//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LinearSolverFactory.hh
 * \author Steven Hamilton
 * \brief  Construct Tpetra_CrsMatrix from ParameterList.
 */
//---------------------------------------------------------------------------//

#ifndef MCREX_LinearSolverFactory_hh
#define MCREX_LinearSolverFactory_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "MCREX_Typedefs.hh"
#include "MCREX_Solver.hh"

namespace mcrex
{

//---------------------------------------------------------------------------//
/*!
 * \class LinearSolverFactory
 * \brief Construct linear solver from MATRIX and ParameterList.
 */
//---------------------------------------------------------------------------//
class LinearSolverFactory
{
  private:

    // Pure static, disallow construction
    LinearSolverFactory(){};

  public:

    static Teuchos::RCP<MCREX_Solver> buildSolver(
        std::string solver_type,
        Teuchos::RCP<const MATRIX> A,
        Teuchos::RCP<Teuchos::ParameterList> pl );

};

}

#endif // MCREX_LinearSolverFactory_hh

