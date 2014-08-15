//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MCREX_Solver.hh
 * \author Steven Hamilton
 * \brief  Base class for MCREX solvers
 */
//---------------------------------------------------------------------------//

#ifndef MCREX_MCREX_Solver_hh
#define MCREX_MCREX_Solver_hh

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "MCREX_Typedefs.hh"

namespace mcrex
{

//---------------------------------------------------------------------------//
/*!
 * \class MCREX_Solver
 * \brief Base class for all linear solvers.
 *
 * This class defines the minimal interface for all MCREX solvers.
 * This class in turn derives from Tpetra::Operator, which allows for very
 * simple syntax and additionally allows any MCREX_Solver to be used directly
 * as a preconditioner with any Belos solver.
 */
//---------------------------------------------------------------------------//
class MCREX_Solver : virtual public OP
{
  public:

    //! Enumeration describing level of output produced by solver
    enum Verbosity_Level
    {
        //! Produce no output.
        NONE,
        //! Minimal output, usually only a final summary.
        LOW,
        //! More output, commonly updates on the iteration process.
        MEDIUM,
        /*! Significant output, possibly extra information about each
            iteration. */
        HIGH,
        //! Highest level of output, may only be meaningful to developers.
        DEBUG
    };

    // Constructor
    MCREX_Solver(Teuchos::RCP<const MATRIX> A,
                 Teuchos::RCP<Teuchos::ParameterList> pl);

    virtual ~MCREX_Solver(){}

    //! Return number of iterations for last call to solve.
    LO getNumIters(){ return b_num_iters; }

    //! Return domain map of solver (this is range map of operator.)
    virtual Teuchos::RCP<const MAP> getDomainMap() const
    {
        return b_A->getRangeMap();
    }
    //! Return range map of solver (this is domain map of operator.)
    virtual Teuchos::RCP<const MAP> getRangeMap() const
    {
        return b_A->getDomainMap();
    }

    //! Is transpose apply available.
    virtual bool hasTransposeApply() const { return false; }

    // Apply
    virtual void apply(const MV &x, MV &y,
                       Teuchos::ETransp mode=Teuchos::NO_TRANS,
                       SCALAR alpha=SCALAR_TRAITS::one(),
                       SCALAR beta=SCALAR_TRAITS::zero()) const;

  protected:

    // Implementation of apply provided by derived class
    virtual void applyImpl(const MV &x, MV &y) const = 0;

    // Set verbosity level from PL
    void setParameters( Teuchos::RCP<Teuchos::ParameterList> pl );

    // Problem matrix
    Teuchos::RCP<const MATRIX> b_A;
    Teuchos::RCP<Teuchos::ParameterList> b_pl;

    std::string b_label;
    mutable LO b_num_iters; // Allow setting this in const Apply call

    LO b_max_iterations;
    MAGNITUDE b_tolerance;
    Verbosity_Level b_verbosity;
};

}

#endif // MCREX_MCREX_Solver_hh

