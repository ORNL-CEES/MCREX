//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   testMonteCarloSolver.cc
 * \author Steven Hamilton
 * \brief  Test of MonteCarloSolver class.
 */
//---------------------------------------------------------------------------//

#include "gtest/gtest.h"

#include <time.h>

#include "LinearSystem.hh"
#include "LinearSystemFactory.hh"
#include "MonteCarloSolver.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "MCREX_Typedefs.hh"

using namespace mcrex;

TEST(MonteCarloSolver, Basic)
{
    Teuchos::TestForException_setEnableStacktrace(true);

    // Create ParameterList
    Teuchos::RCP<Teuchos::ParameterList> pl( new Teuchos::ParameterList() );
    Teuchos::RCP<Teuchos::ParameterList> mat_pl =
        Teuchos::sublist(pl,"Problem");
    Teuchos::RCP<Teuchos::ParameterList> mc_pl =
        Teuchos::sublist(pl,"Monte Carlo");
    Teuchos::RCP<Teuchos::ParameterList> poly_pl =
        Teuchos::sublist(pl,"Polynomial");

    mat_pl->set("matrix_type","laplacian");
    mat_pl->set("matrix_size",10);

    mc_pl->set("estimator","expected_value");
    mc_pl->set("num_histories",20000);
    mc_pl->set("weight_cutoff",0.01);
    mc_pl->set("verbosity","medium");

    poly_pl->set("polynomial_order",100);

    Teuchos::RCP<mcrex::LinearSystem> system =
        mcrex::LinearSystemFactory::buildLinearSystem(pl);
    Teuchos::RCP<const MATRIX> A = system->getMatrix();
    Teuchos::RCP<const MV> b = system->getRhs();

    // Test collision estimator
    mc_pl->set("estimator","collision");
    Teuchos::RCP<mcrex::MonteCarloSolver> solver(
        new mcrex::MonteCarloSolver(A,pl) );
    solver->compute();

    Teuchos::RCP<MV> x( new MV(A->getDomainMap(),1) );

    solver->apply(*b,*x);

    // Compute final residual
    Teuchos::RCP<MV> r( new MV(A->getDomainMap(),1) );
    A->apply(*x,*r);
    r->update(1.0,*b,-1.0);
    Teuchos::ArrayRCP<SCALAR> res_norm(1), b_norm(1);
    r->norm2(res_norm());
    b->norm2(b_norm());
    std::cout << "Final relative residual norm: "
              << res_norm[0]/b_norm[0] << std::endl;

    // This should *almost* always pass
    EXPECT_TRUE( res_norm[0]/b_norm[0] < 0.12 );

    // Explicitly destroy solver before building new one
    // Need old destructor to be called before new constructor
    solver = Teuchos::null;

    // Test expected value estimator
    mc_pl->set("estimator","expected_value");
    solver = Teuchos::RCP<mcrex::MonteCarloSolver>(
        new mcrex::MonteCarloSolver(A,pl) );
    solver->compute();

    solver->apply(*b,*x);

    // Compute final residual
    A->apply(*x,*r);
    r->update(1.0,*b,-1.0);
    r->norm2(res_norm());
    b->norm2(b_norm());
    std::cout << "Final relative residual norm: "
              << res_norm[0]/b_norm[0] << std::endl;

    // This should *almost* always pass
    EXPECT_TRUE( res_norm[0]/b_norm[0] < 0.06 );
}

