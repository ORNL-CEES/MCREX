MCREX Status
============

The MCREX project has 3 objectives for the development of Monte Carlo solver
methods for linear systems:

 * Analyze and improve the iterative performance of Monte Carlo algorithms for
   linear systems
 * Analyze and characterize the resilient aspects of Monte Carlo algorithms
   for linear systems
 * Develop and model highly scalable parallel algorithms for the
   implementation of Monte Carlo algorithms for linear systems

The end goal is to develop a set of resilient, scalable algorithms for the
solution of linear systems.  We have made significant progress in each of
these areas as summarized below.

Algorithmic Performance
-----------------------

We have analyzed the iterative behavior of the basic Monte Carlo Neumann-Ulam
method, and we have developed an accelerated iterative sequence for this
iteration termed Monte Carlo Synthetic Acceleration (MCSA). We have analyzed
and developed a set of polynomial preconditioners for further improving the
efficiency of the MCSA method. Also, we have tested performance with respect
to different termination strategies for the random walk part of the
algorithm.  These improvements have dramatically improved the iterative
performance of MCSA.  To further improve efficiency, our project collaborators
at Emory have been investigating sparse approximate inverse preconditioners
that have the benefit of decreasing the spectral radius of the linear system
without increasing the sparsity, a necessary condition for the application of
Monte Carlo methods to linear systems.  Finally, we have developed metrics for
automatically determining the history and sample length of a given iterative
sequence.

Resiliency
----------

We have developed an open-source version of the MCSA solver (MCLS) that is
available on GitHub (<github.com/ORNL-CEES/Profugus).  We have begun to
implement this solver package in the xSIM MPI fault simulator to test
algorithmic behavior to hard faults. The development of xSIM for this work
featured in Christian Engelmann's DOE Early Career Science award for FY16.

Parallel Performance
--------------------

The MCLS package is designed using a Multiset-Overlapping Domain (MSOD)
parallel strategy that is ammenable to large-scale systems.  We have also
developed a GPU accelerated kernel for on-node performance on heterogeneous
computing systems.  The GPU kernel has resulted in a factor of 80 speedup over
CPU-only use.  We have also developed a performance model for MSOD that will
be subject of a paper that is currently in preparation.

Presentations
-------------

  * "Parallel Algorithms for the Monte Carlo Synthetic Acceleration Linear
    Solver Method" (S. Slattery, T. Evans, S. Hamilton), SIAM CSE, Salt Lake
    City, 2015.
  * "Iterative Performance of Monte Carlo Linear Solver Methods"
    (M.L. Pasini), SIAM CSE, Salt Lake City, 2015.
  * "A multilevel Monte Carlo method for linear systems"
    (S. Slattery, T. Evans, S. Hamilton), Copper Mountain Conference on
    Iterative Solvers, Copper MT, CO, 2014.
  * "Monte Carlo Synthetic Acceleration as approximate polynomial
    preconditioning" (S. Hamilton, T. Evans, S. Slattery), Copper Mountain
    Conference on Iterative Solvers, Copper MT, CO, 2014.
  * "Monte Carlo Linear Solvers" (S. Slattery, T. Evans, S. Hamilton), NC
    State Seminar, Nov. 4, 2014.
  * "Monte Carlo Linear Solvers" (S. Hamilton, T. Evans, S. Slattery), Emory
    Mathematics Dept. Seminar, Sept. 19, 2014.
  * "Monte Carlo Methods for Linear Systems" (T. Evans), MIT Nuclear
    Engineering Dept Seminar Series, Sept. 30, 2014.

Publications
------------

  * S.R. Slattery, T.M. Evans, P.P.H. Wilson. A Spectral Analysis of the
    Domain Decomposed Monte Carlo Method for Linear Systems. *Nuclear
    Engineering Design*,  DOI:10.1016/j.nucengdes.2015.07.054, 2015.
  * T.M. Evans, S.W. Mosher, S.R. Slattery, and S.P. Hamilton. A Monte Carlo
    Synthetic Acceleration Method for Solving the Thermal Radiation Diffusion
    Equation. *J. Comp. Phys.*, **258**, 338-358, 2014.
