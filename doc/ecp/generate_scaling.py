###############################################################################
# File  : ecp/generate_scaling.py
# Author: Thomas M. Evans
# Date  : Mon Jan 25 22:55:47 2016
###############################################################################
from __future__ import (division, absolute_import, print_function, )
#-----------------------------------------------------------------------------#
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('/codes/Exnihilo/Exnihilo/python')

import exnihilotools.matplotlib as ext
ext.screen_style()
###############################################################################

def main():

    str_cores = [256, 1024, 4096, 7744, 10816]
    str_time  = [515.52, 122.76, 27.96, 17.72, 13.72]
    str_eff   = [1.00, 1.05, 1.15, 0.96, 0.89]

    wk_cores = [64, 256, 1024, 4096, 7744, 10816]
    wk_time  = [12.65, 12.86, 14.59, 14.25, 14.75, 13.72]
    wk_eff   = [1.00, 0.98, 0.87, 0.89, 0.86, 0.92]

    # Strong Scaling - time
    plt.loglog(str_cores, str_time, 'o-')
    plt.grid()
    plt.xlabel("cores")
    plt.ylabel("$t$ (sec)")
    plt.savefig('strong_dd_time.pdf', bbox_inches='tight')
    plt.close()

    # Strong Scaling - efficiency
    plt.plot(str_cores, str_eff, 'o-')
    plt.grid()
    plt.xlabel("cores")
    plt.ylabel("$\epsilon$")
    plt.ylim(0.5, 1.5)
    plt.savefig('strong_dd_eff.pdf', bbox_inches='tight')
    plt.close()

    # Weak Scaling - time
    plt.plot(wk_cores, wk_time, 'o-')
    plt.grid()
    plt.xlabel("cores")
    plt.ylabel("$t$ (sec)")
    plt.ylim(6,18)
    plt.savefig('weak_dd_time.pdf', bbox_inches='tight')
    plt.close()

    # Weak Scaling - efficiency
    plt.plot(wk_cores, wk_eff, 'o-')
    plt.grid()
    plt.xlabel("cores")
    plt.ylabel("$\epsilon$")
    plt.ylim(0.5, 1.5)
    plt.savefig('weak_dd_eff.pdf', bbox_inches='tight')
    plt.close()

    # Strong scaling - replication
    str_cores = [1, 2, 3, 4]
    str_time  = [12.65, 6.72, 4.66, 3.71]
    str_eff   = [1.00, 0.92, 0.89, 0.83]

    # Strong Scaling - time
    plt.semilogy(str_cores, str_time, 'o-')
    plt.grid()
    plt.xlabel("Subsets")
    plt.ylabel("$t$ (sec)")
    plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    plt.savefig('strong_rep_time.pdf', bbox_inches='tight')
    plt.close()

    # Strong Scaling - efficiency
    plt.plot(str_cores, str_eff, 'o-')
    plt.grid()
    plt.xlabel("Subsets")
    plt.ylabel("$\epsilon$")
    plt.ylim(0.5, 1.5)
    plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    plt.savefig('strong_rep_eff.pdf', bbox_inches='tight')
    plt.close()

#-----------------------------------------------------------------------------#
if __name__ == '__main__':
    main()

###############################################################################
# end of ecp/generate_scaling.py
###############################################################################
