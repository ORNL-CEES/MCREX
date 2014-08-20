#!/bin/bash

qsub run_profugus_32.pbs
qsub run_profugus_64.pbs
qsub run_profugus_128.pbs
qsub run_profugus_256.pbs
qsub run_profugus_496.pbs
qsub run_profugus_976.pbs