#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

#!/bin/bash

set -e 

SCRIPT_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export GAMMAPYSIMULATOR=$SCRIPT_DIRECTORY

runsimulator.py -conf $1
