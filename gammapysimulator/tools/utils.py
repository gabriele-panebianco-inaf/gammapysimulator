#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

import os

def get_absolute_path(path):
    """Get absolute path from either $ENV or relative path.
    
    Parameter
    ---------
    path : str
        path to file
    
    Return
    ------
    path : str
        absolute path to file
    """
    if "$" in path:
        path = os.path.expandvars(path)
    elif "./" in path:
        path = os.path.abspath(path)
    return path
