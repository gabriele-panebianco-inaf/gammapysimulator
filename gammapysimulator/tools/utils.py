#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

import matplotlib.pyplot as plt
import numpy as np
import os
import shutil

from gammapy.utils.table import table_from_row_data
from scipy import interpolate

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

def clean_directory(directory, logger):
    """
    Make an empty directory. Delete it and re-make it if it already exists.
    """
    if os.path.isdir(directory):
        logger.warning(f"Remove existing directory {directory}")
        shutil.rmtree(directory)
    logger.info(f"Create directory {directory}")
    os.makedirs(directory)
    return None


def info_table(datasets, cumulative=False):
        """Get info table for datasets.

        Parameters
        ----------
        datasets : gammapy.datasets.Datasets
            Datasets.
        cumulative : bool
            Cumulate info across all observations

        Returns
        -------
        info_table : `~astropy.table.Table`
            Info table.
        """
        if not datasets.is_all_same_type:
            raise ValueError("Info table not supported for mixed dataset type.")

        name = "stacked" if cumulative else datasets[0].name
        stacked = datasets[0].to_masked(name=name)
        stacked.models = datasets[0].models

        rows = [stacked.info_dict()]

        for dataset in datasets[1:]:
            if cumulative:
                stacked.stack(dataset)
                row = stacked.info_dict()
            else:
                row = dataset.info_dict()

            rows.append(row)

        return table_from_row_data(rows=rows)
    

def InterpolateFunction(OldValues, OldAxis, NewAxis, scale='lin', kind='linear'):
    """
    Interpolate differential quantities to find their values in a new axis.
    
    Parameters
    ----------
    OldValues : astropy.Quantity
        Quantity evaluated at old axis.
    OldAxis, NewAxis : astropy.Quantity
        Old and New axis where to evaluate the quantity values.
    scale : str
        Either 'lin' or 'log' to evaluate quantities at lin or log axes bins.
    kind : str
        Argument passed to scipy.interpolate.interp1d
        
    Returns
    -------
    Newvalues : astropy.Quantity
        Quantity evaluated at new axis.
    """
    
    # Check New Bounds are within Old Bounds
    if NewAxis[0 ] < OldAxis[0 ]:
        raise ValueError(f"Requested function evaluation at {NewAxis[0 ]}. The minimum is {OldAxis[0 ]}.")
    if NewAxis[-1] > OldAxis[-1]:
        raise ValueError(f"Requested function evaluation at {NewAxis[-1]}. The maximum is {OldAxis[-1]}.")
    
    # Function evaluation at linear or logarithmic bins
    if scale=="lin":
        xold = OldAxis.value
        xnew = NewAxis.value
    elif scale=="log":
        xold = np.log10(OldAxis.value)
        xnew = np.log10(NewAxis.value)
    else:
        raise ValueError("Scale must be either lin or log.")
    yold = OldValues.value
    
    # Perform interpolation
    f = interpolate.interp1d(xold, yold, kind=kind)
    
    # Define new values
    ynew = f(xnew)
    
    # Add unit
    NewValues = ynew * OldValues.unit
    
    return NewValues



def InterpolateMap(OldValues, OldAxis1, OldAxis2, NewAxis1, NewAxis2, scale1='lin', scale2='lin', method='linear'):
    """
    Interpolate 2D differential quantities to find their values in a new axis.
    
    Parameters
    ----------
    OldValues : astropy.Quantity
        Quantity evaluated at old axis.
    OldAxis1, NewAxis1, OldAxis2, NewAxis2 : astropy.Quantity
        Old and New axis where to evaluate the quantity values.
    scale1, scale2 : str
        Either 'lin' or 'log' to evaluate quantities at lin or log axes bins.
    method : str
        Argument passed to scipy.interpolate.RegularGridInterpolator
        
    Returns
    -------
    Newvalues : astropy.Quantity
        Quantity evaluated at new axis.
    """
    
    # Check New Bounds are within Old Bounds
    if NewAxis1[0 ] < OldAxis1[0 ]:
        raise ValueError(f"Requested function evaluation at {NewAxis1[0 ]}. The minimum is {OldAxis1[0 ]}.")
    if NewAxis1[-1] > OldAxis1[-1]:
        raise ValueError(f"Requested function evaluation at {NewAxis1[-1]}. The maximum is {OldAxis1[-1]}.")
    if NewAxis2[0 ] < OldAxis2[0 ]:
        raise ValueError(f"Requested function evaluation at {NewAxis2[0 ]}. The minimum is {OldAxis2[0 ]}.")
    if NewAxis2[-1] > OldAxis2[-1]:
        raise ValueError(f"Requested function evaluation at {NewAxis2[-1]}. The maximum is {OldAxis2[-1]}.")
    
    # Function evaluation at linear or logarithmic bins
    if scale1=="lin":
        xold1 = OldAxis1.value
        xnew1 = NewAxis1.value
    elif scale1=="log":
        xold1 = np.log10(OldAxis1.value)
        xnew1 = np.log10(NewAxis1.value)
    else:
        raise ValueError("Scale must be either lin or log.")
    if scale2=="lin":
        xold2 = OldAxis2.value
        xnew2 = NewAxis2.value
    elif scale2=="log":
        xold2 = np.log10(OldAxis2.value)
        xnew2 = np.log10(NewAxis2.value)
    else:
        raise ValueError("Scale must be either lin or log.")
    yold = OldValues.value
    
    # Perform interpolation
    interp = interpolate.RegularGridInterpolator((xold2, xold1), yold.T, method=method)
    
    # Evaluate new values
    Y, X = np.meshgrid(xnew2, xnew1)
    ynew = interp((Y,X))
    
    # Add unit
    NewValues = ynew * OldValues.unit
    
    return NewValues