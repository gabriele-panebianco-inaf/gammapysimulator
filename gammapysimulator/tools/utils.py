#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

import os
import shutil

from gammapy.utils.table import table_from_row_data

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