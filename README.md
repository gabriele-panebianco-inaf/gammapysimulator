# gammapysimulator
Gamma-ray astrophysical source simulator: produces Event List, Spectra Maps and Lightcurves

# Instructions
In order to run the simulator you need:
1. A `YAML` configuration file for the simulations options.
2. A `YAML` file that contains the models for the astrophysical sources you want to simulate: The format must be compatible with `gammapy`. 
3. One or more IRFs `FITS` files that describe the Response functions of the instrument you want to simulate (Effective Area, Energy Dispersion Matrix, PSF for 3D simulations, Instrumental Background Rate Spectrum).
4. To have created the `gammapysimulator` virtual environment and have installed the `gammapysimulator` python package inside it.
See [Installation](./README.md#installation) section for details on how to do it.

Once you have everything ready, activate the environment and run the simulator with the [runsimulator.sh](./runsimulator.sh) script, giving the path to the configuration file as argument:

    (base)$ conda activate gammapysimulator
    (gammapysimulator)$ ./runsimulator.sh /PATH/TO/CONFIGURATION/FILE.yml

It is a shortcut to run the [runsimulator.py](./gammapysimulator/scripts/runsimulator.py) script, but you can also call it directly if you want as:

    (base)$ conda activate gammapysimulator
    (gammapysimulator)$ runsimulator.py -conf /PATH/TO/CONFIGURATION/FILE.yml
# Installation

### 1 - Environment creation
The first step is to create the virtual environment `gammapysimulator` for the simulator.

You can either do it with `conda`, `mamba` or `conda-lock`:
- [Conda](https://docs.conda.io/projects/conda/en/stable/index.html) is an open-source package management system and environment management system.
To install with conda use:
    ```
    conda env create -f environment.yml
    ```
- [Mamba](https://mamba.readthedocs.io/en/latest/index.html) is a reimplementation of the conda package manager in `C++` (*Gran linguaggio!*).
It is much faster than `conda` because it allows for parallel downloading of repository data and package files using multi-threading, it uses `libsolv` for much faster dependency solving and core parts of mamba are implemented in `C++`.
    ```
    mamba env create -f environment.yml
    ```
- [conda-lock](https://conda.github.io/conda-lock/) is a lightweight library that can be used to generate fully reproducible lock files for conda environments.
`conda-lock` is able to fix every sub-dependency and presolve the environment when creating its lockfile, so the conda solver is not invoked when installing the packages  and creating the environment.
    ```
    conda-lock install -n gammapysimulator conda-lock.yml
    ```
### 2 - Installation of the `gammapysimulator` package 
The first time you activate the environment you must install the `gammapysimulator` package inside it.
You can either do it in a *non-developer mode* (any update or modification to the package source code will not be implemented in the environment) or *developer mode* (updates are automatically installed):
- *Non-developer mode*: activate the environment and run [setup.sh](./setup.sh), it will call `pip install .` and run the code inside [setup.py](./setup.py).
    ```
    (base)$ conda activate gammapysimulator
    (gammapysimulator)$ ./setup.sh
    ```
- *Developer mode*: activate the environment and run [setup_dev.sh](./setup_dev.sh), it will call `pip install -e .` and run the code inside [setup.py](./setup.py) in *editable* mode.
    ```
    (base)$ conda activate gammapysimulator
    (gammapysimulator)$ ./setup_dev.sh
    ```
### 3 - Test installation
You can test the installation is okay by running unit tests with:
    
    $ ./start_test.sh
# Environment development
This is the procedure followed to create the environment files, recommended for any change or update to environment packages.

1. Update [environment.yml](./environment.yml) **by hand** to keep track of the required dependencies and their versions.
2. Update [requirements.txt](./requirements.txt) by copying all the dependencies listed in [environment.yml](./environment.yml) except for python and pip. Every `=` must be changed into `==`.
3. Update [setup.py](./setup.py) by specifying the same python version chosen in [environment.yml](./environment.yml).
4. Create a new version of [conda-lock.yml](./conda-lock.yml) and replace the existing one.
In order to do so the developer must have installed `conda-lock` in one of its local environments, then run:
    ```
    conda-lock -f environment.yml -p linux-64
    ```
    This will create the new lockfile.
5. Create a new version of [requirements.lock](./requirements.lock) and replace the existing one.
In order to do so the developer must have installed `pip-tools` in one of its local environments (it can be done with `pip install pip-tools`), then run:
    ```
    pip-compile --no-annotate --output-file requirements.lock requirements.txt
    ```
