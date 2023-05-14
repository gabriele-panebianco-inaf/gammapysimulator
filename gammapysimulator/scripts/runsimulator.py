#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

from time import time
start = time()

import argparse

from gammapysimulator.simulator import simulator

def Simulate():
    
    # Monitor imports
    Imports_Time = time()-start
    
    # pipeline options
    parser = argparse.ArgumentParser()
    parser.add_argument("-conf", "--configurationfile", help="Name of the analysis configuration YAML file")
    args = parser.parse_args()

    # Instantiate Simulator
    sourcesimulator = simulator.Simulator(args.configurationfile)
        
    # Run Simulation
    datasets = sourcesimulator.RunSimulation()
    
    # Final Time
    sourcesimulator.log.info(f"Runtime Imports = {float(Imports_Time):.3f} s.")
    sourcesimulator.log.info(f"Total Runtime = {float(time()-start):.3f} s.")
    
    return None

if __name__ == '__main__':
    Simulate()
