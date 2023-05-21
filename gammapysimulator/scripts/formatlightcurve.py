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

from gammapysimulator.tools import formatlightcurvemodel

def main():
    
    # Monitor imports
    Imports_Time = time()-start
    
    # pipeline options
    parser = argparse.ArgumentParser()
    parser.add_argument("-file"  , "--lightcurvefile", help="Name of the light curve file")
    parser.add_argument("-format", "--inputformat"   , help="Input light curve format")
    parser.add_argument("-ref"   , "--referencemjd"  , help="Reference time in MJD")
    parser.add_argument("-o"     , "--outputfile"    , help="Name of the output light curve file")
    args = parser.parse_args()

    # Instantiate Formatter: read light curve
    formatter = formatlightcurvemodel.FormatLightCurve(args.lightcurvefile, args.inputformat, args.referencemjd)
    
    # Write light curve in a gammapy-compatible format
    formatter.write(args.outputfile)
    
    # Final Time
    formatter.log.info(f"Runtime Imports = {float(Imports_Time):.3f} s.")
    formatter.log.info(f"Total Runtime = {float(time()-start):.3f} s.")
    
    return None

if __name__ == '__main__':
    main()
