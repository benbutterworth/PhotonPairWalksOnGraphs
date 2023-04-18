# -*- coding: utf-8 -*-
"""
==============
pairphotonwalk
==============

pairphotonwalk is a command line tool that performs a continuous time quantum
    walk (CTQW) of two partially distinguishable photons on a given undirected
    graph or unitary matrix. It calculates a correlation matrix, Γ, at regular 
    time intervals for a specified duration, whose elements Γ(i,j) are the
    coincidence probabilities of detecting a two photon coincidence in modes i     
    and j.

pairphotonwalk takes .csv parameter files as input - a full example is provided      
    in parameter-example.csv. Information detailing the role of these 
    parameters in found in pars.py.

pairphotonwalk requires the scientific computing Python modules numpy, scipy, 
    pandas and networkx to run. These are available via PyPI or the standard 
    Anaconda Scientific Distribution.

"""
# import modules
import argparse
from time import perf_counter
from scipy.linalg import expm
#import all neccessary modules
from calc import *
from save import *


# Start timing
start_runtime = perf_counter()

# =============================================================================
# SETUP ARGUMENT PARSING FROM CMDLINE
# =============================================================================
parser =  argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)

parser.add_argument('parameters', 
                    help="path to parameter file for system")
parser.add_argument('-v','--verbose', 
                    action='store_true', 
                    help='print simulation information after runtime')
parser.add_argument('-p','--print', 
                    action='store_true',
                    help='print last correlation matrix.')
parser.add_argument('-n','--nosave', 
                    action='store_true',
                    help='do not save outputs to file. Reccomended to use with --print') 

args = parser.parse_args()


# =============================================================================
# GET USER PARAMETERS:
# theta, timestep, duration, initial_mode,input_format, input_file, output_file
# =============================================================================
load_parameters(args.parameters)
from pars import *

# =============================================================================
# READ IN ARRAY INFORMATION
# =============================================================================
if not input_format=='multiple':
    if input_format=='unitary' or input_format=='adjacency':
        array = load_array_from_csv(input_files[0])
    else:
        raise BadInputFormat
        
elif input_format=='multiple':
    # get all the big boys in.
    arrays = [load_array_from_csv(path) for path in input_files]


# =============================================================================
# PERFORM CHECKS ON DATA
# =============================================================================
if not input_format=='multiple':
    if can_evolve(array):
        pass
    else:
        raise Exception("Fatal Error Occured! Parameters for system evolution (timestep, duration, input_mode) do not match")
    
    check_square_array(array)

elif input_format=='multiple':
    for array in arrays:
        if can_evolve(array):
            pass
        else:
            raise Exception("Fatal Error Occured! Parameters for system evolution (timestep, duration, input_mode) do not match")
        
        check_square_array(array)

# =============================================================================
# GET THE UNITARY MATRIX & ENSURE VALID.
# =============================================================================

if input_format == 'adjacency':
    # Check that array will create a unitary
    if is_hermitian(array):
        start_unitarytime = perf_counter()
        unitary = expm(1j*array*timestep)
        end_unitarytime = perf_counter()
        unitarytime = end_unitarytime - start_unitarytime
        if unitarytime > 2:
            # ask if you'd like to save as a new unitary
            print(f"\nTIP - it took {unitarytime}s to calculate the unitary matrix for this quantum walk, would you like to save this matrix to refer to again(y) or continue to run the simulation without saving (any key)?\n")
            is_save = input() == 'y'
            if is_save:
                name = input("Enter filename for unitary matrix : ")
                save_array_to_csv(output_path+name, unitary)
                print(f'Unitary matrix saved as {output_path+name}')
            else:
                pass
    else:
        raise Exception("array does not represent adjacency matrix")

elif input_format == 'unitary':
    # Check that array is unitary, throw exception if not
    if is_unitary(array):
        unitary = array
    else:
        check_unitary(array)

elif input_format == 'multiple':
    # check unitarity of all unitaries
    for array in arrays:
        if is_unitary(array):
            unitary = array
        else:
            check_unitary(array)
    # if program progresses beyond this point, every matrix is unitary

else:
    raise BadInputFormat

# =============================================================================
# USER PARAMETER DERIVED SIMULATION PARAMETERS
# =============================================================================
steps_array = np.array(np.array(duration)/timestep,dtype=int)


# =============================================================================
# PERFORM WHOLE EVOLUTION OF SYSTEM
# =============================================================================

if input_format=='multiple':
    data = evolve_state(initial_mode,tuple(steps_array), tuple(arrays),theta=theta)
elif input_format=='unitary' or input_format=='adjacency':
    data = evolve_state(initial_mode, tuple(steps_array), (unitary,),theta=theta)
else:
    raise BadInputFormat


# Raise warnings where correlation matrix is not satisfactory
status_correlations = [is_good_correlations(corr_mat) for corr_mat in data]
loc_bad_correlations= [loc for loc,truthvalue in enumerate(status_correlations) if truthvalue==False]
for step in loc_bad_correlations:
    print(f"WARNING - Bad Correlation Matrix detected at timestep t={step}")


# =========================================================================
# Perform flag operations
# =========================================================================
if not args.nosave:
    save_data_to_csv(output_path,output_name, data)
else:
    pass

if args.verbose:

    end_runtime = perf_counter()
    runtime = end_runtime-start_runtime
    print('Simulation Information\n======================')
    print(f'pairphotonwalk ran successfully, performing the walk from {args.parameters} in {runtime}s')
    print(f'\nSimulation Description:\n======================\n{description}')
    
    print('\nOutput data\n======================')
    if not args.nosave:
        print(f"Save Location: {output_path+output_name}")
    else:
        print('Data not saved to file')
    
else:
    pass

if args.print:
    if data[-1].shape[0] > 20:
        really_print = input('Data is large: print anyway (y)? ')=='y'
        if really_print:
            print(data[-1])
        else:
            print('Data not printed.')
    else:
        print('======================')
        print(data[-1])
else:
    pass