# -*- coding: utf-8 -*-
"""
NAME
    calc

DESCRIPTION
    Calc
    ====
    
    Compute the time evolution of continous time quantum walks (CTQWs) of two
    photons of arbitrary distinguishability by calculating the correlation
    matrix - the matrix whose elements are the liklihood of detecting a two
    photon coincidence in a photon count detector.
    
    Provides:
    1. Tools to compute the correlation matrix of a 2 photon CTQW.
    2. Calculate action of arbitrary sequence of unitary evolutions on a CTQW
        a. throughout a set period of time
        b. for different degrees of distinguishability.
    3. Utility functions to check function inputs

FUNCTIONS
    can_evolve()
        
    check_square_array()
        
    check_unitary()
        
    evolve_state()
        
    form_generator_matrix()
        
    get_correlation_matrix()
        
    is_good_correlations()
        
    is_hermitian()
        
    is_unitary()
        
    is_valid_mode()
        
    is_valid_timestep()
        
    vary_theta_evolutions()
        
"""
# =============================================================================
# IMPORT OFFICIAL & SELF WRITTEN MODULES
# =============================================================================
import numpy as np
import pars # where all shared information is stored


# =============================================================================
# CHECKING FUNCTIONS
# =============================================================================
def check_square_array(obj):
    """Check a given object is a square array, raising Exception for the condiiton it violates."""
    # Check typing
    if not isinstance(obj, np.ndarray):
        raise Exception("array not unitary: array is not of type numpy.ndarray")
    # Check square
    if not len(obj.shape)==2:
        raise Exception("array not unitary: array is not 2 dimensional")
    elif not obj.shape[0]==obj.shape[1]:
        raise Exception("array not unitary: array is not square")

# Check unitarity 
def check_unitary(array):
    """Check an array fulfils the properties of a unitary matrix and raise exceptions relevant to the condition it violates."""
    check_square_array(array)
    # Check det(U)=1
    determinant = np.linalg.det(array) # rounding is causing errors?
    if determinant == 1:
        pass
    else:
        rounded = np.around(determinant, decimals = 5)
        absolute= np.around(abs(determinant),decimals=5)
        if rounded == 1:
            pass
        elif absolute == 1:
            pass
        else:
            raise Exception("array not unitary: det(array) is not unity")
    # Check U U+ = U+ U = I
    herm = array.T.conj()
    eye  = np.eye(array.shape[0])
    U    = array @ herm # np.matmul(array, herm)
    recover_eye = np.all( np.around(U,8) == eye)
    if not recover_eye:
        raise Exception("array not unitary: did not recover identity")

def is_unitary(array):
    """Check if an array is unitary and return bool"""
    try:
        check_unitary(array)
        return True
    except Exception:
        return False

def is_hermitian(array):
    if np.all(array==array.T.conj()):
        return True
    else:
        return False

# Check modes agaisnt dimensions of system
def is_valid_mode(mode_tuple,unitary):
    # Check there are no modes outside the dimensions of the unitary
    shape = unitary.shape[0]
    degree= max(mode_tuple)
    if not degree<=shape:
        return False
    else:
        return True

# Check timesteps is valid
def is_valid_timestep():
    """See if the input timestep is valid"""
    # modify for a list of durations!
    # from save we know all of durations are floats and timestep is float

    any_negative = any([i<0 for i in pars.duration]) or pars.timestep < 0
    if any_negative:
        return False
    if any([i<pars.timestep for i in pars.duration]):
        return False
    else:
        return True

# Combination checker.
def can_evolve(unitary):
    # find success condition - must pass all tests to proceed
    conditions = (is_valid_mode(pars.initial_mode,unitary),is_valid_timestep())
    if all(conditions):
        return True
    else:
        return False

def is_good_correlations(correlations):
    """Perform analysis on correlation matrix. return True if the matrix satisfies the following conditions: """
    # check symmetry
    is_symmetric = is_hermitian(correlations)
    # check sum upper triangle elements = 1.
    is_normalised = round(np.triu(correlations).sum(),10) == 1.0
    # Check all values are greater than 0
    is_positive = np.all(correlations>=0)
    if is_symmetric and is_normalised and is_positive:
        return True
    else:
        return False

def form_generator_matrix(adjacency):
    # Need to make each collumn and row normalised. - check this condition!
    sums = np.sum(adjacency, axis=0)
    # 1. Check not already normalised, if so return itself
    if np.all(sums==1):
        return adjacency
    # 2. If it isn't normalised but is symmetric, normalise everything
    elif np.all(adjacency==adjacency.T):
        for i,s in enumerate(sums):
            adjacency[:,i] *= 1/s
        return adjacency
    # 3. If it isn't normalised or symmetric, throw an error
    else:
        raise Exception("Unable to interpret adjacency as adjacency matrix of graph")
    pass


# =============================================================================
# FUNCTIONS TO GET CORRELATION MATRICES & DATA
# =============================================================================
# Function for creating correlation matrix from a starting state and unitary.
def get_correlation_matrix(state,unitary,distinguishability=0):
    """
    Calculate the correlation matrix for measurements taken after unitary is applied to initial mode state.

    Parameters
    ----------
    state : tuple
        DESCRIPTION.
    unitary : numpy.ndarray
        DESCRIPTION.
    distinguishability: float
        Varies 0-1, equal to the proportion of distinguishable particles 
        (sin^2 theta). The default is 0, meaning perfectly identical photons.

    Returns
    -------
    corr_mat : numpy.ndarray
        DESCRIPTION.

    """
    mode1=state[0]
    mode2=state[1]
    d = distinguishability
    delta = lambda x,y: 0.5 if x==y else 1
    dim = unitary.shape[0]
    corr_mat = np.zeros((dim,dim))
    # If a detector that CANNOT DISTINGUISH BETWEEN PHOTONS measures coincidences at every possible pair of modes (p,q)
    for p in range(dim):
        for q in range(dim):
            pair1 = unitary[p,mode1] * unitary[q,mode2]
            pair2 = unitary[p,mode2] * unitary[q,mode1]
            # pair2 = unitary[mode1,p] * unitary[mode2,q]
            # pair2 = unitary[mode2,p] * unitary[mode1,q]
            identical = abs(pair1+pair2)**2
            distinct  = abs(pair1)**2 + abs(pair2)**2
            # I think this is right
            corr_mat[p,q] = delta(p,q)* (identical*(1-d) +distinct*d)
    return corr_mat

def evolve_state(state, steps_tuple, unitary_tuple, theta=0):
    """
    Calculate the set of correlation matrices for each timestep t of a system 
    evolving via the 1st unitary in unitary_tuple for a number of steps equal
    to the 1st element in steps_tuple, then the ith unitary and step duration.
    

    Parameters
    ----------
    state : tuple
        Tuple containing the initial modes of the system.
    steps_tuple : tuple
        Tuple of integers denoting the number of steps the ith unitary be
        applied to the system for.
    unitary_tuple : tuple
        Tuple of numpy.ndarrays, each representing unitary matrices that evolve
        the system for one timestep.
    theta : float, optional
        The degree of distinguishability between photons in the evolution. The
        default is pars.theta.

    Raises
    ------
    TypeError
        DESCRIPTION.
    ValueError
        DESCRIPTION.
    Exception
        DESCRIPTION.

    Returns
    -------
    data : numpy.ndarray
        A TxNxN array, where T is the sum of all the elements in steps_tuple,
        and N is the dimension of each unitary array. The n-th slice of data
        (data[n]) is the correlation matrix of the system after n timesteps.

    """
    # for a given set of unitarys in unitary_tuple,each to be performed a certain number of times dicttated by steps_tuple, evolve the system accordingly.
    
    distinguishability = np.sin(theta)**2 # proportion of distinct photons
    
    # convert 1...N mode indexing to 0...N python indexing
    state = tuple([i-1 for i in state])
    
    # Error Checking:
    # 1. All items must be ints (1a) or ndarrays (1b)
    # 2. All ndarrays must represent unitary matrices
    # 3. All unitary matrices must represent the same system i.e. same dim.
    is_ints     = all([isinstance(t,(int,np.int32)) for t in steps_tuple])
    is_arrays   = all([isinstance(U,np.ndarray) for U in unitary_tuple])
    is_unitarys = all([is_unitary(U) for U in unitary_tuple])
    is_float    = isinstance(theta,(float, int))
    shape_list = [U.shape for U in unitary_tuple]
    all_same_shape = len(set(shape_list))==1
    
    if not is_ints:
        raise TypeError("evolve_state requires steps_tuple to be tuple of integers. ")
    if not is_arrays:
        raise TypeError("evolve_state requires unitary_tuple to be tuple of numpy.ndarrays.")
    if not is_unitarys:
        raise ValueError("evolve_state requires each numpy.ndarray in unitary_tuple to behave as a unitary matrix")
    if not all_same_shape:
        raise Exception("evolve_states requires each unitary matrix in unitary_tuple to have the same shape in order to act on the same system.")
    if not is_float:
        raise TypeError("evolve_state requires theta to be a float")
    
    # Generate empty dataset to store correlation matrices
    dim = (unitary_tuple[0]).shape[0]
    data = np.zeros( (sum(steps_tuple)+1, dim, dim), dtype=float)
    # Initialise 1st correlation matrix
    data[0,state[0],state[1]] = 1
    data[0,state[1],state[0]] = 1
    # Create an overall increment counter, t, to track whole system evolution
    t = 1
    
    # Make identity as 1st U for t=0 and form pairs (t_i, U_i).
    unitary = np.eye(dim)
    pairset = zip(steps_tuple, unitary_tuple)
    
    # For each i in the tuples, evolve the system w., U_i for t_i steps
    for pair in pairset:
        ti = pair[0]
        Ui = pair[1]
        for step in range(ti): #i.e. do ti times
            unitary  = Ui @ unitary #np.matmul(Ui,unitary)
            # unitary  = unitary @ Ui
            corr_mat = get_correlation_matrix(state, unitary, distinguishability)
            data[t]  = corr_mat
            t+=1 #increment total number of steps so far
    return data

def vary_theta_evolutions(steps_tuple,unitary_tuple,theta_tuple):
    # perform evolutions for N many states after n steps 
    
    # Error Checking:
    # 1. All items must be ints (1a) or ndarrays (1b)
    # 2. All ndarrays must represent unitary matrices
    # 3. All unitary matrices must represent the same system i.e. same dim.
    is_ints     = all([isinstance(t,int) for t in steps_tuple])
    is_arrays   = all([isinstance(U,np.ndarray) for U in unitary_tuple])
    is_unitarys = all([is_unitary(U) for U in unitary_tuple])
    is_floats   = all([isinstance(theta,float) for theta in theta_tuple])
    shape_list = [U.shape for U in unitary_tuple]
    all_same_shape = len(set(shape_list))==1
    
    if not is_ints:
        raise TypeError("evolve_state requires steps_tuple to be tuple of integers. ")
    if not is_arrays:
        raise TypeError("evolve_state requires unitary_tuple to be tuple of numpy.ndarrays.")
    if not is_unitarys:
        raise ValueError("evolve_state requires each numpy.ndarray in unitary_tuple to behave as a unitary matrix")
    if not all_same_shape:
        raise Exception("evolve_states requires each unitary matrix in unitary_tuple to have the same shape in order to act on the same system.")
    if not is_floats:
        raise TypeError("evolve_state requires theta_tuple to be a tuple of floats")
    
    # Create dataset to store correlation matrices
    dim = (unitary_tuple[0]).shape[0]
    data=np.zeros((len(theta_tuple),dim,dim))
    # Construct the unitary for the specified system
    unitary = np.eye(dim)
    for pair in zip(steps_tuple, unitary_tuple):
         unitary = np.linalg.matrix_power(pair[1], pair[0]) @ unitary
    # Find the correlation matrix for the system for every value of theta given
    distinguishabilites = np.sin(np.array(theta_tuple))**2
    for i,d in enumerate(distinguishabilites):
        data[i] = get_correlation_matrix(pars.initial_mode, unitary, d)
    return data