# -*- coding: utf-8 -*-
"""
NAME
    save

DESCRIPTION
    Save
    ====
    
    Bespoke I/O module for data used in an generated from pairphotonwalk.
    
    Provides:
    1. Loading simulation parameters from parameter file into pars module
    2. Loading and saving data generated from pairphotonwalk
    3. Utility for numpy array and networkx Graph file I/O
    
FUNCTIONS
    
    gen_time_prefix()
        Return string of date-time information for use as file prefix.

    load_array_from_csv()
        Recover numpy array from csv file without collumn labels.

    load_data_from_csv()
        Load data from a csv with T rows and M collumns into a numpy array with
        dimensions (T, X1, ... ,Xn) where M = X1*...*Xn.

    load_graph_from_file(filename)
        Recover networkx Graph from adjacency matrix stored in file.

    load_parameters(filename)
        Parse contents of a parameter file and change simulation variables in 
        pars module.

    save_array_to_csv(filename, array)
        Save an array to csv using Pandas Dataframe default options.

    save_data_to_csv()
        Save output data from pairphotonwalk into csv where each row is a
        flattened 2D correlation matrix, and the index row is mode pairs.

    save_graph_to_file(filename, graph)
        Save the adjacency matrix of a networkx Graph to a file.

    save_parameters(filename, parameters)
        Save a dictionary of parameters to a parameter file.

"""
import numpy as np
import networkx as nx
import pandas as pd
import time
import math

import pars


def gen_time_prefix():
    """
    Return string of date-time information for use as file prefix.
    
    Example
    >>> gen_time_prefix()
    '230406-1531-'
    
    >>> gen_time_prefix()+'data.csv'
    '230406-1532-data.csv'
    
    """
    pref = time.gmtime()
    form = time.strftime("%y%m%d-%H%M", pref)+'-'
    return form


def save_parameters(filename, parameters):
    """
    Save a dictionary of parameters to a .csv parameter file.
    
    Parameters
    ----------
    filepath : str
        Path to folder to save parameter file in.
    filename : str
        Name of the parameter file (with extension).
    parameters : dict
        Dictionary with keys as variable names and values of the respective 
        variables 
        

    Returns
    -------
    None
    
    """
    # Save all the parameters to a single file
    df = pd.Series(parameters)
    df.to_csv(filename, header=False)


def load_parameters(path):
    """
    Parse contents of a parameter file and change simulation variables in 
    pars module.
    
    Parameters
    ----------
    path : str
        Path to the parameter file to be read.

    Returns
    -------
    None
    
    Example:
    >>> load_parameters('~\identical_photon_parameters.csv')
    >>> pars.theta
    0.0
    >>> load_parameters('~\distinct_photon_parameters.csv')
    >>> pars.theta / (0.5*pi)
    1.0
    """
    # make sure everything is saved as global variables
    global theta, duration, timestep, initial_mode, input_format, input_file, output_path, output_name, description
    # Read into DataFrame then convert to Series, then to dict w. pairs
    par_series = pd.read_csv(path,index_col=0,
                             header=None).squeeze("columns")
    par_dict = par_series.to_dict()
    
    # retrieve theta
    try:
        pars.theta = float(par_dict["theta"])
    except Exception:
        raise Exception(f"theta cannot be defined from parameter file {path}")
    # retrieve timestep
    try:
        pars.timestep = float(par_dict["time-step"])
    except Exception:
        raise Exception(f"timestep cannot be defined from parameter file {path}")
    # retrieve duration
    try:
        temp = par_dict["duration"]
        pars.duration = tuple([float(i) for i in temp.split(' ')])
    except Exception:
        raise Exception(f"duration cannot be defined from parameter file {path}")
    # retrieve initial mode
    try:
        temp = par_dict["initial mode"]
        pars.initial_mode = tuple([int(i) for i in temp.split(' ')])
    except Exception:
        raise Exception(f"initial_mode cannot be defined from parameter file {path}")
    # retrieve input method
    try:
        pars.input_format = str(par_dict["input format"])
    except Exception:
        raise Exception(f"input_format cannot be defined from parameter file {path}")
    # retrieve input file
    try:
        temp = par_dict["input file path"]
        pars.input_files = tuple([str(i) for i in temp.split(' ')])
    except Exception:
        raise Exception(f"input_file cannot be defined from parameter file {path}")
    # retrieve output file
    try:
        pars.output_path = str(par_dict["output file path"])
    except Exception:
        raise Exception(f"output_path cannot be defined from parameter file {path}")
    # retrieve output path
    try:
        pars.output_name = str(par_dict["output file name"])
    except Exception:
        raise Exception(f"output_name cannot be defined from parameter file {path}")
    # retrieve description
    try:
        pars.description = str(par_dict["description"])
    except Exception:
        raise Exception(f"description cannot be defined from parameter file {path}")
        
    # Check lengths of tuples are the same
    if not len(pars.initial_mode)==2:
        raise Exception('Incorrect number of initial modes')
    lengths = [len(i) for i in (pars.duration, pars.input_files)]
    if not len(set(lengths))==1:
        raise Exception('Number of input matrices and durations does not match')


def save_graph_to_file(path, name, graph):
    """
    Save the adjacency matrix of a networkx Graph to a file.

    Parameters
    ----------
    path : str
        Path to the folder to store adjacency matrix
    name : str
        Name of file
    graph : networkx.classes.graph.Graph
        Networkx Graph object to be stored

    Returns
    -------
    None
    """
    # Convert a graph to its adjacency matrix representation and save to a .csv file
    adj_mat = nx.to_numpy_array(graph)
    size = adj_mat.shape[0]
    col_names = [str(i+1) for i in range(size)]
    df = pd.DataFrame(adj_mat, columns=col_names)
    df.to_csv(path+name)

def load_graph_from_file(path):
    """
    Recover networkx Graph from adjacency matrix stored in file.

    Parameters
    ----------
    path : str
        Path to the file containing adjacency matrix of a graph.

    Returns
    -------
    graph : networkx.classes.graph.Graph
        Graph recovered from adjacency matrix in file.
    """
    # Try to read a csv adjacency matrix from a file and turn variable into the graph object for it.
    # Read input as DataFrame
    df = pd.read_csv(path, index_col=0)
    # Change to numpy array
    adj_mat = df.to_numpy()
    # Create graph object from it
    graph = nx.from_numpy_array(adj_mat)
    # Asign to variable
    return graph


def save_array_to_csv(path, name, array):
    """
    Save an array to csv using Pandas Dataframe default options.
    
    Parameters
    ----------
    path : str
        Path to the folder to store file containing array.
    name : str
        Name of file
    array : numpy.ndarray
        Array to be stored in file
    
    Returns
    -------
    None
    """
    df = pd.DataFrame(array)
    df.to_csv(path, name)

def load_array_from_csv(filename):
    """
    Recover numpy array from csv file without collumn labels.

    Parameters
    ----------
    filename : str
        complete path to the file containing the array.

    Returns
    -------
    data : np.ndarray
        Data recovered from file.

    """
    # Read input as DataFrame
    df = pd.read_csv(filename, index_col=0)
    is_complex = np.any(df.dtypes == object)
    # if there are any complex objects make sure to convert back!
    if is_complex:
        data = df.to_numpy(dtype=complex)
    else:
        data = df.to_numpy()
    return data


def save_data_to_csv(path, name ,data):
    """
    Save output data from pairphotonwalk into csv where each row is a flattened
    2D correlation matrix, and index row is mode pairs.

    Parameters
    ----------
    path : str
        Path to folder to save file containing data in.
    name : str
        Name of the file to have data saved in.
    data : numpy.ndarray
        A TxNxN array where the slice data[i] represents the ith NxN
        correlation matrix of a evolving two photon system.

    Returns
    -------
    None.

    """
    # Get dimensions of unravelled which holds same info as data in fewer dims
    length = data.shape[0]
    size = np.product(data.shape[1:])
    # Unravel all data in each slice to 1D then place into unravelled
    unraveled = np.zeros((length,size),dtype=data.dtype)
    for i in range(length):
        unraveled[i] = data[i].ravel()
    # Make collumn names represent posn in ONE SLICE of the matrix
    col_names = []
    for index, info in np.ndenumerate(data[0]):
        col_names.append(str(index))
    # Save to .csv file
    new_data = pd.DataFrame(unraveled,columns=col_names)
    new_data.to_csv(path+name)

def load_data_from_csv(filename, data_shape=None):
    """
    Load data from a csv with T rows and M collumns into a numpy array with
    dimensions (T, X1, ... ,Xn) where M = X1*...*Xn.

    Parameters
    ----------
    filename : str
        Complete path to the file containing the output data from
        pairphotonwalk.
    data_shape : tuple, optional
        Shape of the matrix slice data[i] = (X1,...,Xn). The default is None,
        which will result matrix slice with dimensions (√M, √M).

    Returns
    -------
    data : numpy.ndarray
        Array containing T many correlation matrices.

    """
    # define own exception
    NoReshape= lambda dim, length: Exception(f"data from {filename} cannot be reshaped into {length} many {(dim,dim)} arrays: specify data_shape.")
    
    # Read input as DataFrame & convert to numpy array
    df = pd.read_csv(filename, index_col=0)
    temp_data = df.to_numpy()
    # Gather information on number of rows and amount of items in each row
    length = temp_data.shape[0]
    size = temp_data.shape[1]
    
    # Try to reshape each row to a given shape, assuming data is square if none given.
    if data_shape == None:
        # See if the data can be put into square format.
        dim = int(math.sqrt(size))
        data = np.zeros((length,dim,dim))
        if not dim**2==size:
            raise NoReshape(dim,length)
        # Reshape data into square format.
        for i in range(length):
            data[i] = temp_data[i].reshape((dim,dim))
    
    else:
        # Proceed only if reshaping into data_shape is possible
        is_right_size = np.product(data_shape)==size
        if is_right_size:
            data = np.zeros([length]+list(data_shape))
            for i in range(length):
                data[i] = temp_data[i].reshape(data_shape)
        else:
            raise NoReshape(dim,length)
    return data