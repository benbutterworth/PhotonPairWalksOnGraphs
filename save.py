# -*- coding: utf-8 -*-
"""
SAVE MODULE
    THIS MODULE IS PURPOSE BUILD FOR SAVING DATA FROM THE PROGRAM CALC_WALK.PY.
    IT TAKES 
@author: ben
"""
import numpy as np
import networkx as nx
import pandas as pd
from time import gmtime, strftime
from math import sqrt

import pars

# Define a path that takes you to where I'd like
graphpath = "D:/02_University/04_MSc/00_Project/04-Data/graphinfo/"
datapath  = "D:/02_University/04_MSc/00_Project/04-Data/"
# Define a generic name
name = "test.txt"


def gen_time_prefix():
    """Generate an identifier prefix string containing date and time and return string"""
    pref = gmtime()
    form = '23'+strftime("%m%d-%H%M", pref)+'-'
    return form


def save_parameters(filename, parameter_dictionary):
    """Save a set of parameters given as dictionary to a file at path filename"""
    # Save all the parameters to a single file
    df = pd.Series(parameter_dictionary)
    df.to_csv(filename, header=False)

def load_parameters(filename):
    """Read a parameter csv file and create the set of known variables"""
    # make sure everything is saved as global variables
    global theta, duration, timestep, initial_mode, input_format, input_file, output_file, description
    # Read into DataFrame then convert to Series, then to dict w. pairs
    par_series = pd.read_csv(filename,index_col=0, header=None).squeeze("columns")
    par_dict = par_series.to_dict()
    
    # retrieve theta
    try:
        pars.theta = float(par_dict["theta"])
    except Exception:
        raise Exception(f"theta cannot be defined from parameter file {filename}")
    # retrieve timestep
    try:
        pars.timestep = float(par_dict["time-step"])
    except Exception:
        raise Exception(f"timestep cannot be defined from parameter file {filename}")
    # retrieve duration
    try:
        pars.duration = float(par_dict["duration"])
    except Exception:
        raise Exception(f"duration cannot be defined from parameter file {filename}")
    # retrieve initial mode
    try:
        temp = par_dict["initial mode"]
        pars.initial_mode = tuple([int(i) for i in temp.split(' ')])
    except Exception:
        raise Exception(f"initial_mode cannot be defined from parameter file {filename}")
    # retrieve input method
    try:
        pars.input_format = str(par_dict["input format"])
    except Exception:
        raise Exception(f"input_format cannot be defined from parameter file {filename}")
    # retrieve input file
    try:
        pars.input_file = str(par_dict["input file path"])
    except Exception:
        raise Exception(f"input_file cannot be defined from parameter file {filename}")
    # retrieve output file
    try:
        pars.output_file = str(par_dict["output file path"])
    except Exception:
        raise Exception(f"theta cannot be defined from parameter file {filename}")
    # retrieve description
    try:
        pars.description = str(par_dict["description"])
    except Exception:
        raise Exception(f"description cannot be defined from parameter file {filename}")


def save_graph_to_file(filename, graph):
    """Convert a given networkx graph object to an adjacency matrix and store in a .csv file at path filename"""
    # Convert a graph to its adjacency matrix representation and save to a .csv file
    adj_mat = nx.to_numpy_array(graph)
    size = adj_mat.shape[0]
    col_names = [str(i+1) for i in range(size)]
    df = pd.DataFrame(adj_mat, columns=col_names)
    df.to_csv(filename)

def load_graph_from_file(filename):
    """Read a .csv file as an adjacency matrix and generate networkx graph object from it and return networkx graph object"""
    # Try to read a csv adjacency matrix from a file and turn variable into the graph object for it.
    # Read input as DataFrame
    df = pd.read_csv(filename, index_col=0)
    # Change to numpy array
    adj_mat = df.to_numpy()
    # Create graph object from it
    graph = nx.from_numpy_array(adj_mat)
    # Asign to variable
    return graph


def save_array_to_csv(filename, array):
    df = pd.DataFrame(array)
    df.to_csv(filename)

def load_array_from_csv(filename):
    # Read input as DataFrame
    df = pd.read_csv(filename, index_col=0)
    is_complex = np.any(df.dtypes == object)
    # if there are any complex objects make sure to convert back!
    if is_complex:
        data = df.to_numpy(dtype=complex)
    else:
        data = df.to_numpy()
    return data


def save_data_to_csv(filename,data):
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
    new_data.to_csv(filename)

def load_data_from_csv(filename,data_shape=None):
    
    # Read input as DataFrame & convert to numpy array
    df = pd.read_csv(filename, index_col=0)
    temp_data = df.to_numpy()
    # Gather information on number of rows and amount of items in each row
    length = temp_data.shape[0]
    size = temp_data.shape[1]
    
    # Try to reshape each row to a given shape, assuming data is square if none given.
    if data_shape == None:
        # See if the data can be put into square format.
        dim = int(sqrt(size))
        data = np.zeros((length,dim,dim))
        if not dim**2==size:
            raise Exception(f"data from {filename} cannot be reshaped into {length} many {(dim,dim)} arrays: specify data_shape.")
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
            raise Exception(f"data from {filename} cannot be reshaped into {length} many {data_shape} arrays")
    return data