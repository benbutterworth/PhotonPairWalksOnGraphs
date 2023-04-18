# PhotonPairWalksOnGraphs

## About
`pairphotonwalk` is a command line tool that performs a continuous time quantum walk (CTQW) of two partially distinguishable photons on a given undirected graph or unitary matrix. It calculates the coincidence probability, $\Gamma_{ij}$, of detecting one photon in mode $i$ and one photon in mode $j$ , at regular time intervals for a specified duration.

## Project Files Description
+  **`pairphotonwalk.py`**
	Command line tool for simulating CTQWs. A user given parameter file is parsed and the CTQW of two photons described therein is performed.
+ **`calc.py`**
	The core simulation module; contains functions for checking inputs to `pairphotonwalk` and for calculating correlation matrices for two photon quantum walks.
+ **`save.py`**
	A utility module providing functions for reading and writing output data from `pairphotonwalk.py`.
+ **`pars.py`**
	Contains the essential information `pairphotonwalk.py` and `save.py` need to run.
+ **`parameter-example.csv`**
	An example parameter file that will perform a CTQW of two identical particles injected into the edge modes of a line graph of degree 10 for 10s in intervals of 0.2s.
+ **`adjacency-example.csv`**
	An example adjacency matrix for the line graph of degree 10 referenced by `parameter-example.csv`.

## Getting Started
 >**Before you begin,** `pairphotonwalk` **requires the following Python modules to be installed:**
>+ [ `numpy`](https://pypi.org/project/numpy/)
>+ [`pandas`](https://pypi.org/project/pandas/)
>+ [`scipy`](https://pypi.org/project/scipy/) 
>+ [`networkx`](https://pypi.org/project/networkx/)

>These modules are all contained within the popular [Anaconda Distribution](https://www.anaconda.com/products/distribution).

### Installation 
You should download this project as a .zip file and unzip it in a directory of your choice.
When you wish to run `pairphotonwalk` you should do so while located in this directory.

### Creating a parameter file

The parameter file specifies all the conditions of the CTQW to be simulated and is provided as an input to `pairphotonwalk.py`. 

The parameter file must be a .csv file with row names matching those provided in the table below. Some fields (`duration`, `initial mode`, `input file path`) take multiple inputs - these must be separated by *spaces* not commas.

This table contains a complete list of the parameters the user controls; these *must* be specified within the parameter file.

| Parameter | Aspect Controlled |
| --- | --- |
| `theta` | The degree of distinguishability between the two photons; varies between fully identical photons ($\theta=0$) to completely distinct photons ($\theta=\pi/2$).|
| `time-step` | The time interval between generating correlation matrices for the CTQWs progression [seconds].|
| `duration` | The duration of the CTQW simulation [seconds]. |
| `initial mode` | The initial modes each photon is injected into.|
| `input format` | Whether the input matrix for the CTQW should be interpreted as an `adjacency` matrix  for a graph, or a `unitary` time evolution matrix.|
| `input file path` | The filepath leading to the CTQW input matrix.|
| `output file path` |The filepath and name for the simulation data. There should be no file already at this path. |
| `description` | A reference description for the user. |
 
 A complete example parameter file for performing a CTQW on a line graph with 10 nodes is provided in `parameter-example.csv`
> **Tip:** You should create a unique parameter file for each simulation you run, especially if it will produce novel data.

### Defining an input matrix
In the parameter file the user specifies a location of a input matrix/matrices - this controls the environment the CTQW evolves in.

This input matrix can have one of two representations:
+ `adjacency` : the input matrix is an adjacency matrix of an undirected graph (it is a real-valued and symmetric matrix, i.e. Hermitian).
+ `unitary` : the matrix is a time evolution operator of one step in the simulation (it is a unitary matrix).    

Additionally, the user can list multiple file locations in the `input file path` box and use the following option.
+ `multiple` : each given matrix  ($M_i$) is a unitary time evolution operator to be applied for $t_i$ steps, where $t_i$ is the $i^\text{th}$ duration divided by the time-step.

>The program is intended to simulate CTQWs on graphs - `unitary` functionality is added to allow the user to calculate the unitary operators for very large adjacency matrices independently, then use these results as inputs.

## Execution Instructions
Once a parameter file and input matrix have been produced, the simulation can be performed by calling `pairphotonwalk.py` from the command line, shown below:
```
$ python pairphotonwalk.py -v ~/test-parameters.csv
	Evolution of system was performed successfully from parameter file stored at Input File. 
	The data was successfully stored as Output file.

	Input File  : ~/test-parameters.csv
	Output File : ~/test-data.csv
	Runtime     : 0.67500s
$ 	
```

The user can use flags to have additional control over how `pairphotonwalk` excecutes.
+ `-h` for help running `pairphotonwalk`.
+ `-p` to print the last correlation matrix to the terminal.
+ `-v` to print runtime information upon completion.
+ `-n` to not save outputs to output file. **NOT RECCOMENDED.**

## Interpreting Outputs
`pairphotonwalk`  outputs a set of correlation matrices, $G=\qty{ \Gamma(t)  \ \forall \ t\in [0, D ] }$, where $\Gamma(t)$ is the correlation matrix associated with the $t^{\text{th}}$ time-step and $D$ is the duration of the simulation. The elements of $\Gamma(t)$ are the coincidence probabilities, $\Gamma_{ij}$, of detecting one photon in mode $i$ and one photon in mode $j$ . $G$ is saved as a .csv file, where the $n^{\text{th}}$  row is the flattened correlation matrix associated with the $n^{\text{th}}$ timestep,  $\Gamma(n \cdot$  `time step` $)$ .

## References 
[1] A. Peruzzo _et al._, ‘Quantum Walks of Correlated Photons’, _Science_, vol. 329, no. 5998, pp. 1500–1503, Sep. 2010, doi: [10.1126/science.1193515](https://doi.org/10.1126/science.1193515). \
[2] D. Aharonov, A. Ambainis, J. Kempe, and U. Vazirani, ‘Quantum walks on graphs’, in _Proceedings of the thirty-third annual ACM symposium on Theory of computing_, Hersonissos Greece: ACM, Jul. 2001, pp. 50–59. doi: [10.1145/380752.380758](https://doi.org/10.1145/380752.380758). \
[3] Z.-Y. J. Ou, _Multi-photon Quantum interference_. New York: Springer, 2007.

## Licensing:  [MIT](https://choosealicense.com/licenses/mit/)

##  Acknowledgements
![UOB logo](https://logos-download.com/wp-content/uploads/2016/12/University_of_Bristol_logo.png)

This program was produced as a part of the authors Physics MSci final year project on _"Simulating Continuous Time Quantum Walks on Complex Graphs"_ under the supervision of J. Matthews and J. Frazer at The University of Bristol.

Credit: Benjamin Butterworth, 2023
[![GitHub Badge](https://img.shields.io/badge/GitHub-100000?style=for-the-badge&logo=github&logoColor=white)](https://github.com/benbutterworth) 
[![LinkedIn Badge](https://img.shields.io/badge/LinkedIn-0077B5?style=for-the-badge&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/bszbutterworth/)
