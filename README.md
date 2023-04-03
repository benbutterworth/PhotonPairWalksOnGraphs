# PhotonPairWalksOnGraphs

## About
`pairphotonwalk` is a command line tool that performs a continuous time quantum walk (CTQW) of two partially distinguishable photons on a given undirected graph or unitary matrix. It calculates the coincidence probability, $\Gamma_{ij}$, of detecting one photon in mode $i$ and one photon in mode $j$ , at regular time intervals for a specified duration.

## Project Files Description
+  **`calc.py`**
	Reads a user defined parameter file and simulates a CTQW of two photons as described therein.
+ **`save.py`**
	A utility module that provides functions for reading and writing output data from `calc.py`.
+ **`pars.py`**
	Contains the essential information `calc.py` and `save.py` need to run.
+ **`parameter-example.csv`**
	An example parameter file that will perform a CTQW of two identical particles injected into the edge modes of a line graph of degree 10 for 10s in intervals of 0.2s.
+ **`adjacency-example.csv`**
	An example adjacency matrix for the line graph of degree 10 referenced by `parameter-example.csv`.

## Getting Started
 >:warning: **Before you begin,** `pairphotonwalk` **requires the following Python modules to be installed:**
>+ [ `numpy`](https://pypi.org/project/numpy/)
>+ [`pandas`](https://pypi.org/project/pandas/)
>+ [`scipy`](https://pypi.org/project/scipy/) 
>+ [`networkx`](https://pypi.org/project/networkx/)

>These modules are all contained within the popular [Anaconda Distribution](https://www.anaconda.com/products/distribution).

### Installation 
Within pars.py you should 

### Creating a parameter file

The parameter file specifies all the conditions of the CTQW to be simulated and is provided as an input to `calc.py`. Below is a list of the parameters the user controls; these *must* be specified within the parameter file.

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
 
> :bulb: **Tip:** You should create a unique parameter file for each simulation you run, especially if it will produce novel data.

### Defining an input matrix
In the parameter file the user specifies a location of a input matrix - this controls the environment the CTQW evolves in.

This input matrix can have one of two representations:
+ `adjacency` : the input matrix is an adjacency matrix of an undirected graph (it is a real-valued and symmetric matrix, i.e. Hermitian).
+ `unitary` : the matrix is a tim evolution operator of one step in the simulation (it is a unitary matrix).

>The program is intended to simulate CTQWs on graphs - `unitary` functionality is added to allow the user to calculate the unitary operators for very large adjacency matrices and use these as inputs.

## Excecution Instructions
Once a parameter file and input matrix have been produced, the simulation can be performed by calling `calc.py` from the command line, shown below:
```
$ python pairphotonwalk.py ~/test-parameters.csv 
	Evolution of system was performed successfully from parameter file stored at Input File. 
	The data was successfully stored as Output file.

	Input File  : ~/test-parameters.csv
	Output File : ~/test-data.csv
	Runtime     : 0.67500s
$ 	
```

## Interpreting Outputs
`pairphotonwalk`  outputs a set of correlation matrices, $G= { \Gamma(t)  \ \forall \ t\in [0, D ] \}$, where $\Gamma(t)$ is the correlation matrix associated with the $t^{\text{th}}$ time-step and $D$ is the duration of the simulation. The elements of $\Gamma(t)$ are the coincidence probabilities, $\Gamma_{ij}$, of detecting one photon in mode $i$ and one photon in mode $j$ . $G$ is saved as a `.csv` file, where the $n^{\text{th}}$  row is the flattened correlation matrix associated with the $n^{\text{th}}$ timestep,  $\Gamma(n \cdot$  `time step` $)$ .

## :books: References 
[1] A. Peruzzo _et al._, ‘Quantum Walks of Correlated Photons’, _Science_, vol. 329, no. 5998, pp. 1500–1503, Sep. 2010, doi: [10.1126/science.1193515](https://doi.org/10.1126/science.1193515). \
[2] D. Aharonov, A. Ambainis, J. Kempe, and U. Vazirani, ‘Quantum walks on graphs’, in _Proceedings of the thirty-third annual ACM symposium on Theory of computing_, Hersonissos Greece: ACM, Jul. 2001, pp. 50–59. doi: [10.1145/380752.380758](https://doi.org/10.1145/380752.380758). \
[3] Z.-Y. J. Ou, _Multi-photon Quantum interference_. New York: Springer, 2007.

## Licensing:  [MIT](https://choosealicense.com/licenses/mit/)

##  :scroll: Acknowledgements
![enter image description here](https://logos-download.com/wp-content/uploads/2016/12/University_of_Bristol_logo.png)


This program was produced as a part of the authors Physics MSci final year project on _"Simulating Continuous Time Quantum Walks on Complex Graphs"_ under the supervision of J. Matthews and J. Frazer at The University of Bristol.

Credit: Benjamin Butterworth, 2023
[![GitHub Badge](https://img.shields.io/badge/GitHub-100000?style=for-the-badge&logo=github&logoColor=white)](https://github.com/benbutterworth) 
[![LinkedIn Badge](https://img.shields.io/badge/LinkedIn-0077B5?style=for-the-badge&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/bszbutterworth/)
