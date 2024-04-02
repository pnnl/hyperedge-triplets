# Hyperedge Triplets
This repository contains implementations of the top-k hyperedge triplet retrieval algorithm from this [paper](https://arxiv.org/pdf/2311.07783.pdf).

If you use this code, please cite

@article{niu2023size,
  title={Size-Aware Hypergraph Motifs},
  author={Niu, Jason and Amburg, Ilya D and Aksoy, Sinan G and Sar{\i}y{\"u}ce, Ahmet Erdem},
  journal={arXiv preprint arXiv:2311.07783},
  year={2023}
}

There are three implementations: (1) HyperNetX, (2) standalone Python, and (3) C++ with a Python wrapper. Each implementation has a Jupyter notebook tutorial. 

## HyperNetX

The HyperNetX folder contains code that will be embeded in the hypergraph analytics library [HyperNetX](https://pnnl.github.io/HyperNetX/) (hnx). Steps to add max_triplet to HyperNetX:

(1) Clone HyperNetX from https://github.com/pnnl/HyperNetX.
    - We refer to this as HyperNetX-master

(2) Add max_triplet.py to HyperNetX-master/hypernetx/algorithms

(3) Add
        from .max_triplet import *
    to HyperNetX-master/hypernetx/algorithms/__init__.py

(4) Add test_max_triplet.py to HyperNetX-master/hypernetx/algorithms/tests

(5) Add Tutorial 14 - Maximum Hyperedge Triplets.ipynb to HyperNetX-master/tutorials
    - Change Tutorial 14 to a later number if necessary

(6) Add images/ShadedTriplet.png to HyperNetX-master/tutorials/images

## Standalone in Python

The python_standalone folder contains code that you could run independently of HyperNetX. 

## C++ with Python Wrapper

The python_wrapper folder contains code in C++ with a Python wrapper. This is the fastest implementation. 

Portable memory mapping class from https://github.com/stbrumme/portable-memory-mapping.

Parallel hashmap from https://github.com/greg7mdp/parallel-hashmap.

This code contains a Python wrapper built on top of C++ and was tested using Python 3.8.10 and GCC 7.3.0.
Tested on Linux and Windows operating systems.

### Compilation

To use the Python wrapper, we first need to compile the C++ code into a shared library as follows:

    1) g++ -fPIC -shared -o max_triplet.so max_triplet.cpp portable-memory-mapping-master/MemoryMapped.cpp -std=c++17 -O3

Make sure the shared library max_triplet.so is in the same directory as your code which imports the python wrapper.
An example of how to run the code is in max_triplet.ipynb.

Dataset format: <br />
    |E| |U| |V| <br />
    u1 v1 <br />
    u2 v1 <br />
    u2 v2 <br />

Example dataset: <br />
    3 2 2 <br />
    0 0 <br />
    0 1 <br />
    1 0




