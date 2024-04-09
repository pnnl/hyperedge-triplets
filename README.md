# Hyperedge Triplets
This repository contains implementations of the top-k hyperedge triplet retrieval algorithm.

Note: All author-related information has been removed from this temporary repository.
This repository should only be used for reproducibility in the double-blinded review process.
The final published paper will have the full version of this repository with complete copyright information.

## C++ with Python Wrapper

The python_wrapper folder contains code in C++ with a Python wrapper.

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




