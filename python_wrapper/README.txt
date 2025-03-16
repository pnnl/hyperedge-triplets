Here we explain how to compile the codes.

If you use this code, cite the following paper:

Niu, J., Amburg, I. D., Aksoy, S. G., & Sarıyüce, A. E. (2024, December). Retrieving Top-k Hyperedge Triplets: Models and Applications. In 2024 IEEE International Conference on Big Data (BigData) (pp. 630-639). IEEE.

Portable memory mapping class from https://github.com/stbrumme/portable-memory-mapping.

Parallel hashmap from https://github.com/greg7mdp/parallel-hashmap.

This code contains a Python wrapper built on top of C++ and was tested using Python 3.8.10 and GCC 7.3.0.
Tested on Linux and Windows operating systems.

Compilation
-----------

To use the Python wrapper, we first need to compile the C++ code into a shared library as follows:

    1) g++ -fPIC -shared -o max_triplet.so max_triplet.cpp portable-memory-mapping-master/MemoryMapped.cpp -std=c++17 -O3

Make sure the shared library max_triplet.so is in the same directory as your code which imports the python wrapper.
An example of how to run the code is in max_triplet.ipynb.

Dataset format:
    |E| |U| |V|
    u1 v1
    u2 v1
    u2 v2

Example dataset:
    3 2 2
    0 0
    0 1
    1 0
