Steps to add max_triplet to HyperNetX:

(1) Clone HyperNetX from https://github.com/pnnl/HyperNetX.
    - We refer to this as HyperNetX-master

(2) Add max_triplet.py to HyperNetX-master/hypernetx/algorithms

(3) Add
        from .max_triplet import *
    to HyperNetX-master/hypernetx/algorithms/__init__.py

(4) Add test_max_triplet.py to HyperNetX-master/hypernetx/tests/algorithms

(5) Add Advanced 8 - Maximum Hyperedge Triplets.ipynb to HyperNetX-master/tutorials/advanced
    - Change Advanced 8 to a later number if necessary

(6) Add images/ShadedTriplet.png to HyperNetX-master/tutorials/images
