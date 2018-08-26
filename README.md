DistLib
=======

Small library for distributing an underlying uniform uint32 distribution into various 
other 1D, 2D, and 3D distributions. Features:

* Generation of int, float, 2D and 3D samples over various shapes
* Cheap routines that only take one input number, but may be less accurate or controllable.
* Accurate routines that take a full number of input samples.
* Pre-warping of distribution into triangular, gaussian-like, and weighted forms in a
  way that is composable with

To build and run the test app:

    c++ --std=c++11 Splines.cpp SplinesTest.cpp -o splines && ./splines


Examples
--------

Circle and Disc

![points](images/circle.gif "Splines from Points")

