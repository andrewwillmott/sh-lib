sh-lib
======

Spherical/zonal harmonics library

Provides the usual utilities for projecting environments into SH coefficients, 
evaluating them, and performing cheap(ish) operations like symmetric rotation 
and mirroring. Also contains a set of ZH routines for constructing rotationally 
symmetric BRDFs that can then be rotated into the full 3D SH space, or convolved
with a set of SH coefficients. 

To build and run the test app:

    c++ SHLib.cpp SHLibTest.cpp -o shtest && ./shtest

Or add those files to your favourite IDE.

See http://www.andrewwillmott.com/tech-notes#TOC-Spherical-Harmonics for more 
info.
