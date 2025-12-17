sh-lib
======

A library for creating and manipulating spherical and zonal harmonics.

Provides the usual utilities for projecting environments into SH coefficients,
evaluating them, and performing cheap(ish) operations such as symmetric rotation
and mirroring. Also contains a set of zonal harmonics routines for constructing
rotationally-symmetric BRDFs which can then be rotated into the full 3D
spherical harmonics space along a particular direction, or convolved with a set
of spherical harmonics.

The focus is on fast(ish) code for low band counts, so a number of routines rely
on pregenerated code, and only support up to a certain band count, usually at
least 5 (SH) or 7 (ZH). This is noted in the comments where it applies.

To build and run the test app:

    c++ ZHLib.cpp SHLib.cpp SHLibTest.cpp -o test && ./test

Or add those files to your favourite IDE.

There is also a makefile supplied -- run `make` to build release and debug
versions of the library, and `make test` for the test app. `make install` will
install the library and headers under /usr/local, or use `PREFIX=/my/path make
install` for other locations.

There is an [interactive demo](https://andrewwillmott.github.io/app/DemoApp.html?test=Spherical%20Harmonics&ui.hue=210)
of some of this functionality.

See the shaders directory for some of these routines recast in GLSL. There is
a corresponding [shader toy demo](https://www.shadertoy.com/view/dl3Bz8).

https://github.com/andrewwillmott/sh-lib/assets/4959327/fe11cffe-49fa-455c-9a6b-21de4d04e900
