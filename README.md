![](docs/images/kiva-logo.png)

[![Documentation Status](https://readthedocs.org/projects/kiva/badge/?version=latest)](http://kiva.readthedocs.org/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/bigladder/kiva.svg?branch=develop)](https://travis-ci.org/bigladder/kiva)
[![Build Status](https://ci.appveyor.com/api/projects/status/pv2c4no2mv4uds26/branch/develop?svg=true)](https://ci.appveyor.com/project/bigladder/kiva/branch/develop) [![codecov](https://codecov.io/gh/bigladder/kiva/branch/develop/graph/badge.svg)](https://codecov.io/gh/bigladder/kiva)

Kiva
====

Kiva is a free and open source ground heat transfer calculation tool written in
C++. Specifically, Kiva is used to calculate heat loss and gain on a timestep
basis from building foundations. The goal is to create a tool that can integrate
the multi-dimensional heat transfer into standard building energy simulation
engines.

Documentation
-------------

See the [online documentation](http://kiva.readthedocs.org/en/latest/) for information on using Kiva and creating Kiva input files.

Contributing
------------

Kiva is configured as a cross-platform CMake project. To build Kiva, you'll need to clone the git repository and use CMake (pointing to the kiva root directory).

Pre-requisites:

1. A C++ compiler (e.g., Clang, GCC, MSVC)
2. CMake

Building Kiva from source
-------------------------

1. Clone the git repository, or download and extract the source code (`tar.gz` or `zip`).
2. Make a directory called `build` inside the top level of your source.
3. Open a console in the `build` directory.
4. Type `cmake ..`.
5. Type `cmake --build . --config Release`.
6. The Kiva executable (`kiva` or `kiva.exe`) will appear in your build directory.

If you'd like to contribute to this code or if you have questions, send an email to Neal
Kruis (neal.kruis AT bigladdersoftware DOT com).
