Kiva
====

Kiva is a ground heat transfer calculation tool written in C++. Specifically, Kiva is used
to calculate heat loss and gain on a timestep basis from building foundations. The goal is
to create a tool that can integrate the multi-dimensional heat transfer into standard building energy simulation engines.

Contributing
------------

Kiva is configured as an Eclipse project built using the GCC toolchain.

To build Kiva, you'll need to import the Git repository into your copy of Eclipse. I prefer using [eGit](http://www.eclipse.org/egit/) to interface with Git from within Eclipse. For the best development experience, I suggest you develop on an Ubuntu 64-bit machine.

Pre-requisites:

1. GCC
2. CMake (for library dependencies -- if there is a demand for a build system outside of Eclipse, I will consider creating a CMake build for Kiva).

Optional:

1. Eclipse
2. eGit


Kiva has a number of library dependencies which you will need in order to compile the code.

1. [Boost](http://www.boost.org/) (with compiled binaries)
2. [Boost Numeric Bindings](http://mathema.tician.de/software/boost-numeric-bindings)
3. [SuiteSparse](http://www.cise.ufl.edu/research/sparse/umfpack/) (with openBLAS set in SuiteSparse_config.mk)
 - [OpenBLAS](http://xianyi.github.io/OpenBLAS/)
4. [MathGL](http://mathgl.sourceforge.net/) (with `-D enable-gif=on` passed to cmake)
 - [giflib](http://sourceforge.net/projects/giflib/)
 - [libpng](http://www.libpng.org/pub/png/libpng.html)
     - [zlib](http://www.zlib.net/)
5. [yaml-cpp](https://code.google.com/p/yaml-cpp/)

If you'd like to contribute to this code or if you have questions, send an email to Neal 
Kruis (neal.kruis AT bigladdersoftware DOT com).
