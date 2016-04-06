********
PSTRUDEL
********

Calculates distances between unlabeled phylogenetic trees and alignments of
different sizes.

Dependencies
============

This project has the following internal dependencies:

    * platypus-phyloinformary
      https://github.com/jeetsukumaran/platypus-phyloinformary

    * colugo-utilities
      https://github.com/jeetsukumaran/colugo-utilities

If this project was downloaded as an archive, then these dependencies should
already be present in the working tree. If this project was checked out, then
you will have to initialize and checkout the dependencies as submodules:

    $ git submodule init
    $ git submodule update

Building
========

Standard System-Wide Install
----------------------------

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make
    $ sudo make install

Testing
-------

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make check

Options
-------

The following options can be passed to the invocation of `cmake`:

    -DCMAKE_INSTALL_PREFIX="</path/to/install/dir>"

        Sets the prefix (top directory) for the installation.

    -DNCL_PREFIX="</path/to/NCL/installation>"

        Use the NCL (NEXUS Class Library) installation available at the
        specified path instead of building and using the bundled version.

    -DNCL_DEFAULT="yes"

        Use the default system NCL (NEXUS Class Library) installation instead
        of building and using the bundled version.

    -DINSTALL_NCL="yes"

        Install the locally-built NCL to ${CMAKE_INSTALL_PREFIX} (not needed
        run any PSTRUDEL products: this is only if you want the NCL library
        available for your own use).

For example:

    mkdir -p build/release
    cd build/release
    cmake -DCMAKE_INSTALL_PREFIX=${PWD}/install ../..

Copyright and License
=====================

   Copyright 2013 Jeet Sukumaran.
   All rights reserved.

   With code contributions from: Mark T. Holder.
   See "LICENSE.txt" for terms and conditions of usage.
