PSTRUDEL
========

Calculates distances between unlabeled phylogenetic trees and alignments of
different sizes.

Dependencies
------------

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
--------

Standard System-Wide Install
----------------------------

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make
    $ sudo make install

Development Testing and Install
-------------------------------

    $ mkdir -p build/debug
    $ cd build/debug
    $ cmake .. -DCMAKE_INSTALL_PREFIX="$(PWD)/install"
    $ make check
    $ make install

Copyright and License
---------------------

   Copyright 2013 Jeet Sukumaran.
   All rights reserved.

   With code contributions from: Mark T. Holder.
   See "LICENSE.txt" for terms and conditions of usage.
