PSTRUDEL
========

Calculates distances between unlabeled phylogenetic trees and alignments of
different sizes.

Dependencies
------------

Internal Dependencies
.....................

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

External Dependencies
.....................

This project has the following external dependencies:

    * The NEXUS Class Library (NCL)
      http://sourceforge.net/projects/ncl/

You can download and install these yourself separately, or else run the helper
scripts below (e.g., see ``Developer Setup and Installation``).

Developer Setup and Installation
--------------------------------

The following steps will download, build and (locally) install all
dependencies, as well as build and install the project in the local source
tree.

If this source code was cloned from a repository, then make sure to run:

    $ git submodule init
    $ git submodule update

Invoke the following script from the top-level directory to download, build,
and install the external dependencies in the local source tree:

    $ ./install_prerequisites.sh

Run the following script in the top-level directory set up the autotools build
framework:

    $ ./bootstrap.sh

Create a build directory and configure the build framework (assuming
'install_prerequisites.sh' was run):

    $ mkdir -p build/debug
    $ cd build/debug
    $ ../../cfgcommand.sh debug

Or, for release/production build:

    $ mkdir -p build/release
    $ cd build/release
    $ ../../cfgcommand.sh release

Build and install the project in the build directory (assuming 'install_prerequisites.sh'
was run, and in build subdirectory, e.g. 'build/release'):

    $ make
    $ make install

Following this recipe will result in all products being installed in
subdirectory 'installed/' of the build directory.

To run the tests:

    $ source ../../pstrudel_deps_env.sh
    $ installed/opt/pstrudel/test/run-tests.py

If you get an error message regarding missing libraries during the build or
run, just make sure that the environmental variables are correctly set:

    $ source ../../pstrudel_deps_env.sh

Copyright and License
---------------------

   Copyright 2013 Jeet Sukumaran.
   All rights reserved.

   With code contributions from: Mark T. Holder.
   See "LICENSE.txt" for terms and conditions of usage.
