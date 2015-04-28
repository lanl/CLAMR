===============
Getting Started
===============

Download the *MallocPlus* package, currently MallocPlus_v2.0.7.tgz. Untar
the *MallocPlus* files::

   tar -xzvf MallocPlus_v2.0.7.tgz

The build processs uses cmake. Run the cmake command as follows::

   cmake .                         // in-tree build
   cmake <path-to-src>             // out-of-tree build
   cmake -DCMAKE_BUILD_TYPE=debug <path-to-src>
   cmake -DCMAKE_BUILD_TYPE=release <path-to-src>
   cmake -DCMAKE_INSTALL_PREFIX=/usr/local // set the install directory

Then build the package::

   make

And install::

   make install

To run the unit tests as a subpackage of an applications

   cd tests
   make MallocPlus_check

Running unit tests as a standalone package.

   make check

To build the documentation, in the main *MallocPlus* directory:

   make MallocPlus_doc

or for the standalone package:

   make doc


