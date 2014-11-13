===============
Getting Started
===============

Download the *EZCL* package, currently ezcl_v2.0.7.tgz. Untar the EZCL files::

   tar -xzvf ezcl_v2.0.7.tgz

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

