===============
Getting Started
===============

Download the *PowerParser* package, currently PowerParser_v2.0.7.tgz. The tests use
the genmalloc package, genmalloc_v2.0.7.tgz, to allocate contiguous 2D arrays for testing
the *PowerParser* code. So download that if you want to run the tests.

Untar the genmalloc and PowerParser files::

   tar -xzvf genmalloc_v2.0.7.tgz
   tar -xzvf PowerParser_v2.0.7.tgz

The build processs uses cmake and is similar for both packages. If building with genmalloc,
modify the CFLAGS and LDFLAGS variable to include the path to the genmalloc build
products::

   export CFLAGS="${CFLAGS} -I<genmalloc_path>/include
   export LDFLAGS="${LDFLAGS} -L<genmalloc_path> -lgenmalloc

Run the cmake command as follows::

   cmake .                         // in-tree build
   cmake <path-to-src>             // out-of-tree build
   cmake -DCMAKE_BUILD_TYPE=debug <path-to-src>
   cmake -DCMAKE_BUILD_TYPE=release <path-to-src>
   cmake -DCMAKE_INSTALL_PREFIX=/usr/local // set the install directory

Then build the package::

   make

And install::

   make install
