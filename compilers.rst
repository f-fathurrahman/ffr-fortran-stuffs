Compilers
=========

Some useful Intel Fortran Compiler options
------------------------------------------

The following options of ``ifort`` may be useful:

- ``-warn -nogen-interfaces``

- ``-check bounds`` or ``-CB``

- ``-check uninit``

- ``-openmp``

- ``-fpp``: call the Fortran preprocessor

- ``-mkl``: automatic linking to Intel MKL

- ``-static-intel``

- ``-xHost``


Some useful GNU Fortran Compiler options
----------------------------------------

- ``-Wall``

- ``-x f95-cpp-input``


G95 error: ``ld: cannot find crt1.o``
-------------------------------------


I found the following error when trying to compile with ``g95``

.. code-block:: shell-session

  $ ~/mysoftwares/g95-install/bin/x86_64-unknown-linux-gnu-g95 hello.f90
  ld: cannot find crt1.o: No such file or directory
  ld: cannot find crti.o: No such file or directory

Solution::

  export LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LIBRARY_PATH

My computer details

.. code-block:: shell-session

  $ uname -a
  Linux kusanone2 3.8.0-35-generic #50-Ubuntu SMP Tue Dec 3 01:24:59 UTC 2013
  x86_64 x86_64 x86_64 GNU/Linux

  $ gcc -v
  Using built-in specs.
  COLLECT_GCC=gcc
  COLLECT_LTO_WRAPPER=/usr/lib/gcc/x86_64-linux-gnu/4.7/lto-wrapper
  Target: x86_64-linux-gnu
  Thread model: posix
  gcc version 4.7.3 (Ubuntu/Linaro 4.7.3-1ubuntu1)

  $ ld -v
  GNU ld (GNU Binutils for Ubuntu) 2.23.2


