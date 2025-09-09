.. _ch:getting_started:

***************
Getting Started
***************

Requirements
===============
Suggested packages:

- hdf5-1.10.9
- fftw-3.3.6
- mpich-3.2

Getting Started
===============

.. note:: Set up ``Makefile`` and ``.bashrc`` properly. Make sure you have libraries ``fftw3``, ``hdf5``, and ``mpich`` installed and know the path.


Quick Start
==========
Move to the code directory,

#. Clean the code

   .. code-block:: bash

      make clean

#. Compile the code

   .. code-block:: bash

      make

#. Run the code

   .. code-block:: bash

      mpiexec -n 4 ./Scorpio

   4 means 4 processors (cpu cores).

.. warning:: Don't use too many cpu cores!!


- Suggested Compiler: F90=mpif90

- Makefile: choose optimization flag -O2 for the first time, -O3 for faster run



Create a new environment (optional)
========================
#. Initialize conda

   .. code-block:: bash    

      conda init

#. Activate conda base environment

   .. code-block:: bash    

      conda activate

   You are now in the base environment.

#. Create new environment by cloning base

   .. code-block:: bash    

      conda create -n my_new_env --clone base

#. Create new environment from scratch with Python 3

   .. code-block:: bash    

      conda create -n my_new_env python=3

#. Activate the new environment

   .. code-block:: bash    

      conda activate my_new_env

#. Install packages

   .. code-block:: bash    

      conda install <package_name>

#. Update conda and all packages

   .. code-block:: bash    

      conda update --all

For more information, visit: https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html


Platforms
==================
- tianhexy

- cluster2

- hbli-s2

- hbli-s1

- scorpio

- gpu2/mike

- paulsr

- stor2

- nas3
