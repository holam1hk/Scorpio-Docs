.. _ch:parallelization:

***************
Parallelization
***************

Introduction
============
Scorpio uses Message Passing Interface (MPI) (``MPICH`` library)

Questions
=========
oct-tree-block based mesh refinement.?
Athena used a grid structure like Berger & Colella (1989),?

Future Development
==================
Plan to adopt the hybrid MPI + OpenMP approach. OpenMP parallelization is shared-memory parallelization used within a node while MPI distributed-memory parallelization used for messages exchange between nodes.

From Athena++ wiki
OpenMP parallelization is not very scalable in general. Generally speaking, we recommend to use flat-MPI or hybrid parallelization with a few (2-4) threads per physical core, except when you are using only one shared-memory node on which full OpenMP parallelization works fine. For CPUs with rich cores, such as Intel Xeon processors, flat-MPI parallelization often gives the best performance. Although modern CPUs support HyperThreading which allows you to launch more than one thread per core, usually you get no performance gain. On the other hand, for architectures with simpler cores but with fast hardware threads, including Intel Xeon Phi (Knights Landing; KNL) and IBM Blue Gene/Q (possibly newer POWER8/9 as well but we have not tested them), it is probably the best to launch a few threads per core (usually 4 for KNL and BG/Q) to keep the cores busy. Also, some modern processors such as AMD EPYC internally consist of multiple chiplets. In this case, you need to set assign threads and processes according to the internal structure in order to achieve the best performance. Because OpenMP threads can share some data, especially the MeshBlock tree, it saves some memory. When you are running extremely large parallel simulations, this will be helpful. In short, hybrid parallelization is recommended only when one or more of the following conditions are met:

    You are using CPUs with fast hardware threads (e.g. KNL and BG/Q).
    You need to save memory foot print by sharing some resources.
    Your system vendor requires it.
