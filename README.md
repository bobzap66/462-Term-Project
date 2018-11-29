# 462-Term-Project
Written by:

Jonathan Walker
Nick Greene


part1.c is a single processor matrix-matrix multiplication program. compile with mpicc and run with mpirun -n 1

part2.c is a multi-processor matrix-matrix multiplication program. compile with mpicc and run with mpirun -n <4, 16, 64, or 256>

part3.cpp is a multi-processor matrix-matrix multiplication program. compile with mpicxx and run with mpirun -n <4, 16, 64, or 256>
part3.cpp is having errors running with 64 or more processors. Not sure what the issue is.

Also included is "Run Time Part 2".png. It's a picture of the plotting of part1.c and part2.c running against 1, 4, 16, 64, and 256 processors. It was run a dozen times on each, and the times were averaged to produce the graph.

Finally, also included is Term Project Part 2.xlsx, which contains the data used to create the graph as well as a copy of the grpah itself.

