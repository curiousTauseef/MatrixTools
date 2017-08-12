/////////////////////////////////////////////////////////////////////////////////////////////////
> Matrix Tools
/////////////////////////////////////////////////////////////////////////////////////////////////

This is a tiny generic matrix library for solving
different matrix tasks including linear equations, determinants, computing inverse matrices
as well as decompositions. The 2D matrices can either have a dynamic (runtime) imposed number of rows and columns 
or in order to make things a lot safer and catch many of the possible errors at compile time the number of rows and
columns can be directly imposed as template parameters. The matrices are also designed to be STL friendly
up to a certain point thanks to the cusom iterators that are included.

The problems are divided according to the specific tasks which need to be acomplished:
-solving linear equations (Gauss with partial pivoting implemented)
-computing inverse matrices (Gauss Jordan method implemented)
-LU decomposition (Crout and Doolittle algorithms)
-computing determinants (using either Gauss or LU decomposition<Crout/Doolittle>)

The plan is to add new algorithms for solving each problem and create a complete library which
contains the basic algorithms for performing matrix calculations. The library can be used without building
it, the .hpp files can be included right away for maximum flexibility. The actual declarations are
kept in the .hpp files and the definitions for the methods are provided in the inline files, just for the 
sake of keeping things nice and clean. For usage examples the library contains a .cpp file
which can be compiled into an executable just to see how things work from the start.