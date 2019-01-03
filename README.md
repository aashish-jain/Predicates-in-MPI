# Predicates-in-MPI
A header file for processing predicates in parallel using MPI


<p>This project implements efficient parallel C++/OpenMPI function <code>mpi_extract_if</code>. The function works as follows: given a sequence <code>X</code> of some objects of type <code>T</code>, and a predicate that returns <code>true or false</code>, if an object satisfies some condition, the function creates sequence <code>Y</code> with a copy of only these objects in <code>X</code> for which the predicate is true in all the processes running and then evenly distributes the number of elements with a possible difference of at max p-1 elements. 

