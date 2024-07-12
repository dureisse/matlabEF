# matlabEF
A prototype of a finite element code written in matlab/octave without need for toolbox

The goal is to test some computational mechanics algorithms using matrix analysis mainly with finite elements. This is therefore not indeed cost-efficient but provides kind of a toolbox for testing and feasibility purposes.
There are two different data structures: one coming from reading AVS UCD ASCII standards input files and inspired by the finite element academic-open code Cast3M (https://www-cast3m.cea.fr/), the other one being a translation into matrix-like data structure.
Therefore, there are routines used to read files and provide the data in the first structure, then routines for translating into the second structure.
After computing, there are routines to come back to the first data structure, and to write output files, either of the same kind of the input files, or partially info the open code gmsh format (https://gmsh.info/).

The repository contains 
matlabEF, matlabEF2 and matlabUtils: contain the subroutine files
doc: contains a pdf document file, discussing the data format used
d-tests and d-test2: debug and non-regression tests, a script LaunchTests.csh can be used to check all of them. Since only source code are provided, the input data files should be generated first by using Cast3M code on *.dgibi files to generate some of the *.inp required files, before running the *.m codes.
They have been used on a MacBook Pro, macOS Monterey 12.5.1 with Cast3M version 21.0.1 and tested both with matlab2024a update 4, and GNU Octave version 9.1.0

This was originally and mainly designed by David Dureisseix (ENS Cachan, then University Montpellier 2, now INSA Lyon).
