# Implementation of Ant Colony Optimisation for the Travelling Salesman Problem #
 
The algorithm is described in [1, page 8].

See AntTsp.java for more details.

## Usage ##

- Compile: javac AntTsp.java
- Run: java AntTsp <TSP file>

## TSP file format ##

Full adjacency matrix. Columns separated by spaces, rows by newline.
Weights parsed as doubles, must be greater than or equal to zero.

## References ##

[1] M. Dorigo, The Ant System: Optimization by a colony of cooperating agents -- ftp://iridia.ulb.ac.be/pub/mdorigo/journals/IJ.10-SMC96.pdf
