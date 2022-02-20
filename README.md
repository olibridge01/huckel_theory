# Part II Chemistry Programming Course
Oli Bridge, <ob344@cam.ac.uk>

St Catharine's College, Cambridge

## Exercise 1: A general Hückel solver

In this exercise, a program was written that calculates the energies and degeneracies for various molecules of sp2 hybridised systems - linear molecules, cyclic molecules and the platonic solids. It does this by generating the Hückel Hamiltonian matrix for a given molecular system, diagonalising it to find the energies of the molecular orbitals, and tallying up the repeated eigenvalues to determine the degeneracies.

All energies in the program are computed given that α = 0 (energy of the AOs) and β = -1 (interaction energy between adjacent AOs).


## Libraries Required
The program runs on Python3 with the following libraries:
- numpy
- scipy
- sys
- networkx

## User Input
To run the program, execute the following command in the directory containing the python file:
```
python3 huckel.py [linear/cyclic/platonic] [number of atoms]
```
NB. n = 4, 6, 8, 12 or 20 for the platonic solids.

## Test Example
```
python3 huckel.py linear 4

Energy    Degeneracy
---------------------
 1.618         1         
 0.618         1         
-0.618         1         
-1.618         1         
---------------------
System has 4 orbitals
```
```
python3 huckel.py cyclic 6

 Energy    Degeneracy
---------------------
 2.000         1         
 1.000         2         
-1.000         2         
-2.000         1         
---------------------
System has 6 orbitals
```
## Notes 
The program extensively uses the numpy library for creation of matrix arrays and using linear algebra to find eigenvalues etc. 

Additionally, due to the close link between the connectivity of molecules with graph theory, I used the networkx library, which is often used for network analysis. In the case of the platonic solids, their atomic connectivities can be described very neatly by a graph/network, and their adjacency matrix (matrix where the ijth component is 1 if atom i is connected to atom j) is the hamiltonian matrix multiplied by -1.
