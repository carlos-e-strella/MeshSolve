#  MeshSolve
A c rewrite of my python pde solver. Will be implemented as a python package with a new interface. 
I am not a mathematician nor an engineer, so this program is probably packed with errors.

## How does this program solve PDEs
I have implemented (barely) the finite element method, which is a numerical method to aproximate solutions for PDEs.
The method requires you to specify a mesh where your equations occurr and calculates a stiffness matrix for it.
The solution then involves solving a linear system involving said matrix and your initial conditions (and your boundary conditions).

## How does this program implement the finite element method?
This program provides functions to calculate the stiffness matrix for a specific mesh.
Once its finished, you will be able to declare your partial differential equation as a sum of different terms.
Each term has a weak form that is used to calculate its contribution to the stiffness matrix.
Finally, the program will calculate the resulting solution by solving a system of equations involving said matrix.