#  MeshSolve
A c rewrite of my python pde solver. Will be implemented as a python package with a new interface. 
I am not a mathematician nor an engineer, so this program is probably packed with errors.

## How does this program solve PDEs?
I have implemented (barely) the finite element method, which is a numerical method to aproximate solutions for PDEs.
The method requires you to specify a mesh where your equations occurr and calculates a stiffness matrix for it.
The solution then involves solving a linear system involving said matrix and your initial conditions (and your boundary conditions).

## How does this program implement the finite element method?
This program provides functions to calculate the stiffness matrix for a specific mesh.
Once its finished, you will be able to declare your partial differential equation as a sum of different terms.
Each term has a weak form that is used to calculate its contribution to the stiffness matrix.
Finally, the program will calculate the resulting solution by solving a system of equations involving said matrix.

## How do I use this program?
In the future I plan to implement this as a python package to let it work with pyvista, but right now the ``matrix_factory`` function handles global stiffness matrix construction with regular c arrays so you can call that.
Also, the cmake files *should* include openblas and lapack if you have vcpkg installed and I also provided the clapack fortran to c interface because I refuse to build lapacke from source. If you want yo build this yourself you should be able to call cmake configure with the provided presets, debug for executable and default for dll.