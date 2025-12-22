#  MeshSolve
A c rewrite of my python pde solver. In the future I will implement this as a python module. 
I am not a mathematician nor an engineer, so this program is probably packed with errors.

This project uses the lapack and openblas linear algebra libraries for c, as well as cmake and vcpkg for build and package management respectively.

If you are running the python demo you will have to install pyvista and numpy in your environment.
I believe ctypes comes with python though.

## How do I use this program?
Solutions are computed through calls to the `equation_assembler` function. Each call requires an `Equation` struct to be provided, which has to include relevant information on the mesh, the equation terms and initial conditions if you are doing a time dependant equation.
Currently there are only three terms that can be declared as type `Term`:
- `SOURCE`
- `TRANSIENT` 
- `DIFFUSION`

You can build the meshsolve dll and use it in other c programs (I geniuenly have no idea how you would do that) *or* you can go through the python wrapper implementation in `demo.py`, although there will be some performance issues if you try it with bigger meshes (I'm working on that). Here is the example provided in the demo file, where MeshSolve is used to solve the heat equation:

![Heat equation demo gif](images/heat_demo.gif)

More information is (or will be added) to the code itself.

### Building
The cmake files *should* include openblas and lapack if you have vcpkg installed and I also provided the clapack fortran to c interface because I refuse to build lapacke from source. If you want to build this yourself you should be able to call cmake configure with the provided presets, debug for executable and default for dll.

## The Finite Element Method

### How does this program solve PDEs?
I have implemented (barely) the finite element method, which is a numerical method to aproximate solutions for PDEs.
The method requires you to specify a mesh where your equations occurr and calculates a stiffness matrix for it.
The solution then involves solving a linear system involving said matrix and your initial conditions (and your boundary conditions).

### How does this program implement the finite element method?
This program provides functions to calculate the stiffness matrix for a specific mesh.
Once its finished, you will be able to declare your partial differential equation as a sum of different terms.
Each term has a weak form that is used to calculate its contribution to the stiffness matrix.
Finally, the program will calculate the resulting solution by solving a system of equations involving said matrix.