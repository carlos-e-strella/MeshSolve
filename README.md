#  MeshSolve
A c rewrite of my python pde solver. In the future I will implement this as a python module. 
I am not a mathematician nor an engineer, so this program is probably packed with errors.

This project uses the lapack and openblas linear algebra libraries for c, as well as cmake and vcpkg for build and package management respectively.

## How do I use this program?
Solutions are computed through calls to the `evaluate` function. Each call requires an `Equation` struct to be provided, which has to include relevant information on the mesh and the equation terms.
Currently there are only three terms that can be declared as type `Term` (with only two of them actually working):
- `SOURCE`
- `TRANSIENT` (there is no time stepping yet so this one does not work)
- `DIFFUSION`

I will add more info on these terms and the finite element method in general but for now I have provided `test_function` really just as a showcase of some *actual* solutions that can be computed by meshsolve.
```c
void test_function(Region region) {
    Equation equation;
    Term terms[2] = {SOURCE, DIFFUSION};

    // Initializing the equation struct members
    equation.region = region;
    equation.force_function = test_force;
    equation.terms = terms;
    equation.term_num = 2;

    Matrix solution = evaluate(equation);
    print_matrix(solution, "Solution", 5);
    free(solution.value);
    array_destructor;
}
```
More details are in the code as comments.

### Building
The cmake files *should* include openblas and lapack if you have vcpkg installed and I also provided the clapack fortran to c interface because I refuse to build lapacke from source. If you want yo build this yourself you should be able to call cmake configure with the provided presets, debug for executable and default for dll.

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