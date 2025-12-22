"""
An actually worked through example of a MeshSolve implementation.
Here we are solving the heat equation for a 40 x 40 square mesh.
Please forgive my very messy python code.
"""

from ctypes import *
import pyvista as pv
import numpy as np

# You have to actually build the dll first
library_path = './build/meshsolve.dll'

try:
    lib = windll.LoadLibrary(library_path)
except OSError:
    lib = cdll.LoadLibrary(library_path)
except FileNotFoundError:
    print("Cannot find meshsolve.dll. Make sure you have actually built it.")

# This function provides a square mesh with n x n elements
def square_mesh(n):
    return pv.Plane(i_resolution=n, j_resolution=n)


"""
ctypes expects the functions in the dll to be declared here in python with the respective
type of their arguments. Except python does not have the same types or structures as the
dll, so we have to declare them here.
"""

FUNCTION2D = CFUNCTYPE(c_double, c_double, c_double)

class REGION(Structure):
    _fields_ = [
        ("element_map", POINTER(c_double)),
        ("face_map", POINTER(c_int)),
        ("element_map_size", c_int),
        ("element_num", c_int),
        ("shape_num", c_int),
        ("point_num", c_int)
    ]

    def __repr__(self):
        return f"Region with {self.point_num} points and {self.element_num} elements."

class MATRIX(Structure):
    _fields_ = [
        ("value", POINTER(c_double)),
        ("coefficient", c_double),
        ("upper_band", c_int),
        ("lower_band", c_int),
        ("sizex", c_size_t),
        ("sizey", c_size_t)
    ]

    def __repr__(self):
        return f"Matrix of size {self.sizex} x {self.sizey}"

class EQUATION(Structure):
    _fields_ = [
        ("region", REGION),
        ("terms", POINTER(c_int)),
        ("term_num", c_int),
        ("force_function", FUNCTION2D),
        ("f", MATRIX),
        ("k", MATRIX),
        ("m", MATRIX),
        ("ic", MATRIX)
    ]


# Now we declare the functions with their arguments and return types
# THESE WILL APPEAR LATER IN THE CODE
get_region = lib[5]
get_region.argtypes = [
    POINTER(c_double),
    POINTER(c_int),
    c_size_t,
    c_size_t
]
get_region.restype = REGION

assembler = lib[2]
assembler.argtypes = [POINTER(EQUATION)]
assembler.restype = None

time_step = lib[7]
time_step.argtypes = [
    POINTER(EQUATION),
    c_double,
]
time_step.restype = None

equation_init = lib[3]
equation_init.argtypes = None
equation_init.restype = EQUATION


"""
This is the part you can actually customize! Several things can and *will* brake though.
A couple of pointers: 
  - A mesh bigger than 50 will begin causing performace issues. 100 onward is undefined 
    behaviour territory since C sometimes is not able to find enough contiguous space in 
    memory. 
  - When increasing mesh size you also have to decrease the step size since the heat
    equation is time dependant, and the solver will scream if dt is too high.
  - When running non time dependant equations the solution is stored in equation.f and not
    in the initial conditions vector.
"""

"""
You can define sources for your equations in this way. First you declare the function
instance and then you can set the equation.force_function member to it. I do not actually
use this in the demo but i thought it would be good to know.
"""
def test_callback(a, b):
    return a * b
function_instance = FUNCTION2D(test_callback)

test_mesh = square_mesh(40)

# This is just doing type conversion of the mesh members to a format the C code can digest
points = test_mesh.points
x_points = points[:, 0]
y_points = points[:, 1]
c_points = np.stack((x_points, y_points), axis=1).ravel().astype(np.float64)
len_points = len(c_points)
points_instance = (c_double * len_points)(*c_points)

faces = test_mesh.faces # More info on this particular member in the MeshSolve source code.
len_faces = len(faces)
faces_instance = (c_int * len_faces)(*faces)


"""
NOW THE ACTUAL EQUATION
Here we define the heat equation:
    dT/dt = ∇²T <- Transient term = Diffusion term
    (The checking which side of the equation they are on part is quite messed up though)
Those are the terms that will form our terms array.

Next is the region: We pass the processed members points_instance and faces_instance and
their lenghts to the get_region function to get our Region struct.

After that we define the initial conditions (ic). I just picked arbitrary values for this
demo. Keep in mind that I actually pass those values to the Matrix struct constructor that
was defined earlier, since MeshSolve uses it to simplify passing information about
matrices and vectors.

Finally, all of these members get converted to C types and they are passed to the equation
struct. Note that I first called the equation_init function to make sure no funny business
is going on with undefined terms picking random values.
"""

terms = [1, 2] # 0: SOURCE  1: TRANSIENT  2: DIFFUSION
len_terms = len(terms)
terms_instance = (c_int * len_terms)(*terms)

new_region = get_region(points_instance, faces_instance, len_faces, len_points)

ic = (np.sin(x_points * 20) + np.cos(y_points * 26.3)) * 0.5
ic_len = len(ic)
ic_instance = (c_double * ic_len)(*ic)
ic_matrix = MATRIX(ic_instance, 0, 0, 0, 1, ic_len)

equation = equation_init()
equation.region = new_region
equation.terms = terms_instance
equation.term_num = len_terms
equation.force_function = function_instance
equation.ic = ic_matrix


"""
Now we assemble the linear form of our differential equation! By calling the assembler
function we get either the full solution or the ODE ready for time stepping. Note that for
any other manipulation we do to our equation we will be passing it by reference by using 
the byref function ctypes provides. I am using the Runge Kutta 4 method for the ODE part.
"""

assembler(byref(equation)) 

# Plotting
test_mesh['temp'] = ic
warped_mesh = test_mesh.warp_by_scalar(scalars='temp', factor=0.5)
warped_mesh['temp'] = ic

plotter = pv.Plotter()
plotter.add_mesh(warped_mesh, show_edges=True, scalars='temp', cmap='curl')
plotter.show(interactive_update=True, auto_close=False)

# Once again, going too high on the step size will make the solver SCREAM
step_size = c_double(0.00000001)
steps = 150


"""
Finally the full display of the solution! Since we are doing time stepping we will collect
the values of the system at any given time from the ic vector returned by the equation.
Once again, you would get your complete solution from the force vector if you were doing
non time dependant equations. Im 
"""
for j in range(steps):
    time_step(byref(equation), step_size)
    matrix = equation.ic

    # Here we convert from ctypes back to numpy arrays
    # THIS IS WHERE THE PERFORMANCE ISSUES HAPPEN BY THE WAY
    new_list_type = POINTER(c_double * matrix.sizey)
    new_list_buffer = cast(matrix.value, new_list_type).contents
    new_list = np.frombuffer(new_list_buffer, dtype=np.float64)

    test_mesh['temp'] = new_list
    warped_mesh.points = test_mesh.warp_by_scalar(scalars='temp', factor=0.5).points
    warped_mesh['temp'] = new_list
    plotter.update() # Finally updating every step

# This is the end of the demo!