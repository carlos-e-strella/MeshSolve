// Using the clapack fortran to c interface plus some lines to stop conficts with c definitions
// I REFUSE to build lapacke from source
#include "lapack_interface.h"
#include <cblas.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "meshsolve.h"

typedef double (*Function2d)(double, double);

typedef struct Region {
    double* element_map;
    int* face_map;
    int element_map_size;
    int element_num;
    int shape_num;
    int point_num;
} Region;

typedef enum Term {
    SOURCE,
    TRANSIENT,
    DIFFUSION
} Term;

typedef struct Equation {
    Region region;
    Term* terms;
    int term_num;
    Function2d force_function;
} Equation;

// Bandwidth is required for lapack computations
// Also just helpful in general
typedef struct Matrix {
    double* value;
    double coefficient;
    int upper_band;
    int lower_band;
    size_t sizex;
    size_t sizey;
} Matrix;

// Reference for arrays that will be freed later
double* array_reference[128];
int array_index = 0;

// Printing for matrices in column major order (what cblas uses)
void print_matrix(Matrix matrix, char* name, int precision) {
    (name) ? printf("%s:\n", name) : printf("M:\n");
    size_t sizex = matrix.sizex; size_t sizey = matrix.sizey;

    for (int i = 0; i < sizey; i++) {
        printf(" | ");

        for (int j = 0; j < sizex; j++) {
            double number = matrix.value[i + j*sizey];

            // Sometimes happens???
            if (number == -0.0) { number = 0.0; }

            int truncation = (int) number;
            int digits = 0;
            if (number >= 0.0) { printf(" "); }

            if (truncation == 0) {
                digits = 1;
            }
            else {
                truncation = abs(truncation);
                while (truncation != 0) {
                    truncation /= 10;
                    digits++;
                }
            }

            // Max of 4 decimal places + 1 units place
            digits = precision + 1 - digits;
            if (digits < 0) { digits = 0; }
            printf("%.*f", digits, number);
            if (j < sizex - 1) { printf(" "); }
        }

        printf("  |\n");
    }
}

Matrix band_store(Matrix matrix) {
    size_t sizex = matrix.sizex; size_t sizey = matrix.sizey;
    int KL = matrix.lower_band; int KU = matrix.upper_band;

    int new_sizey = 2*KL + KU + 1;
    double* band_matrix_value = array_constructor(sizex * new_sizey);
    for (int i = 0; i < sizey; i++) {
        for (int j = 0; j < sizex; j++) {
            int index = i + j*sizey;
            double number = matrix.value[index];
            int upper_index = j - i; int lower_index = i - j;
            if (upper_index == 0 && lower_index == 0) {
                band_matrix_value[KU + j*new_sizey + KL] = number;
            }
            else if (upper_index > 0 && KU >= upper_index) {
                band_matrix_value[KU - upper_index + j*new_sizey + KL] = number;
            }
            else if (lower_index > 0 && KL >= lower_index) {
                band_matrix_value[KU + lower_index + j*new_sizey + KL] = number;
            }
        }
    }
    Matrix band_matrix;
    band_matrix.value = band_matrix_value;
    band_matrix.lower_band = KL; band_matrix.upper_band = KU;
    band_matrix.sizex = sizex; band_matrix.sizey = new_sizey;
    return band_matrix;
}

double* create_array(size_t size) {
    double* memory = (double*) malloc(size * sizeof(double));
    if (!memory) {
        printf("Error during memory allocation\n");
        return NULL;
    }
    for (int i = 0; i < size; i++) {
        memory[i] = 0;
    }
    return memory;
}

double* array_constructor(size_t size) {
    double* new_array = create_array(size);
    array_reference[0] = new_array;
    array_index++;
    return new_array;
}

void array_destructor() {
    for (int i = 0; i <= 128; i++) {
        free(array_reference[i]);
        array_reference[i] = NULL;
    }
}

// Apparently openblas does not have this?
double determinant_2d(double* matrix) {
    return fabs(matrix[0] * matrix[3] - matrix[1] * matrix[2]);
}

// Used for conversion from pyvista meshes to elements
Region get_region(double* points, int* face_map, size_t face_map_size, size_t points_size) {
    /*
    Face map is the input for the pyvista face definitions. A number N of elements is defined
    and then the N following members of the array are the index of the points at the points
    array for said element. Of course in python points is a matrix so you can access individual
    points, but here once you map those points for each element you have to double the number
    of members in the resulting array slice. Since realistically you would only imput tri or
    quad meshes there is no need to keep checking the number of points per element very time
    */

    int i = 1;
    int j = 0;
    int point_count = 0;
    Region region;
    
    double* elements = array_constructor(face_map_size * 2); 
    int shape_num = face_map[0];

    while (i < face_map_size) {
        int index = face_map[i] * 2;
        if (point_count < shape_num) {
            elements[j*2] = points[index];
            elements[j*2+1] = points[index+1];
            j++;
            point_count++;
        }  
        else {
            point_count = 0;
        }
        i++;
    }

    //Names here are actually really bad lol
    // j is the counter of total points described in the face map

    region.element_map = elements; // List of elements by their points
    region.element_map_size = j * 2; // Number of total x and y points in said list
    region.element_num = face_map_size / (shape_num + 1); // Number of actual elements
    region.face_map = face_map;
    region.shape_num = shape_num; // The number of points per element AKA number of shape functions
    region.point_num = points_size / 2;

    return region;
}

// Shape function interpolations for a quadrilateral reference element
double* shape_funcs_quad(double xi, double eta) {
    static double shape_funcs[4];
    shape_funcs[0] = 0.25 * (1-xi) * (1-eta);
    shape_funcs[1] = 0.25 * (1+xi) * (1-eta);
    shape_funcs[2] = 0.25 * (1-xi) * (1+eta);
    shape_funcs[3] = 0.25 * (1+xi) * (1+eta);
    return shape_funcs;
}

// Gradients of the shape functions for quadrilateral reference element
// dN_dxi are the even numbered elements and dN_deta the odd ones
double* grad_shape_funcs_quad(double xi, double eta) {
    static double grad_shape_funcs[8];
    grad_shape_funcs[0] = -0.25 * (1-eta); grad_shape_funcs[1] = -0.25 * (1-xi);
    grad_shape_funcs[2] = 0.25 * (1-eta); grad_shape_funcs[3] = -0.25 * (1+xi);
    grad_shape_funcs[4] = -0.25 * (1+eta); grad_shape_funcs[5] = 0.25 * (1-xi);
    grad_shape_funcs[6] = 0.25 * (1+eta); grad_shape_funcs[7] = 0.25 * (1+xi);
    return grad_shape_funcs;
}

// Jacobian of the transformation from quadrilateral reference element
double* jacobian_quad(double* element, double* gradients) {
    static double jacobian[4];
    double* dxi = &(gradients[0]);
    double* deta = &(gradients[1]);
    double* x = &(element[0]);
    double* y = &(element[1]);

    jacobian[0] = cblas_ddot(4, x, 2, dxi, 2);
    jacobian[1] = cblas_ddot(4, x, 2, deta, 2);
    jacobian[2] = cblas_ddot(4, y, 2, dxi, 2);
    jacobian[3] = cblas_ddot(4, y, 2, deta, 2);

    return jacobian;
}

// Integration with gauss-legendre quadrature
// TODO add more order options
double* integral_quad(double* (*integrand)(double), int shape_num) {
    int order = 2;
    double gauss_points[2] = {-1/sqrt(3), 1/sqrt(3)};
    // Weights for 2 points are 1, 1 soo...
    
    // This here will definitely exceed 128 arrays
    double* quadrature = create_array(shape_num);
    for (int i = 0; i < order; i++) {
        double* integral = integrand(gauss_points[i]);
        cblas_daxpy(shape_num, 1, integral, 1, quadrature, 1);
    }
        
    return quadrature;
}

// Double integration
double* dbl_integral_quad(double* (*integrand)(double, double), int shape_num) {
    int order = 2;
    double gauss_points[2] = {-1/sqrt(3), 1/sqrt(3)};
    // Weights for 2 points are 1, 1 soo...

    double* quadrature = create_array(shape_num);
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            double* integral = integrand(gauss_points[i], gauss_points[j]);
            cblas_daxpy(shape_num, 1, integral, 1, quadrature, 1);
            free(integral); integral = NULL;
        }
    }

    return quadrature;
}

// The actually useful part
// TODO add triangular elements
double quad_element_buffer[8];
Function2d function_2d_buffer;

// Integral for source term f on quadrilateral element
// This integral returns a vector, and is not dependent on the unknown parameters u
double* source_integrand(double x, double y) {
    // Size four arrays
    double* shape_funcs = shape_funcs_quad(x, y);
    double* gradients = grad_shape_funcs_quad(x, y); 

    // 2x2 matrix
    double* jacobian = jacobian_quad(quad_element_buffer, gradients);
    double determinant = determinant_2d(jacobian);

    // Size 4 element vector
    // Element VECTOR: Fi = Ni * f
    double* integrand = create_array(4);
    if (determinant) {
        double f = function_2d_buffer(x, y) * determinant;
        cblas_daxpy(4, f, shape_funcs, 1, integrand, 1);
    }
    else {
        printf("WARN: Singular jacobian. Check your mesh geometry\n");
    }
    return integrand;
}

// Integral for transient term ∂T/∂t on a quadrilateral element
double* transient_integrand(double x, double y) {
    // Size four arrays
    double* shape_funcs = shape_funcs_quad(x, y);
    double* gradients = grad_shape_funcs_quad(x, y); 
    double* jacobian = jacobian_quad(quad_element_buffer, gradients);
    double determinant = determinant_2d(jacobian);

    // 4x4 element matrix
    // Element matrix: Kij = Ni * Nj
    double* integrand = create_array(16);
    if (determinant) {
        cblas_dger(CblasColMajor, 4, 4, determinant, shape_funcs, 1, shape_funcs, 1, integrand, 4);
    }
    else {
        printf("WARN: Singular jacobian. Check your mesh geometry\n");
    }
    return integrand;
}

// Integral for diffusion term ∇²T on a quadrilateral element
double* diffusion_integrand(double x, double y) {
    // Size four arrays
    double* shape_funcs = shape_funcs_quad(x, y);
    double* gradients = grad_shape_funcs_quad(x, y);
    double* jacobian = jacobian_quad(quad_element_buffer, gradients);
    double determinant = determinant_2d(jacobian);

    // Element matrix: Kij ∇Ni · ∇Nj 
    double* integrand = create_array(16);
    integer ipiv[2]; // Pivot matrix required for decomposition
    integer size = 2;  
    integer info;

    double work[64]; integer work_size = 64;

    if (determinant) {
        dgetrf_(&size, &size, jacobian, &size, ipiv, &info);
        dgetri_(&size, jacobian, &size, ipiv, work, &work_size, &info);
        // jacobian pointer now holds the inverse jacobian

        // Generate gradients in stiffness coordinates for 4 shape functions
        // Multiplies inverse transpose jacobian to each shape function gradient
        double stiffness_gradients[8];
        for (int i = 0; i < 4; i++) {
            double* gradient = &(gradients[i*2]);
            double* stiffness_gradient = &(stiffness_gradients[i*2]);
            cblas_dgemv(CblasColMajor, CblasTrans, 2, 2, 1, jacobian, 2, gradient, 1, 0, stiffness_gradient, 1);
        }

        double* x = &(stiffness_gradients[0]);
        double* y = &(stiffness_gradients[1]);

        cblas_dger(CblasColMajor, 4, 4, determinant, x, 2, x, 2, integrand, 4);
        cblas_dger(CblasColMajor, 4, 4, 1, y, 2, y, 2, integrand, 4);
    }
    else {
        printf("WARN: This element's jacobian is singular. Bad mesh geometry much?\n");
    }
    return integrand;
}

// stiffness matrix creation
// TODO support multiple dof 
Matrix matrix_factory(Region region, Term term, int dof) {
    double* element_map = region.element_map;
    int* face_map = region.face_map;
    int element_map_size = region.element_map_size; // Not used here
    int element_num = region.element_num;
    int shape_num = region.shape_num;
    int point_num = region.point_num;

    // Each point contributes with its own degrees of freedom
    int total_dof = point_num * dof;
    Matrix matrix;
    double* stiffness_matrix = array_constructor(total_dof * total_dof);
    int lower_band = 0; int upper_band = 0;

    int is_vector = 0; // Source integral returns a vector instead of a matrix
                       // This would cause issues when assembeling without this flag

    // Assembly into stiffness matrix per element
    // Currently only supports quad elements
    for (int i = 0; i < element_num; i++) {
        double* element = &(element_map[i*shape_num*2]);
        memcpy(quad_element_buffer, element, 8 * sizeof(double));

        int shape_array_num = shape_num * shape_num;
        double* element_matrix;

        switch (term) {
            case SOURCE:
                element_matrix = dbl_integral_quad(source_integrand, shape_num);
                is_vector = 1;
                break;
            case TRANSIENT: 
                element_matrix = dbl_integral_quad(transient_integrand, shape_array_num);
                break;
            case DIFFUSION:
                element_matrix = dbl_integral_quad(diffusion_integrand, shape_array_num);
        }

        int* map = &(face_map[i*(shape_num+1)+1]);

        // Adds the value of the source integral to its corresponding point on the stiffness vector
        if (is_vector) {
            for (int j = 0; j < shape_num; j++) {
                int index = map[j];
                stiffness_matrix[index] += element_matrix[j];
            }
        }
        // Adds the value of integrals to their corresponding point on the stiffness matrix
        else {
            for (int j = 0; j < shape_num; j++) {
                for (int k = 0; k < shape_num; k++) {
                    double element = element_matrix[shape_num * j + k];
                    int index = total_dof * map[j] + map[k];
                    stiffness_matrix[index] += element;

                    if (element) {
                        int new_upper = map[j] - map[k]; int new_lower = map[k] - map[j];
                        upper_band = (new_upper > upper_band) ? new_upper : upper_band;
                        lower_band = (new_lower > lower_band) ? new_lower : lower_band;
                    }
                }
            }
        }
        free(element_matrix); element_matrix = NULL;
    }

    matrix.value = stiffness_matrix;
    matrix.lower_band = lower_band;
    matrix.upper_band = upper_band;
    matrix.sizex = (is_vector) ? 1 : total_dof;
    matrix.sizey = total_dof;
    return matrix;
}

// Assumes 1 dof per point
Matrix evaluate(Equation equation) {
    int term_num = equation.term_num;
    Region region = equation.region;

    // Currently only one force function can be specified
    // Geniuenly, why would you have more?
    function_2d_buffer = equation.force_function;
    Matrix force_vector;
    Matrix stiffness_matrix; 

    for (int i = 0; i < term_num; i++) {
        Term term = equation.terms[i];
        Matrix matrix = matrix_factory(region, term, 1);
        if (term == SOURCE) {
            // Given that matrix factory returned a vector sizex would be 1
            force_vector = matrix;
        }
        // TODO add more terms
        else {
            stiffness_matrix = matrix;
        }
    }

    if (!stiffness_matrix.value) {
        printf("FATAL ERROR: No stiffness matrix definition.");
    }
    if (!force_vector.value) {
        force_vector.value = create_array(stiffness_matrix.sizey);
        force_vector.sizex = 1;
        force_vector.sizey = stiffness_matrix.sizey;
    }

    // Actually solving the linear system
    Matrix band_matrix = band_store(stiffness_matrix);
    
    integer size = stiffness_matrix.sizex;
    integer ldab = band_matrix.sizey;
    integer lower_band = stiffness_matrix.lower_band;
    integer upper_band = stiffness_matrix.upper_band;
    integer* ipiv = (integer *) malloc(size * sizeof(integer)); 
    integer info;
    char no_transpose = 'N';
    integer nrhs = 1;

    dgbtrf_(&size, &size, &lower_band, &upper_band, band_matrix.value, &ldab, ipiv, &info);
    if (info > 0) {
        printf("FATAL ERROR: Stiffness matrix is singular! Probably due to unconstrained equation.");
        exit(1);
    }
    dgbtrs_(&no_transpose, &size, &lower_band, &upper_band, &nrhs, band_matrix.value, &ldab, ipiv, 
            force_vector.value, &size, &info);
    free(ipiv);
    free(stiffness_matrix.value);

    return force_vector;
}

double test_force(double x, double y) {
    return x * y;
}

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

// Test function with values collected from a pyvista mesh
// The mesh is a square subdivided into 4 faces and 9 points

int main() {
    int ispec = 1;
    int face_map[20] = {4, 0, 1, 4, 3, 4, 1, 2, 5, 4, 4, 3, 4, 7, 6, 4, 4, 5, 8, 7};
    double points[18] = {-0.5, -0.5, 0.0, -0.5, 0.5, -0.5, -0.5, 0.0, 0.0, 
        0.0, 0.5, 0.0, -0.5, 0.5, 0.0, 0.5, 0.5, 0.5};

    Region region = get_region(points, face_map, 20, 18);
    test_function(region);
    return 0;
}