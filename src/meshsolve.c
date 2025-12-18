// Using the clapack fortran to c interface plus some lines to stop conficts with c definitions
// I REFUSE to build lapacke from source
#include "lapack_interface.h"
#include <cblas.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "meshsolve.h"

typedef struct Region {
    double* element_map;
    int* face_map;
    int element_map_size;
    int element_num;
    int shape_num;
    int point_num;
} Region;

typedef enum Term {
    TRANSIENT,
    DIFFUSION
} Term;

// Reference for arrays that will be freed later
double* array_reference[128];
int array_index = 0;

// Printing for matrices in column major order (what cblas uses)
void print_matrix(double* matrix, int sizex, int sizey) {
    for (int i = 0; i < sizey; i++) {
        for (int j = 0; j < sizex; j++) {
            printf("%.4f ", matrix[i + j*sizey]);
        }
        printf("\n");
    }
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

// Integral for transient term ∂T/∂t on a quadrilateral element
double* transient_integrand(double x, double y) {
    // Size four arrays
    double* shape_funcs = shape_funcs_quad(x, y);
    double* gradients = grad_shape_funcs_quad(x, y); 

    // 2x2 matrix
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

        // Generate gradients in global coordinates for 4 shape functions
        // Multiplies inverse transpose jacobian to each shape function gradient
        double global_gradients[8];
        for (int i = 0; i < 4; i++) {
            double* gradient = &(gradients[i*2]);
            double* global_gradient = &(global_gradients[i*2]);
            cblas_dgemv(CblasColMajor, CblasTrans, 2, 2, 1, jacobian, 2, gradient, 1, 0, global_gradient, 1);
        }

        double* x = &(global_gradients[0]);
        double* y = &(global_gradients[1]);

        cblas_dger(CblasColMajor, 4, 4, determinant, x, 2, x, 2, integrand, 4);
        cblas_dger(CblasColMajor, 4, 4, 1, y, 2, y, 2, integrand, 4);
    }
    else {
        printf("WARN: This element's jacobian is singular. Bad mesh geometry much?\n");
    }
    return integrand;
}

// Global matrix creation
// TODO support multiple dof and multiple terms
double* matrix_factory(Region region, Term term, int dof) {
    double* element_map = region.element_map;
    int* face_map = region.face_map;
    int element_map_size = region.element_map_size;
    int element_num = region.element_num;
    int shape_num = region.shape_num;
    int point_num = region.point_num;

    // Each point contributes with its own degrees of freedom
    int total_dof = point_num * dof;
    double* global_matrix = array_constructor(total_dof * total_dof);

    // Assembly into global matrix per element
    // Currently only supports quad elements
    for (int i = 0; i < element_num; i++) {
        double* element = &(element_map[i*shape_num*2]);
        memcpy(quad_element_buffer, element, 8 * sizeof(double));

        int shape_array_num = shape_num * shape_num;
        double* element_matrix;
        switch (term) {
            case TRANSIENT: 
                element_matrix = dbl_integral_quad(transient_integrand, shape_array_num);
                break;
            case DIFFUSION:
                element_matrix = dbl_integral_quad(diffusion_integrand, shape_array_num);
        }
        element_matrix = transient_wrapper(element, shape_num);

        int* map = &(face_map[i*(shape_num+1)+1]);

        for (int j = 0; j < shape_num; j++) {
            for (int k = 0; k < shape_num; k++) {
                int index = total_dof * map[j] + map[k];
                global_matrix[index] += element_matrix[shape_num * j + k];
            }
        }
        free(element_matrix); element_matrix = NULL;
    }

    return global_matrix;
}

void test_function(double* points, int* face_map, size_t face_map_size, size_t points_size) {
    Region region = get_region(points, face_map, face_map_size, points_size);
    double* global_matrix = matrix_factory(region, TRANSIENT, 1);
    print_matrix(global_matrix, region.point_num, region.point_num);
    array_destructor;
}

// Test function with values collected from a pyvista mesh
// The mesh is a square subdivided into 9 parts

int main() {
    int ispec = 1;
    int face_map[20] = {4, 0, 1, 4, 3, 4, 1, 2, 5, 4, 4, 3, 4, 7, 6, 4, 4, 5, 8, 7};
    double points[18] = {-0.5, -0.5, 0.0, -0.5, 0.5, -0.5, -0.5, 0.0, 0.0, 
        0.0, 0.5, 0.0, -0.5, 0.5, 0.0, 0.5, 0.5, 0.5};

    test_function(points, face_map, 20, 18);
    return 0;
}