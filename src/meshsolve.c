#include <cblas.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Reference for arrays that will be freed later
double* array_reference[128];
int array_index = 0;

// Printing for matrices in column major order (what cblas uses)
void print_matrix(double* matrix, int sizex, int sizey) {
    for (int i = 0; i < sizey; i++) {
        for (int j = 0; j < sizex; j++) {
            printf("%.2f ", matrix[i + j*sizey]);
        }
        printf("\n");
    }
}

double* create_array(size_t size) {
    double* memory = (double*) malloc(size * sizeof(double));
    if (!memory) {
        printf("Error during memory allocation");
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
double* get_2d_elements(int* points, int* face_map) {
    int i = 0;
    size_t face_num = sizeof(face_map) / sizeof(int);
    
    double* elements = array_constructor(face_num * 2);

    while (i < face_num) {
        int index = face_map[i] * 2;
        elements[i*2] = points[index];
        elements[i*2+1] = points[index+1];
        i++;
    }

    return elements;
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
double* jacobian_quad(double* element, double* gradient) {
    static double jacobian[4];
    double* dxi; dxi = &(gradient[0]);
    double* deta; deta = &(gradient[1]);
    double* x; x = &(element[0]);
    double* y; y = &(element[1]);

    jacobian[0] = cblas_ddot(4, x, 2, dxi, 2);
    jacobian[1] = cblas_ddot(4, x, 2, deta, 2);
    jacobian[2] = cblas_ddot(4, y, 2, dxi, 2);
    jacobian[3] = cblas_ddot(4, y, 2, deta, 2);

    return jacobian;
}

// Integration with gauss-legendre quadrature
// TODO add more order options
double* quad_integral(double* (*integrand)(double), int shape_num) {
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
double* dbl_quad_integral(double* (*integrand)(double, double), int shape_num) {
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
double quad_element_buffer[4];

// Integral for transient term ∂T/∂t on a quadrilateral element
double* transient_integrand(double x, double y) {
    // Size four arrays
    double* shape_funcs = shape_funcs_quad(x, y);
    double* gradient = grad_shape_funcs_quad(x, y); 

    // 2x2 matrix
    double* jacobian = jacobian_quad(quad_element_buffer, gradient);
    double determinant = determinant_2d(jacobian);

    // 4x4 element matrix
    // N ⊗ N for each quadrature point
    double* integrand = create_array(16);
    cblas_dger(CblasColMajor, 4, 4, determinant, shape_funcs, 1, shape_funcs, 1, integrand, 4);
    return integrand;
}

// Wrapper for transient elements, only for quadrilateral elements
double* transient_wrapper(double* element, int shape_num) {
    memcpy(quad_element_buffer, element, 4 * sizeof(double));
    double* element_stiffness = dbl_quad_integral(transient_integrand, shape_num * shape_num);
    return element_stiffness;
}

// TODO ensure safety when passing the number of shape functions
int main() {
    double test_element[8] = {8.0, 2.0, 5.0, 6.0, 1.0, 2.0, 8.0, 9.0};
    double* matrix = transient_wrapper(test_element, 4);
    print_matrix(matrix, 4, 4);
    free(matrix); matrix = NULL;

    return 0;
}