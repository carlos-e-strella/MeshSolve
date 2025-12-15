#include <cblas.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Reference for arrays that will be freed later
double* array_reference[128];
int array_index = 0;


// TODO What happens when you have more than 128 arrays???
double* array_constructor(size_t size) {
    double* new_array = (double*) malloc(size * sizeof(double));
    for (int i = 0; i < size; i++) {
        new_array[i] = 0;
    }
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

// Integration with gauss-legendre quadrature
// TODO add more order options
double* quad_integral(double* (*integrand)(double), int shape_num) {
    int order = 2;
    double gauss_points[2] = {-1/sqrt(3), 1/sqrt(3)};
    // Weights for 2 points are 1, 1 soo...
    
    // This here will definitely exceed 128 arrays
    double* quadrature = array_constructor(shape_num);
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

    double* quadrature = array_constructor(shape_num);
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            double* integral = integrand(gauss_points[i], gauss_points[j]);
            cblas_daxpy(shape_num, 1, integral, 1, quadrature, 1);
        }
    }

    return quadrature;
}

// TODO ensure safety when passing the number of shape functions
