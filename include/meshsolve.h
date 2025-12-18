#ifndef MESHSOLVE_H
#define MESHSOLVE_H

#ifdef _WIN32
    #define API __declspec(dllexport)
#else
    #define API __attribute__((visibility("default")))
#endif

#include <stdlib.h>

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

API void print_matrix(double* matrix, int sizex, int sizey);

API double* create_array(size_t size);

API double* array_constructor(size_t size);

API void array_destructor();

API double determinant_2d(double* matrix);

API Region get_region(double* points, int* face_map, size_t face_map_size, size_t points_size);

API double* shape_funcs_quad(double xi, double eta);

API double* grad_shape_funcs_quad(double xi, double eta);

API double* jacobian_quad(double* element, double* gradient);

API double* quad_integral(double* (*integrand)(double), int shape_num);

API double* dbl_quad_integral(double* (*integrand)(double, double), int shape_num);

API double* transient_integrand(double x, double y);

API double* transient_wrapper(double* element, int shape_num);

API double* matrix_factory(Region region, Term term, int dof);

API void test_function(double* points, int* face_map, size_t face_map_size, size_t points_size);

#endif