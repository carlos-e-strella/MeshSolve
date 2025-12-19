#ifndef MESHSOLVE_H
#define MESHSOLVE_H

#ifdef _WIN32
    #define API __declspec(dllexport)
#else
    #define API __attribute__((visibility("default")))
#endif

#include <stdlib.h>

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

typedef struct Matrix {
    double* value;
    double coefficient;
    int upper_band;
    int lower_band;
    size_t sizex;
    size_t sizey;
} Matrix;

API void print_matrix(Matrix matrix, char* name, int precision);

API Matrix band_store(Matrix matrix);

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

API Matrix matrix_factory(Region region, Term term, int dof);

API Matrix evaluate(Equation equation);

API void test_function(Region region);

#endif