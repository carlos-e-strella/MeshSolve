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

typedef struct Matrix {
    double* value;
    double coefficient;
    int upper_band;
    int lower_band;
    size_t sizex;
    size_t sizey;
} Matrix;

typedef struct Equation {
    Region region;
    Term* terms;
    int term_num;
    Function2d force_function;
    Matrix f;
    Matrix k;
    Matrix m;
    Matrix ic;
} Equation;

API void print_matrix(Matrix matrix, char* name, int precision);

API Region get_region(double* points, int* face_map, size_t face_map_size, size_t points_size);

API void evaluate(Equation* equation);

API void step_rtk(Equation* equation, double step_size);

API void equation_assembler(Equation* equation);

API void clear_equation(Equation* equation);

API Equation equation_init();

#endif