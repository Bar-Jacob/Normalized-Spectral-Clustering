#ifndef spkmeans
#define spkmeans
#define PY_SSIZE_T_CLEAN
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <Python.h>
#include "kmeans.c"

typedef struct Eigenvector
{
    double* vector;
    double value;
} Eigenvector;

/*goals functions*/
void wam_goal(double** points, int num_of_points, int dimension);
void ddg_goal(double** points, int row, int dimension);
void lnorm_goal(double** points, int row, int dimension);
void jacobi_goal(double** sym_matrix, int row);
void spk_goal(double** points, int row, int col, int k);
void spk_goal_python(double** points, int row, int col, int k);

/*jacobi related functions*/
Eigenvector* jacobi(double** points, int row);
double* s_c_t_calculation(double** lp_matrix, int row);
double** rotation_matrix(double** lp_matrix, int row, double* sct_val);
double* largest_off_digonal_value(double** lp_matrix, int row);
double** jacobi_calculations(double** lp_matrix, int row,
                                        double* sct_vals);
int is_convergence(double** lp_matrix, double** lp_matrix_tag, int row);
Eigenvector* creating_eigenvector_array(double** eigenvectors_matrix, 
                                double** eigenvalues_matrix, int row);
double sos_off(double** lp_matrix, int row);

/*matrix algebra functions*/
double** zero_matrix(int row, int col);
double** identity_matrix(int row, int col);
double** matrix_multiplication(double** matrix1,
    double** matrix2, int row, int col);
double** three_matrix_multiplication(double** matrix1,
    double** matrix2, double** matrix3, int row, int col);
void matrix_substraction(double** matrix1, double** matrix2, int row, int col);

/*spk related functions*/
int eigengap_heuristic(Eigenvector* eigenvector, int row);
double** creating_U(Eigenvector* eignvector, int k, int row);
void renormalizing_U(double** U, int k, int row);
double** adjacency_matrix(double** points, int row, int col);
double** diagonal_degree_matrix(double** points, int row, int dimension);
double** normalized_graph_laplacian(double** points, int row, int dimension);
void manipulated_diagonal(double** diag_matrix, int row);

/*sorting functions*/
void merge(Eigenvector* eigenvector, Eigenvector* R, Eigenvector* L, int l, int m, int r);
void merge_sort(Eigenvector* eigenvector, int l, int r);

/*print functions*/
void print_array(double* array, int row);
void print_matrix(double** matrix, int row, int col);

/*free memory function*/
void free_memory(double** matrix, int row);

/*module related function*/
double** convert_python_to_c(PyObject* data_points_p, int dimension_p, int num_of_points_p);
PyObject* cToPyObject(double** T, int dimension, int num_of_points);

#endif