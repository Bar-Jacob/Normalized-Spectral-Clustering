#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

double** zero_matrix(int row, int col);
void print_matrix(double** matrix, int row, int col);
double** identity_matrix(int row, int col);
double** adjacency_matrix(double** points, int row, int col);
double** diagonal_degree_matrix(double** points, int row, int col);
double** normalized_graph_laplacian(double** points, int row, int col);
double*** jacobi(double** points, int row, int col);
double* largest_off_digonal_value(double** lp_matrix, int row, int col);
double** rotation_matrix(double** lp_matrix, int row, int col);
double** jacobi_calculations(double** lp_matrix, int row,
                             int col, double* sct_vals);
double sos_off(double** lp_matrix, int row, int col);
int is_convergence(double** lp_matrix, double** lp_matrix_tag, 
                                            int row, int col, int first_iter);
double** matrix_multiplication(double** matrix1,
                               double** matrix2, int row, int col);
double** three_matrix_multiplication(double** matrix1,
                                     double** matrix2, double** matrix3, int row, int col);

void manipulated_diagonal(double** diag_matrix, int row);
void matrix_substraction(double** matrix1, double **matrix2, int row, int col);
double** transpose_matrix(double** matrix, int row, int col);
void free_memory(double** matrix, int row);
void wam_goal(double** points, int num_of_points, int dimension);
void ddg_goal(double** points, int row, int col);
void lnorm_goal(double** points, int row, int col);
void jacobi_goal(double** sym_matrix, int row, int col);
double* s_c_t_calculation(double** lp_matrix, int row, int col);