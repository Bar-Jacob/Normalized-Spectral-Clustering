#include "spkmeans.h"

int main(int argc, char const *argv[])
{
    double **res;
    double **res2;
    double ***result;

    double* sorted;

    res = zero_matrix(3, 3);
    res2 = zero_matrix(3, 3);
    res[0][0] = 3;
    res[0][1] = 2;
    res[0][2] = 4;
    res[1][0] = 2;
    res[1][1] = 0;
    res[1][2] = 2;
    res[2][0] = 4;
    res[2][1] = 2;
    res[2][2] = 3;

    res2[0][0] = 1;
    res2[0][1] = 2;
    res2[0][2] = 3;
    res2[1][0] = 4;
    res2[1][1] = 5;
    res2[1][2] = 6;
    res2[2][0] = 7;
    res2[2][1] = 8;
    res2[2][2] = 9;

    sorted = (double*) calloc(6, sizeof(double));
    for (int i = 0; i < 6; i++)
    {
        sorted[i] = 6-i;
    }

    for (int i = 0; i < 6; i++)
    {
        printf("%lf", sorted[i]);
    }


    // jacobi_goal(res,3,3);
    // print_matrix(result[0], 3, 3);
    // printf("\n");
    // print_matrix(result[1], 3, 3);
    // print_matrix(result[0], 3,3);
    // printf("\n");
    // print_matrix(result[1],3,3);
    return 0;
}

void wam_goal(double** points, int num_of_points, int dimension)
{
    double** result;
    result = adjacency_matrix(points, num_of_points, dimension);
    print_matrix(result, num_of_points, dimension);
    free_memory(result, num_of_points);
}

void ddg_goal(double** points, int row, int col)
{
    double** result;
    result = diagonal_degree_matrix(points, row, col);
    print_matrix(result, row, col);
    free_memory(result, row);
}

void lnorm_goal(double** points, int row, int col)
{
    double** result;
    result = normalized_graph_laplacian(points, row, col);
    print_matrix(result, row, col);
    free_memory(result, row);
}

void jacobi_goal(double** sym_matrix, int row, int col)
{
    double*** result;
    double** trans_eigenvectors;
    int i;

    result = jacobi(sym_matrix, row, col);
    for (i = 0; i < row; i++)
    {
        if(i == row-1){
            printf("%.4f", result[1][i][i]);
        }else{
            printf("%.4f,", result[1][i][i]);
            }
    }
    printf("\n");

    trans_eigenvectors = transpose_matrix(result[0], row, col);
    print_matrix(trans_eigenvectors, row, col);
    free_memory(result[1], row);
    free_memory(result[0], row);
    free_memory(trans_eigenvectors, row);
    free(result);
}

// void spk_goal()

/*
computing adjacency matrix
*/
double** adjacency_matrix(double** points, int num_of_points, int dimension)
{
    double** result;
    int i = 0;
    int j = 0;
    int d = 0;
    double weight = 0.0;
    result = zero_matrix(num_of_points, num_of_points);
    /*
    for each point we will go through all of the other points
    only the upper triangle, since the matrix is symmetrical
    */
    for (i = 0; i < num_of_points; i++)
    {
        for (j = i + 1; j < num_of_points; j++)
        {
            for (d = 0; d < dimension; d++)
            {
                weight += pow((points[i][d] - points[j][d]), 2);
            }
            weight = sqrt(weight);
            weight = exp(-weight / 2);
            /*
            filling the upper and lower triangle in parallel
            */
            result[i][j] = weight;
            result[j][i] = weight;
        }
    }
    return result;
}

double** diagonal_degree_matrix(double** points, int row, int col)
{
    double** result;
    double** adj_matrix;
    int i = 0;
    int j = 0;
    result = zero_matrix(row, col);
    adj_matrix = adjacency_matrix(points, row, col);
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            result[i][i] += adj_matrix[i][j];
        }
    }
    free_memory(adj_matrix, row);
    return result;
}

/*
computing laplacian by substract the multiplication of 3 matrixes
(D^-0.5, W, D^-0.5, whereas D is the diagonal degree and W is the adjacency), 
from the identity matrix
*/
double** normalized_graph_laplacian(double** points, int row, int col)
{
    double** adj_matrix;
    double** diag_matrix;
    double** multip_matrix;
    double** result;

    adj_matrix = adjacency_matrix(points, row, col);
    diag_matrix = diagonal_degree_matrix(points, row, col);
    /*
    changes the matrix itself, since we give a pointer to the matrix 
    */
    manipulated_diagonal(diag_matrix, row);
    multip_matrix = three_matrix_multiplication(diag_matrix, adj_matrix,
                                                diag_matrix, row, col);
    free_memory(adj_matrix, row);
    free_memory(diag_matrix, row);
    result = identity_matrix(row, col);
    /*
    changes the first matrix given, so the result will be on the result matrix
    */
    matrix_substraction(result, multip_matrix, row, col);
    free_memory(multip_matrix, row);

    return result;
}

/*
The jacobi eigenvalue algorithm
*/
double*** jacobi(double** lp_matrix, int row, int col)
{
    double** lp_matrix_tag;
    double** eigenvectors_matrix;
    double** tmp_matrix;
    double** tmp_matrix2;
    double** rotate_matrix;
    double* sct_vals;
    double*** result;
    int cnt = 0;
    int first_iter = 1;
    int is_diag = 0;

    /*in order to avoid free function on a pointer that has'nt been
    defined yet, we define tmp_matrix2*/
    tmp_matrix2 = identity_matrix(row, col);
    eigenvectors_matrix = identity_matrix(row, col);

    /*
    The algorithm will stop after 100 iterations or if there is a convergence,
    whatever happens first
    */
    while (cnt != 100 && is_convergence(tmp_matrix2, lp_matrix_tag, 
                        row, col, first_iter) == 0 && is_diag == 0)
    {
        rotate_matrix = rotation_matrix(lp_matrix, row, col);
        sct_vals = s_c_t_calculation(lp_matrix, row, col);
        lp_matrix_tag = jacobi_calculations(lp_matrix, row, col, sct_vals);

        /*points to the last updated eigenvectors matrix*/
        tmp_matrix = eigenvectors_matrix;

       /*eigenvectors_matrix will now point to a new updated eigenvectors matrix*/
        eigenvectors_matrix = matrix_multiplication(tmp_matrix, rotate_matrix, row, col);

        /*release tmp_matrix and tmp_matrix2*/
        free_memory(tmp_matrix, row);
        free_memory(tmp_matrix2, row);
        
        if (is_diagonal(lp_matrix_tag, row, col) == 1) {
            is_diag = 1;
        }
        else {
            tmp_matrix2 = lp_matrix;
            lp_matrix = lp_matrix_tag;
        }
        cnt++;
        first_iter = 0;
    }

    free_memory(lp_matrix, row);
    free_memory(rotate_matrix, row);
    free(sct_vals);

    /*returns an array with two matrixes - eigenvectors and eigenvalues */
    result = (double***)calloc(2,sizeof(double**));
    result[0] = eigenvectors_matrix;
    result[1] = lp_matrix_tag;
    return result;
}

/*
builds the rotation matrix
*/
double** rotation_matrix(double** lp_matrix, int row, int col)
{
    double** result = identity_matrix(row, col);
    double* sct_val = s_c_t_calculation(lp_matrix, row, col);
    int i = (int)sct_val[3];
    int j = (int)sct_val[4];

    result[i][i] = sct_val[1];
    result[j][j] = sct_val[1];
    result[i][j] = sct_val[0];
    result[j][i] = -sct_val[0];

    free(sct_val);

    return result;
}
/*
calculates s, c, t
*/
double* s_c_t_calculation(double** lp_matrix, int row, int col)
{
    double* pivot_arr = largest_off_digonal_value(lp_matrix, row, col);
    double* result_arr = (double *)calloc(5, sizeof(double));
    int i = (int)pivot_arr[1];
    int j = (int)pivot_arr[2];
    double pivot_max = pivot_arr[0];
    double teta = (lp_matrix[j][j] - lp_matrix[i][i]) / (2 * pivot_max);

    double t = 0.0;
    if (teta < 0)
    {
        t = -1 / (fabs(teta) + sqrt(pow(teta, 2) + 1));
    }
    else
    {
        t = 1 / (fabs(teta) + sqrt(pow(teta, 2) + 1));
    }

    double c = 1 / (sqrt(pow(t, 2) + 1));
    double s = c * t;

    result_arr[0] = s;
    result_arr[1] = c;
    result_arr[2] = t;
    result_arr[3] = i;
    result_arr[4] = j;

    free(pivot_arr);
    return result_arr;
}

/*
calculating A_tag (only the cells that changed)
*/
double** jacobi_calculations(double** lp_matrix, int row,
                             int col, double* sct_vals)
{
    double** A_tag = zero_matrix(row, col);
    int n = 0;
    int m = 0;
    int i = sct_vals[3];
    int j = sct_vals[4];
    double c = sct_vals[1];
    double s = sct_vals[0];
    /*
    copy the laplacian matrix, in order to use the input for the calculations ahead
    */
    for (n = 0; n < row; n++)
    {
        for (m = 0; m < col; m++)
        {
            A_tag[n][m] = lp_matrix[n][m];
        }
    }
    /*
    calculating the changed cells in A_tag
    */

    for (n = 0; n < row; n++)
    {
        if (n != i && n != j)
        {

            A_tag[n][i] = (c * lp_matrix[n][i]) - (s * lp_matrix[n][j]);
            A_tag[n][j] = (c * lp_matrix[n][j]) + (s * lp_matrix[n][i]);
            /*A_tag is symmetrical*/
            A_tag[i][n] = A_tag[n][i];
            A_tag[j][n] = A_tag[n][j];
        }
    }

    A_tag[i][i] = (pow(c, 2) * lp_matrix[i][i]) + 
                    (pow(s, 2) * lp_matrix[j][j]) - (2 * s * c * lp_matrix[i][j]);
    A_tag[j][j] = (pow(s, 2) * lp_matrix[i][i]) + 
                    (pow(c, 2) * lp_matrix[j][j]) + (2 * s * c * lp_matrix[i][j]);

    A_tag[i][j] = 0.0;
    /*A_tag is symmetrical*/
    A_tag[j][i] = 0.0;

    return A_tag;
}

/*
calculates the sum of squares of all off-diagonal elements of the matrix
*/
double sos_off(double** sym_matrix, int row, int col)
{
    int i = 0;
    int j = 0;
    double sos = 0.0;
    for (i = 0; i < row; i++)
    {
        /*the input is a symmetrical matrix 
        computing only the upper triangle without diagonal*/
        for (j = i + 1; j < col; j++)
        {
            sos = sos + pow(sym_matrix[i][j], 2);
        }
    }
    /*from symmetry*/
    sos = sos*2;
    return sos;
}

/*
checks if there is a convergence by calculating:
off(A)^2 - off(A')^2 < = 0.001
*/
int is_convergence(double** lp_matrix, double** lp_matrix_tag, int row,
                                                 int col, int first_iter)
{
    /*lp_matrix_tag isn't ready before we start the jacobi algorithm*/
    if(first_iter == 1){
        return 0;
    }
    double sos_lp = sos_off(lp_matrix, row, col);
    double sos_lp_tag = sos_off(lp_matrix_tag, row, col);

    if (sos_lp - sos_lp_tag <= 0.001)
    {
        return 1;
    }
    return 0;
}

/*checks if the given matrix is diagonal by calculating 
the matrix's sum of square off diagonal*/
int is_diagonal(double** matrix, int row, int col)
{
    double matrix_sos = sos_off(matrix, row, col);
    if(matrix_sos == 0.0){
        return 1;
    }
    return 0;
}

/*
finds the largest off diagonal value, searches only the upper triangle 
since the matrix is symmetrical. 
returns an array of the value, and it's indices
*/
double* largest_off_digonal_value(double** lp_matrix, int row, int col)
{
    int i = 0;
    int j = 0;
    double max = 0.0;
    int index_r = 0;
    int index_c = 0;

    double* result = (double*)calloc(3, sizeof(double));
    assert(result != NULL && "An Error Has Occured");

    for (i = 0; i < row; i++)
    {
        for (j = i + 1; j < col; j++)
        {
            if (fabs(lp_matrix[i][j]) > max)
            {
                max = lp_matrix[i][j];
                index_c = j;
                index_r = i;
            }
        }
    }
    result[0] = max;
    result[1] = index_r;
    result[2] = index_c;
    return result;
}

/*
creating a rowxcol matrix filled with zeros
*/
double** zero_matrix(int row, int col)
{
    int i = 0;
    double** result;
    result = (double**)calloc(row, sizeof(double*));
    assert(result != NULL && "An Error Has Occured");

    for (i = 0; i < row; i++)
    {
        result[i] = (double*)calloc(col, sizeof(double));
        assert(result[i] != NULL && "An Error Has Occured");
    }

    return result;
}

/*
creating a rowxcol identity matrix
*/
double** identity_matrix(int row, int col)
{
    int i = 0;
    int j = 0;
    double** result;
    result = (double**)calloc(row, sizeof(double*));
    assert(result != NULL && "An Error Has Occured");

    for (i = 0; i < row; i++)
    {
        result[i] = (double*)calloc(col, sizeof(double));
        assert(result[i] != NULL && "An Error Has Occured");
    }

    /*
    filling the diagonal with 1's
    */
    for (i = 0; i < row; i++)
    {
        result[i][i] = 1;
    }

    return result;
}

/*
recieve a matrix and make a new transposed one
*/
double** transpose_matrix(double** matrix, int row, int col)
{
    double** result;
    int i = 0;
    int j = 0;

    result = zero_matrix(row, col);
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            result[i][j] = matrix[j][i];
        }
    }
    return result;
}

/*
matrix multiplication, only works for square matrixes
taken from: https://www.geeksforgeeks.org/c-program-multiply-two-matrices/
*/
double** matrix_multiplication(double** matrix1,
                               double** matrix2, int row, int col)
{
    double** result;
    int i = 0;
    int j = 0;
    int n = 0;
    result = zero_matrix(row, col);
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            for (n = 0; n < col; n++)
                result[i][j] += matrix1[i][n] * matrix2[n][j];
        }
    }
    return result;
}

double** three_matrix_multiplication(double** matrix1,
                                     double** matrix2, double** matrix3, int row, int col)
{
    double** result;
    double** tmp;
    tmp = matrix_multiplication(matrix1, matrix2, row, col);
    result = matrix_multiplication(tmp, matrix3, row, col);
    free_memory(tmp, row);
    return result;
}

void matrix_substraction(double** matrix1, double** matrix2, int row, int col)
{
    int i = 0;
    int j = 0;
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            matrix1[i][j] = matrix1[i][j] - matrix2[i][j];
        }
    }
}

/*
gets a diagonal matrix and raise each cell on the diagonal by the power of -0.5
*/
void manipulated_diagonal(double** diag_matrix, int row)
{
    int i = 0;
    for (i = 0; i < row; i++)
    {
        if (diag_matrix[i][i] != 0)
        {
            diag_matrix[i][i] = 1 / sqrt(diag_matrix[i][i]);
        }
    }
}

/*
print a matrix (row after row)
*/
void print_matrix(double** matrix, int row, int col)
{
    int i = 0;
    int j = 0;
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            if (j == col - 1)
            {
                printf("%.4f", matrix[i][j]);
            }
            else
            {
                printf("%.4f,", matrix[i][j]);
            }
        }
        printf("\n");
    }
}

/*
gets a list of lists (could be matrix or list of points)
free the inner lists and then the outer one
*/
void free_memory(double** matrix, int row)
{
    int i = 0;
    for (i = 0; i < row; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}