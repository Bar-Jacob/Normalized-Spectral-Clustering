//#include <spkmeans.h>
//transfer includes to the header
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

double** zero_matrix(int row, int col);
void print_matrix(double** matrix, int row, int col);
double** identity_matrix(int row, int col);
double** transpose_matrix(double** matrix, int row, int col);
double** adjacency_matrix(double** points, int row, int col);
double** diagonal_degree_matrix(double** points, int row, int col);
double** normalized_graph_laplacian(double** points, int row, int col);
double* largest_off_digonal_value(double** lp_matrix, int row, int col);
double** rotation_matrix(double** lp_matrix, int row, int col);
double** jacobi_calculations(double** lp_matrix, int row, 
                                                int col, double* sct_vals);

double** matrix_multiplication(double** matrix1, 
                                double** matrix2, int row, int col);
double** three_matrix_multiplication(double** matrix1, 
                                double** matrix2, double** matrix3, int row, int col);

void manipulated_diagonal(double** diag_matrix, int row);
void matrix_substraction(double** matrix1, double** matrix2, int row, int col);
void free_memory(double** matrix, int row);
void wam_goal(double** points, int num_of_points, int dimension);
void ddg_goal(double** points, int row, int col);
void lnorm_goal(double** points, int row, int col);
double* s_c_t_calculation(double** lp_matrix, int row, int col);


int main(int argc, char const *argv[])
{
    double** res;
    double** res2;
    double* high;

    double** dig;
    double** adj;
    double** i;

    i = identity_matrix(3,3);
    res = zero_matrix(3,3);
    res2 = zero_matrix(3,3);
    res[0][0] = 3;
    res[0][1] = 2;
    res[0][2] = 4;
    res[1][0] = 2;
    res[1][1] = 0;
    res[1][2] = 2;
    res[2][0] = 4;
    res[2][1] = 2;
    res[2][2] = 3;

    // res2[0][0] = 1;
    // res2[0][1] = 2;
    // res2[0][2] = 3;
    // res2[1][0] = 4;
    // res2[1][1] = 5;
    // res2[1][2] = 6;
    // res2[2][0] = 7;
    // res2[2][1] = 8;
    // res2[2][2] = 9;
    // for (int i = 0; i < 3; i++)
    // {
    //      for (int j = 0;j < 3; j++)
    //      {
    //          res[i][j] = i;
    //      }
        
    // }
    // res = adjacency_matrix(res, 3, 3);
    // dig = diagonal_degree_matrix(res, 3, 3);
    // manipulated_diagonal(dig, 3);
    // matrix_substraction(res2,res,3,3);
    // print_matrix(res2,3,3);
    // printf("\n");
    // print_matrix(res,3,3);
    // adj = adjacency_matrix(res, 3, 3);
    // dig = diagonal_degree_matrix(res, 3, 3);
    // manipulated_diagonal(dig, 3);
    // matrix_substraction(i,three_matrix_multiplication(dig, adj, dig, 3, 3),3,3);
    // print_matrix(i, 3, 3);
    // printf("\n");
    // res2 = normalized_graph_laplacian(res, 3, 3);
    // print_matrix(res2, 3, 3);
    high = s_c_t_calculation(res,3,3);
    res2 = jacobi_calculations(res,3,3, high);
    print_matrix(res2,3,3);

    return 0;
    
}

void wam_goal(double** points, int num_of_points, int dimension){
    double** result;
    result = adjacency_matrix(points, num_of_points, dimension);
    print_matrix(result, num_of_points, dimension);
    free_memory(result, num_of_points);
}

void ddg_goal(double** points, int row, int col){
    double** result;
    result = diagonal_degree_matrix(points, row, col);
    print_matrix(result, row, col);
    free_memory(result, row);
}

void lnorm_goal(double** points, int row, int col){
    double** result;
    result = normalized_graph_laplacian(points, row, col);
    print_matrix(result, row, col);
    free_memory(result, row);
}
// void jacobi_goal
// void spk_goal


/*
computing adjacency matrix
*/
double** adjacency_matrix(double** points, int num_of_points, int dimension){
    double** result;
    int i = 0;
    int j = 0;
    int d = 0;
    double weight = 0.0;
    result = zero_matrix(num_of_points,num_of_points);
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
                weight += pow((points[i][d] - points[j][d]),2); 
            }
            weight = sqrt(weight);
            weight = exp(-weight/2);
            /*
            filling the upper and lower triangle in parallel
            */
            result[i][j] = weight;
            result[j][i] = weight;
        } 
    }
    return result;
}

double** diagonal_degree_matrix(double** points, int row, int col){
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
double** normalized_graph_laplacian(double** points, int row, int col){
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

double** rotation_matrix(double** lp_matrix, int row, int col){
    double** result = identity_matrix(row, col);
    double* sct_val = s_c_t_calculation(lp_matrix, row, col);
    int i = (int) sct_val[3];
    int j = (int) sct_val[4]; 

    result[i][i] = sct_val[1];
    result[j][j] = sct_val[1];
    result[i][j] = sct_val[0];
    result[j][i] = -sct_val[0];
    
    free(sct_val);

    return result;
}

double* s_c_t_calculation(double** lp_matrix, int row, int col){
    double* pivot_arr = largest_off_digonal_value(lp_matrix, row, col);
    double* result_arr = (double*) calloc(5, sizeof(double));
    int i = (int) pivot_arr[1];
    int j = (int) pivot_arr[2];
    double pivot_max = pivot_arr[0];
    double teta = (lp_matrix[j][j] - lp_matrix[i][i])/(2*pivot_max);
    
    double t = 0.0;
    if(teta < 0){
        t = -1/(fabs(teta)+sqrt(pow(teta,2) + 1));
    }else{
        t = 1/(fabs(teta)+sqrt(pow(teta,2) + 1));
    }

    double c = 1/(sqrt(pow(t,2) + 1));
    double s = c*t;

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
                                                int col, double* sct_vals){
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
        if(n != i && n != j){
            A_tag[n][i] = (c*lp_matrix[n][i]) - (s*lp_matrix[n][j]);
            A_tag[n][j] = (c*lp_matrix[n][j]) + (s*lp_matrix[n][i]);
        }
    }

    A_tag[i][i] = (pow(c,2)*lp_matrix[i][i]) + (pow(s,2)*lp_matrix[j][j])
                                             - (2*s*c*lp_matrix[i][j]);

    A_tag[j][j] = (pow(s,2)*lp_matrix[i][i]) + (pow(c,2)*lp_matrix[j][j])
                                             + (2*s*c*lp_matrix[i][j]);
    
    A_tag[i][j] = 0.0;

    return A_tag;
}

/*
finds the largest off diagonal value, searches only the upper triangle 
since the matrix is symmetrical. 
returns an array of the value, and it's indices
*/
double* largest_off_digonal_value(double** lp_matrix, int row, int col){
    int i = 0;
    int j = 0;
    double max = 0.0;
    int index_r = 0;
    int index_c = 0;

    double* result = (double*) calloc(3, sizeof(double));
    assert(result != NULL && "An Error Has Occured");

    for (i = 0; i < row; i++)
    {
        for (j = i+1; j < col; j++)
        {
            if(fabs(lp_matrix[i][j]) > max){
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
double** zero_matrix(int row, int col){
    int i = 0;
    double** result;
    result = (double**) calloc(row, sizeof(double*));
    assert(result != NULL && "An Error Has Occured");

    for(i = 0; i < row; i++){
        result[i] = (double*) calloc(col, sizeof(double));
        assert(result[i] != NULL && "An Error Has Occured");
    }

    return result;
}

/*
creating a rowxcol identity matrix
*/
double** identity_matrix(int row, int col){
    int i = 0;
    int j = 0;
    double** result;
    result = (double**) calloc(row, sizeof(double*));
    assert(result != NULL && "An Error Has Occured");

    for(i = 0; i < row; i++){
        result[i] = (double*) calloc(col, sizeof(double));
        assert(result[i] != NULL && "An Error Has Occured");
    }

    /*
    filling the diagonal with 1's
    */
    for(i = 0; i < row; i++){
        result[i][i] = 1;
    }

    return result;
}

/*
recieve a matrix and make a new transposed one
*/
double** transpose_matrix(double** matrix, int row, int col){
    double** result;
    int i = 0;
    int j = 0;

    result = zero_matrix(row,col);
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
                                double** matrix2, int row, int col){
    double** result;
    int i = 0;
    int j = 0;
    int n = 0;
    result = zero_matrix(row,col);
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            for (n = 0; n < col; n++)
                result[i][j] += matrix1[i][n] * matrix2[n][j];
        }
    }
    return result;
}

double** three_matrix_multiplication(double** matrix1, 
                                double** matrix2, double** matrix3, int row, int col){
    double** result;
    double** tmp;
    tmp = matrix_multiplication(matrix1, matrix2, row, col);
    result = matrix_multiplication(tmp, matrix3, row, col);
    free_memory(tmp, row);
    return result;
}

void matrix_substraction(double** matrix1, double** matrix2, int row, int col){
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
void manipulated_diagonal(double** diag_matrix, int row){
    int i = 0;
    for (i = 0; i < row; i++)
    {
        if(diag_matrix[i][i] != 0){
            diag_matrix[i][i] = 1/sqrt(diag_matrix[i][i]);
        }
    }
}

/*
print a matrix (row after row)
*/
void print_matrix(double** matrix, int row, int col){
    int i = 0;
    int j = 0;
    for(i = 0; i < row; i++){
        for(j = 0; j < col; j++){
            if(j == col-1){
                printf("%.4f", matrix[i][j]);
            }else{
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
void free_memory(double** matrix, int row){
    int i = 0;
    for(i = 0; i < row; i++){
        free(matrix[i]);
    }
    free(matrix);
}