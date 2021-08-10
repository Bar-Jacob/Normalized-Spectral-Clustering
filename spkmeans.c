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
double** diagonal_degree_matrix(double** adj_matrix, int row, int col);
double** matrix_multiplication(double** matrix1, 
                                double** matrix2, int row, int col);
void manipulated_diagonal(double** dig_matrix, int row);
void matrix_substraction(double** matrix1, double** matrix2, int row, int col);
void free_memory(double** matrix, int row);
void wam_goal(double** points, int num_of_points, int dimension);

int main(int argc, char const *argv[])
{
    double** res;
    double** res2;

    double** dig;
    res = zero_matrix(3,3);
    res2 = zero_matrix(3,3);
    res[0][0] = 1;
    res[0][1] = 2;
    res[0][2] = 3;
    res[1][0] = 4;
    res[1][1] = 5;
    res[1][2] = 6;
    res[2][0] = 7;
    res[2][1] = 8;
    res[2][2] = 9;

    res2[0][0] = 1;
    res2[0][1] = 2;
    res2[0][2] = 3;
    res2[1][0] = 4;
    res2[1][1] = 5;
    res2[1][2] = 6;
    res2[2][0] = 7;
    res2[2][1] = 8;
    res2[2][2] = 9;
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
    wam_goal(res, 3, 3);
    
    return 0;
    
}

void wam_goal(double** points, int num_of_points, int dimension){
    double** result;
    result = adjacency_matrix(points, num_of_points, dimension);
    print_matrix(result, num_of_points, dimension);
    free_memory(result, num_of_points);
}

// void ddg_goal()
// void lnorm_goal
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

double** diagonal_degree_matrix(double** adj_matrix, int row, int col){
    double** result;
    int i = 0;
    int j = 0;
    result = zero_matrix(row, col);
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            result[i][i] += adj_matrix[i][j];
        }
    }
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
void manipulated_diagonal(double** dig_matrix, int row){
    int i = 0;
    for (i = 0; i < row; i++)
    {
        if(dig_matrix[i][i] != 0){
            dig_matrix[i][i] = 1/sqrt(dig_matrix[i][i]);
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