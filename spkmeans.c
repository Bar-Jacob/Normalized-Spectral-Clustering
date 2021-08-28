#include "spkmeans.h"

int main(int argc, char const *argv[])
{
    FILE* file;
    int num_of_points = 0;
    int first_point = 1;
    int i = 0;
    int dimension = 1;
    char* token;
    double result;
    int dimension_cnt = 0;
    int k = atoi(argv[1]);
    char goal = argv[2][0];
    file = fopen(argv[3], "r");
    assert(file != NULL && "An Error Has Occurred");
    /*we have a limit of 10 features, 9 commas,
    and each number has max of 4 digits before decimal point and after it*/
    char* str_point = (char *)calloc(109, sizeof(char));
    assert(str_point != NULL && "An Error Has Occurred");
    double** data_points = (double **)calloc(50, sizeof(double *));
    assert(data_points != NULL && "An Error Has Occurred");

    while (fgets(str_point, 109, file) != NULL)
    {
        if (first_point == 1)
        {
            for (i = 0; i < 109; i++)
            {
                /*calculating the dimension*/
                if (str_point[i] == ',')
                {
                    dimension++;
                }
            }
            /*memory for the points*/
            for (i = 0; i < 1000; i++)
            {
                double* point = (double*)malloc(dimension * sizeof(double));
                assert(point != NULL && "An Error Has Occurred");
                data_points[i] = point;
            }
        }
        token = strtok(str_point, ",");

        /* walk through other tokens */
        while (token != NULL)
        {
            result = strtod(token, NULL);
            data_points[num_of_points][dimension_cnt] = result;
            dimension_cnt++;
            token = strtok(NULL, ",");
        }
        dimension_cnt = 0;
        first_point = 0;
        num_of_points++;
    }
    fclose(file);
    data_points = (double**)realloc(data_points, (num_of_points) * sizeof(double *));
    assert(data_points != NULL && "An Error Has Occurred");

    switch (goal)
    {
    case 'w':
        wam_goal(data_points, num_of_points, dimension);
        break;
    case 'd':
        ddg_goal(data_points, num_of_points, dimension);
        break;
    case 'l':
        lnorm_goal(data_points, num_of_points, dimension);
        break;
    case 'j':
        jacobi_goal(data_points, num_of_points);
        break;
    case 's':
        spk_goal(data_points, num_of_points, dimension, k);
        break;
    }
    return 0;
}

void wam_goal(double** points, int num_of_points, int dimension)
{
    double** result;
    result = adjacency_matrix(points, num_of_points, dimension);
    print_matrix(result, num_of_points, num_of_points);
    free_memory(result, num_of_points);
    free_memory(points, num_of_points);
}

void ddg_goal(double** points, int row, int dimension)
{
    double** result;
    result = diagonal_degree_matrix(points, row, dimension);
    print_matrix(result, row, row);
    free_memory(result, row);
    free_memory(points, row);
}

void lnorm_goal(double** points, int row, int dimension)
{
    double** result;
    result = normalized_graph_laplacian(points, row, dimension);
    print_matrix(result, row, row);
    free_memory(result, row);
    free_memory(points, row);
}

void jacobi_goal(double** sym_matrix, int row)
{
    Eigenvector* result;
    int i = 0;

    result = jacobi(sym_matrix, row);
    for (i = 0; i < row; i++)
    {
        if (i == row - 1)
        {
            printf("%.4f", result[i].value);
        }
        else
        {
            printf("%.4f,", result[i].value);
        }

    }
    printf("\n");

    for (i = 0; i < row; i++)
    {
        print_array(result[i].vector, row);
        free(result[i].vector);
        printf("\n");
    }
    free(result);
}

void spk_goal(double** points, int row, int col, int k)
{
    double** U;
    double** laplacian = normalized_graph_laplacian(points, row, col);
    Eigenvector* eignvectors = jacobi(laplacian, row);

    if(k == 0){
        k = eigengap_heuristic(eignvectors, row);
    }else{
        merge_sort(eignvectors, 0, row-1);

    }
    U = creating_U(eignvectors, k, row);
    renormalizing_U(U, k, row);
    kmeans(U, k, k, row);
    free_memory(U, row);
    free_memory(points, row);
}

PyObject* spk_goal_python(double** points, int row, int col, int k)
{
    double** U;
    double** laplacian = normalized_graph_laplacian(points, row, col);
    PyObject* points_py;
    Eigenvector* eignvectors = jacobi(laplacian, row);

    if(k == 0){
        k = eigengap_heuristic(eignvectors, row);
    }else{
        merge_sort(eignvectors, 0, row-1);

    }
    U = creating_U(eignvectors, k, row);
    renormalizing_U(U, k, row);
    points_py = cToPyObject(U, k, row, k);
    free_memory(U, row);
    free_memory(points, row);
    return points_py;
}

/*takes the first k eignvectors (after being sorted by their values)
put the k vectors as columns*/
double** creating_U(Eigenvector* eignvector, int k, int row)
{
    double** U = zero_matrix(row, k);
    int i = 0;
    int j = 0;

    for (i = 0; i < k; i++)
    {
        for (j = 0; j < row; j++)
        {
            U[j][i] = eignvector[i].vector[j];
        }
    }
    return U;
}

/*renormalizing each of U’s rows to have unit length*/
void renormalizing_U(double** U, int k, int row)
{
    int i = 0;
    int j  = 0;
    for (i = 0; i < row; i++)
    {
        /*computing sum of squares of each row*/
        double sos_row = 0.0;
        for (j = 0; j < k; j++)
        {
            sos_row += pow(U[i][j],2);
        }
        sos_row = sqrt(sos_row);

        /*updating U[i][j] to be renormalized*/
        for (j = 0; j < k; j++)
        {
            U[i][j] = U[i][j]/sos_row;
        }
    }
     
}

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
            weight = 0.0;
        }
    }
    return result;
}

double** diagonal_degree_matrix(double** points, int row, int dimension)
{
    double** result;
    double** adj_matrix;
    int i = 0;
    int j = 0;
    result = zero_matrix(row, row);
    adj_matrix = adjacency_matrix(points, row, dimension);
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < row; j++)
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
double** normalized_graph_laplacian(double** points, int row, int dimension)
{
    double** adj_matrix;
    double** diag_matrix;
    double** multip_matrix;
    double** result;

    adj_matrix = adjacency_matrix(points, row, dimension);
    diag_matrix = diagonal_degree_matrix(points, row, dimension);
    /*
    changes the matrix itself, since we give a pointer to the matrix 
    */
    manipulated_diagonal(diag_matrix, row);
    multip_matrix = three_matrix_multiplication(diag_matrix, adj_matrix,
                                                diag_matrix, row, row);
    free_memory(adj_matrix, row);
    free_memory(diag_matrix, row);
    result = identity_matrix(row, row);
    /*
    changes the first matrix given, so the result will be on the result matrix
    */
    matrix_substraction(result, multip_matrix, row, row);
    free_memory(multip_matrix, row);

    return result;
}

/*
The jacobi eigenvalue algorithm
*/
Eigenvector* jacobi(double** lp_matrix, int row)
{
    double** lp_matrix_tag;
    double** eigenvectors_matrix;
    double** tmp_matrix;
    double** rotate_matrix;
    double* sct_vals;
    Eigenvector* result;
    int cnt = 0;
    int converged = 0;

    eigenvectors_matrix = identity_matrix(row, row);

    /*
    The algorithm will stop after 100 iterations or if there is a convergence,
    whatever happens first
    */
    while (cnt != 100 && converged == 0)
    {
        sct_vals = s_c_t_calculation(lp_matrix, row);
        rotate_matrix = rotation_matrix(lp_matrix, row, sct_vals);
        lp_matrix_tag = jacobi_calculations(lp_matrix, row, sct_vals);

        /*points to the last updated eigenvectors matrix*/
        tmp_matrix = eigenvectors_matrix;

        /*eigenvectors_matrix will now point to a new updated eigenvectors matrix*/
        eigenvectors_matrix = matrix_multiplication(eigenvectors_matrix, rotate_matrix, row, row);

        /*release tmp_matrix*/
        free_memory(tmp_matrix, row);
        tmp_matrix = lp_matrix;

        if (is_convergence(lp_matrix, lp_matrix_tag, row) == 1)
        {
            free_memory(tmp_matrix, row);
            converged = 1;
        }
        else
        {
            lp_matrix = lp_matrix_tag;
            free_memory(tmp_matrix, row);
        }
        cnt++;
        free(sct_vals);
        free_memory(rotate_matrix, row);
    }


    /*returns an array with two matrixes - eigenvectors and eigenvalues */
    result = creating_eigenvector_array(eigenvectors_matrix, lp_matrix_tag, row);
    free_memory(eigenvectors_matrix, row);
    free_memory(lp_matrix_tag, row);
    return result;
}

/*
creating an array of eigenvector structs
*/
Eigenvector* creating_eigenvector_array(double** eigenvectors_matrix, 
                                    double** eigenvalues_matrix, int row)
{
    Eigenvector* result = (Eigenvector*)malloc(row * (sizeof(struct Eigenvector)));
    assert(result != NULL && "An Error Has Occured");
    int i = 0;
    int m = 0;
    for (i = 0; i < row; i++)
    {
        Eigenvector tmp_struct;
        tmp_struct.value = eigenvalues_matrix[i][i];
        tmp_struct.vector = (double*)malloc(row * sizeof(double));
        assert(tmp_struct.vector != NULL && "An Error Has Occured");
        for (m = 0; m < row; m++)
        {
            tmp_struct.vector[m] = eigenvectors_matrix[m][i];
        }
        result[i] = tmp_struct;
    }
    return result;
}

/*
In order to determine the number of clusters k, we will use eigengap heuristic
*/
int eigengap_heuristic(Eigenvector* eigenvector, int row)
{
    double max = 0;
    int i = 0;
    int index = 0;
    int lim = row / 2;
    /*sorting the eigenvector structs by comparing the eigenvalues*/
    merge_sort(eigenvector, 0, row-1);
    for (i = 0; i < lim; i++)
    {
        if (fabs(eigenvector[i].value - eigenvector[i + 1].value) > max)
        {
            max = fabs(eigenvector[i].value - eigenvector[i + 1].value);
            index = i;
        }
    }
    return index + 1;
}

/*
builds the rotation matrix
*/
double** rotation_matrix(double** lp_matrix, int row, double* sct_val)
{
    double** result = identity_matrix(row, row);
    int i = (int)sct_val[3];
    int j = (int)sct_val[4];

    result[i][i] = sct_val[1];
    result[j][j] = sct_val[1];
    result[i][j] = sct_val[0];
    result[j][i] = -sct_val[0];

    return result;
}
/*
calculates s, c, t
*/
double* s_c_t_calculation(double** lp_matrix, int row)
{
    double* pivot_arr = largest_off_digonal_value(lp_matrix, row);
    double* result_arr = (double *)calloc(5, sizeof(double));
    int i = (int)pivot_arr[1];
    int j = (int)pivot_arr[2];
    double pivot_max = pivot_arr[0];
    double teta = (lp_matrix[j][j] - lp_matrix[i][i]) / (2 * pivot_max);
    double tmp_calc = fabs(teta) + sqrt(pow(teta, 2) + 1);
    
    double t = 0.0;
    if (teta < 0)
    {
        t = -1 / (tmp_calc);
    }
    else
    {
        t = 1 / (tmp_calc);
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
                             double* sct_vals)
{
    double** A_tag = zero_matrix(row, row);
    int n = 0;
    int m = 0;
    int i = (int) sct_vals[3];
    int j = (int) sct_vals[4];
    double c = sct_vals[1];
    double s = sct_vals[0];
    /*
    copy the laplacian matrix, in order to use the input for the calculations ahead
    */
    for (n = 0; n < row; n++)
    {
        for (m = 0; m < row; m++)
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
double sos_off(double** sym_matrix, int row)
{
    int i = 0;
    int j = 0;
    double sos = 0.0;
    for (i = 0; i < row; i++)
    {
        /*the input is a symmetrical matrix 
        computing only the upper triangle without diagonal*/
        for (j = i + 1; j < row; j++)
        {
            sos += pow(sym_matrix[i][j], 2);
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
int is_convergence(double** lp_matrix, double** lp_matrix_tag, int row)
{
    double sos_lp = sos_off(lp_matrix, row);
    double sos_lp_tag = sos_off(lp_matrix_tag, row);
    double epsilon = pow(10, -15);

    if (sos_lp - sos_lp_tag <= epsilon || (sos_lp_tag == 0.0))
    {
        return 1;
    }
    return 0;
}

/*
finds the largest off diagonal value, searches only the upper triangle 
since the matrix is symmetrical. 
returns an array of the value, and it's indices
*/
double* largest_off_digonal_value(double** lp_matrix, int row)
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
        for (j = i + 1; j < row; j++)
        {
            if (fabs(lp_matrix[i][j]) > fabs(max))
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
    result = zero_matrix(row, row);
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
                if (matrix[i][j] < 0 && matrix[i][j] > -0.00001)
                {
                    matrix[i][j] = matrix[i][j] * -1;
                }
                printf("%.4f", matrix[i][j]);
            }
            else
            {
                if (matrix[i][j] < 0 && matrix[i][j] > -0.00004)
                {
                    matrix[i][j] = matrix[i][j] * -1;
                }
                printf("%.4f,", matrix[i][j]);
            }
        }
        printf("\n");
    }
}

void print_array(double* array, int row)
{
    int i = 0;
    for (i = 0; i < row; i++)
    {
        if (i == row - 1)
        {
        if (array[i] < 0 && array[i] > -0.00004)
        {
            array[i] = array[i] * -1;
        }
            printf("%.4f", array[i]);
        }
        else
        {
        if (array[i] < 0 && array[i] > -0.00004)
        {
            array[i] = array[i] * -1;
        }
            printf("%.4f,", array[i]);
        }
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

/*https://www.geeksforgeeks.org/merge-sort/*/
// Merges two subarrays of arr[].
// First subarray is arr[l..m]
// Second subarray is arr[m+1..r]
void merge(Eigenvector* eigenvector, Eigenvector* R, Eigenvector* L, int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++)
        L[i] = eigenvector[l + i];
    for (j = 0; j < n2; j++)
        R[j] = eigenvector[m + 1 + j];

    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2)
    {
        if (L[i].value <= R[j].value)
        {
            eigenvector[k] = L[i];
            i++;
        }
        else
        {
            eigenvector[k] = R[j];
            j++;
        }
        k++;
    }

    /* Copy the remaining elements of L[], if there
    are any */
    while (i < n1)
    {
        eigenvector[k] = L[i];
        i++;
        k++;
    }

    /* Copy the remaining elements of R[], if there
    are any */
    while (j < n2)
    {
        eigenvector[k] = R[j];
        j++;
        k++;
    }
}

/* l is for left index and r is right index of the
sub-array of arr to be sorted */
void merge_sort(Eigenvector* eigenvector, int l, int r)
{
    if (l < r)
    {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        int m = l + (r - l) / 2;
        int n1 = m - l + 1;
        int n2 = r - m;

        /* create temp arrays */
        Eigenvector* L = (Eigenvector*)malloc(n1 * sizeof(struct Eigenvector));
        Eigenvector* R = (Eigenvector*)malloc(n2 * sizeof(struct Eigenvector));

        // Sort first and second halves
        merge_sort(eigenvector, l, m);
        merge_sort(eigenvector, m + 1, r);

        merge(eigenvector, R, L, l, m, r);

        free(L);
        free(R);
    }
}

/*
converting all data points from python into data points in C
*/
double** convert_python_to_c(PyObject* data_points_p, int dimension, int num_of_points)
{
    int cnt = 0;
    int i = 0;
    int j = 0;
    double** data_points = (double**)calloc(num_of_points, sizeof(*data_points));
    assert(data_points != NULL && "An Error Has Occured");

    for (i = 0; i < num_of_points; i++)
    {
        data_points[i] = (double *)calloc(dimension, sizeof(*data_points[i]));
        assert(data_points[i] != NULL && "An Error Has Occured");
        for (j = 0; j < dimension; j++)
        {
            data_points[i][j] = PyFloat_AsDouble(PyList_GetItem(data_points_p, cnt));
            cnt++;
        }
    }
    return data_points;
}

/*
after finishing running kmeans algorithm we want to return the results to python 
converting types from C to python
*/
/****add int k to input**/
PyObject* cToPyObject(double** T, int dimension, int num_of_points, int k)
{
    PyObject* points_py;
    int i = 0;
    int j = 0;
    PyObject* value;

    /***adding another cell for k***/
    points_py = PyList_New(num_of_points+1);
    for (i = 0; i < num_of_points; i++)
    {
        PyObject* curr_vector;
        curr_vector = PyList_New(dimension);
        for (j = 0; j < dimension; j++)
        {
            value = Py_BuildValue("d", T[i][j]);
            PyList_SetItem(curr_vector, j, value);
        }
        /*
        adding the PyObject centroid to the PyList clusters
        */
        PyList_SetItem(points_py, i, curr_vector);
    }
    /***getting k as last var***/
    value = Py_BuildValue("i", k);
    PyList_SetItem(points_py, num_of_points, value);
    return points_py;
}