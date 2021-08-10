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