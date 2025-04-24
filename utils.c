#include <stdlib.h>
#include <stdio.h>

double **continuous_matrix_creation(int m, int n) {
    /* Creates a continuous matrix via method shown in class 
    Input:
        - int m: Number of rows in matrix
        - int n: Number of columns in matrix
    Returns:
        - Continuous mxn matrix, all elements are zero instantiated by default due to use of calloc
    */
    int i;
    double *flattened_matrix;
    double **matrix;

    flattened_matrix = calloc(m * n, sizeof(double));
    matrix = calloc(m, sizeof(double *));
    for (i = 0; i < m; i++) {
        matrix[i] = flattened_matrix + i * n;
    }

    return matrix;
}

double **matrix_deep_copy(double **matrix_to_copy, int m, int n) {
    /* Deep copies a 2D matrix in order to preserve immutability
    Input:
        - double matrix_to_copy[][]: Matrix whose values we're copying, by copying just values instead of pointers we are essentially deep-copying it.
        - int m: Number of rows in matrix we're copying
        - int n: Number of columns in matrix we're copying
    Returns:
        - Deep copy of given matrix
    */
    double **copy = continuous_matrix_creation(m, n);
    int i;
    int j;

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            copy[i][j] = matrix_to_copy[i][j];
        }
    }
    return copy;
}

double **matrix_scaling(double **matrix, double constant, int m, int n) {
    /* Creates new matrix by multiplying (coordinate wise) given matrix by a scalar
    Input:
        - double matrix[][]: Matrix we are scaling.
        - double constant: Constant to multiply matrix entries by.
        - int m: Number of rows in the matrix.
        - int n: Number of columns in the matrix.
    Returns:
        mxn result of matrix scaling.
    */
    int i;
    int j;
    double **result_matrix = continuous_matrix_creation(m, n);

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            result_matrix[i][j] = matrix[i][j] * constant;
        }
    }
    return result_matrix;
}

double **matrix_multiplication(double **matrix, double **other_matrix, int m, int s, int n) {
    /* Multiplies two matrices together.
    Input:
        - double matrix[][]: Left matrix we are multiplying by.
        - double other_matrix[][]: Right matrix we are multiplying by.
        - int m: Number of rows in left matrix.
        - int s: Number of columns in left matrix / Number of rows in right matrix.
        - int n: Number of columns in right matrix.
    Returns:
        mxn result of multiplying the matrices.
    */
    int i;
    int j;
    int k;
    double **result_matrix;
    result_matrix = continuous_matrix_creation(m, n);
    
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < s; k++) {
                result_matrix[i][j] += matrix[i][k] * other_matrix[k][j];
            }
        }
    }
    return result_matrix;
}

double **matrix_transpose(double **matrix, int m, int n) {
    /* Transposes given matrix (immutable operation, returns new matrix and doesn't modify old one)
    Input:
        - double matrix[][]: Matrix we calculate transpose for.
        - int m: Number of rows in the matrix.
        - int n: Number of columns in the matrix.
    Returns:
        mxn result of transpose on the matrix.
    */
    int i;
    int j;
    double **result_matrix;
    result_matrix = continuous_matrix_creation(n, m);

    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            result_matrix[i][j] = matrix[j][i];
        }
    }
    return result_matrix;
}

double matrix_trace(double **matrix, int n) {
    /* Calculates trace of square matrix
    Input:
        - double matrix[][]: Matrix we calculate trace for.
        - int n: Size of matrix.
    Returns:
        Trace of given matrix.
    */
    int i;
    double trace = 0;

    for (i = 0; i < n; i++) {
        trace += matrix[i][i];
    }
    return trace;
}

void free_matrix(double **matrix, int num_rows) {
    /* Frees up matrix memory by freeing every pointer to a row and then freeing the pointer to the array of row pointers.
    Input:
        - double matrix[][]: Matrix whose memory we are freeing
        - int num_rows: Number of rows in the matrix
    */
    int i;
    if (matrix == NULL) return;  /* Safety check */

    for (i = 0; i < num_rows; i++) {
        free(matrix[i]);  /* Free each row */
    }
    free(matrix);  /* Free the row pointer array */
}

void free_continuous_matrix(double **continuous_matrix) {
    /* Frees up continuous matrix memory by freeing the pointer to the flattened array, at continuous_matrix[0] and 
       then freeing continuous_matrix: the array(pointer) of "pseudo-row" pointers(which are part of the flattened array allocation and therefore dont need to be freed)
       which point to the appropriate row positions in the flattened array.
    Input:
        - double matrix[][]: Matrix whose memory we are freeing
        - int num_rows: Number of rows in the matrix
    */
    if (continuous_matrix == NULL) return;
    free(continuous_matrix[0]);  /* Free the flattened array */
    free(continuous_matrix);     /* Free the array of row pointers */
}

void print_matrix(double **matrix, int m, int n) {
    /* Prints matrix as per project specifications.
    Input:
        - double **matrix: The matrix to be printed
        - int m: Number of rows
        - int n: Number of columns
    */
    int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("%.4f", matrix[i][j]);  /* Print each element as an integer (no decimal) */
            if (j < n-1) {
                putchar(',');
            }
        }
        printf("\n");
    }
}

double **matrix_subtraction(double **matrix, double **other_matrix, int m, int n) {
    /* Subtracts one matrix from another (coordinate wise) (immutable operation, returns new matrix and doesn't modify old ones)
    Input:
        - double matrix[][]: Matrix we are subtracting form.
        - double other_matrix[][]: Matrix we are using to subtract.
        - int m: Number of rows in both matrices.
        - int n: Number of columns in both matrices.
    Returns:
        mxn result of subtracting other_matrix from matrix
    */
    double **subtracted_matrix;
    int i;
    int j;
    
    subtracted_matrix = continuous_matrix_creation(m, n);
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            subtracted_matrix[i][j] = matrix[i][j] - other_matrix[i][j];
        }
    }
    return subtracted_matrix;
}


