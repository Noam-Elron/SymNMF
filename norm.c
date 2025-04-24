#include <math.h>
#include "utils.h"
#include "sym.h"
#include "diagonal.h"

double **diagonal_matrix_multiplication(double **matrix, double **diagonal_matrix, int num_points, int multiplication_direction) {
    /* Helper function for norm to calculate result of multiplying matrix by a diagonal matrix, supports both left and right multiplication. Does not modify data passed in.
    Input:
        - double matrix[][]: Square matrix we are multiplying diagonal matrix by.
        - double diagonal_matrix[][]: Square diagonal matrix.
        - int num_points: Dimension of the square matrices (they must have the same dimensions).
        - int multiplication_direction: If 0, multiplies the diagonal matrix from the left (D * M), if 1 multiplies the diagonal matrix from the right(M * D).
    Returns:
        mxn result of the multiplication
    */
    int i, j;
    double **result_matrix = continuous_matrix_creation(num_points, num_points);

    for (i = 0; i < num_points; i++) {
        for (j = 0; j < num_points; j++) {
            result_matrix[i][j] = matrix[i][j] * diagonal_matrix[multiplication_direction ? j : i][multiplication_direction ? j : i];
        }
    }
    return result_matrix;
}



double **diagonal_matrix_exponentiation(double **diagonal_matrix, double exponent, int matrix_dimension) {
    /* Helper function for norm to calculate result of calculating exponent of diagonal matrix. Does not modify data passed in.
    Input: 
        - double diagonal matrix[][] D: Square matrix where each entry other than diagonal is zero, diagonal is free to be any value.
        - double exponent: Power we are raising diagonal matrix D by.
        - int matrix_dimension: Size of the square diagonal matrix D.
    Returns:
        2D Matrix D^(exponent)
    */
   int i;
   for (i = 0; i < matrix_dimension; i++) {
        diagonal_matrix[i][i] = pow(diagonal_matrix[i][i], exponent);
   }
   return diagonal_matrix;
}

double **norm_matrix(double **similarity_matrix, double **diagonal_matrix, int num_points) {
    /* Creates norm matrix as per project instructions
    Input: 
        - double Similarity Matrix[][]: Matrix where each entry corresponds to similarity between points as described in PDF.
        - double Diagonal Matrix[][]: Matrix where each entry other than diagonal is zero, i'th diagonal entry is sum of the i'th row in the similarity matrix.
    Returns:
        2D Norm Matrix W: W = D^(-1/2) * A * D^(-1/2)
    */
    double **norm_matrix;
    double **diag_exponentiated;
    double **temp_result;
    diag_exponentiated = matrix_deep_copy(diagonal_matrix, num_points, num_points);
    diag_exponentiated = diagonal_matrix_exponentiation(diag_exponentiated, -0.5, num_points);
    temp_result = diagonal_matrix_multiplication(similarity_matrix, diag_exponentiated, num_points, 0); /* Left Multiplication: D^(-1/2) * A  */
    norm_matrix = diagonal_matrix_multiplication(temp_result, diag_exponentiated, num_points, 1); /* Right Multiplication by prev result: (D^(-1/2) * A) * D^(-1/2) */

    free_continuous_matrix(diag_exponentiated);
    free_continuous_matrix(temp_result);

    return norm_matrix;
}