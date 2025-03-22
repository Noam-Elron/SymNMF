#include "utils.h"
#include "sym.h"

double matrix_row_sum(double row[], int num_points) {
    /* Helper function for diagonal matrix creation, calculates the sum of the i'th row in the similarity matrix.
    Input:
        - double row[]: i'th row in the similarity matrix we are currently calculating the sum of.
        - int num_points: Length of the i'th row
    Returns:
        Sum of the i'th row in the similarity matrix
    */
    int i;
    double sum = 0.0;
    for (i = 0; i < num_points; i++) {
        sum += row[i];
    }
    return sum;
}

double **diagonal_matrix(double **similarity_matrix, int num_points) {
    /* 
    Input: 
        - double Similarity Matrix[][]: Matrix where each entry corresponds to similarity between points as described in PDF.
        - int num_points: Number of points in Datapoints, this is also the size of the diagonal matrix as the matrix has a diagonal 
          of length n (each diagonal entry corresponds to a row in the similarity matrix)
    Returns:
        2D Square Diagonal Matrix, Diagonal entry i equals sum of row i in Similarity Matrix, all other entries are 0 (taken care of by continuous_matrix_creation_function who uses calloc)
    */
    int i;
    double **diagonal_matrix;
    
    diagonal_matrix = continuous_matrix_creation(num_points, num_points);
    for (i = 0; i < num_points; i++) {
        /* Calloc instantiates all elements to zero, therefore no need to set d_ij = 0 manually for i != j */
        diagonal_matrix[i][i] = matrix_row_sum(similarity_matrix[i], num_points);
    }
    return diagonal_matrix;
}