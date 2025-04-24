#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "sym.h"
#include "diagonal.h"
#include "norm.h"
#define beta 0.5
#define epsilon 0.0001 
#define max_iter 300


double **update_H(double **prev_H, double **W, int n, int k) {
    /* Updates H to next iteration as per project instructions
    Input: 
        - double prev_H[][]: Previous iteration of H we are trying to update.
        - double W[][]: Norm matrix we are using to calculate next iteration of H.
        - int n: Size of norm matrix, number of rows in H.
        - int k: Number of columns in H.
    Returns:
        Next iteration of H
    */
    double **w_h_mult;
    double **h_t;
    double **h_h_t_mult;
    double **h_h_t_h_mult;
    double **next_H;
    int i;
    int j;
    w_h_mult = matrix_multiplication(W, prev_H, n, n, k);
    h_t = matrix_transpose(prev_H, k, n);
    h_h_t_mult = matrix_multiplication(prev_H, h_t, n, k, n);
    h_h_t_h_mult = matrix_multiplication(h_h_t_mult, prev_H, n, n, k);
    
    next_H = continuous_matrix_creation(n, k);
    for (i = 0; i < n; i++) {
        for (j = 0; j < k; j++) {
            next_H[i][j] = prev_H[i][j] * (1 - beta + beta*(w_h_mult[i][j]/h_h_t_h_mult[i][j]));
        }
    }

    free_continuous_matrix(w_h_mult);
    free_continuous_matrix(h_t);
    free_continuous_matrix(h_h_t_mult);
    free_continuous_matrix(h_h_t_h_mult);
    return next_H;
}   

double frobenius_norm(double **matrix, int n, int k) {
    /* Calculates frobenius norm of matrix.
    Input: 
        - double matrix[][]: Matrix we are trying to calculate frobenius norm of.
        - int n: Number of rows in matrix.
        - int k: Number of columns in matrix.
    Returns:
        Frobenius norm of given matrix
    */
    double **transpose;
    double **m_t_mult;
    double trace;
    transpose = matrix_transpose(matrix, n, k);
    m_t_mult = matrix_multiplication(transpose, matrix, k, n, k);
    trace = matrix_trace(m_t_mult, k);
    return sqrt(trace);
}


double **converge_H(double **initial_H, double **W, int n, int k) {
    /* Continuously updates H until either convergence or until reaching max iterations, as per project instructions
    Input: 
        - double Initial_H[][]: Initial H matrix we received from Python.
        - double W[][]: Norm matrix.
        - int n: Size of norm matrix, number of rows in H.
        - int k: Number of columns in H.
    Returns:
        Final iteration of H.
    */
    double **prev_H;
    double **cur_H;
    double **distance_matrix;
    double frobenius_distance;
    int iteration = 0;
    
    prev_H = matrix_deep_copy(initial_H, n, k);
    while (iteration < max_iter) {
        iteration++;
        cur_H = update_H(prev_H, W, n, k);
        distance_matrix = matrix_subtraction(cur_H, prev_H, n, k);
        frobenius_distance = frobenius_norm(distance_matrix, n, k);
        free_continuous_matrix(prev_H);
        if (frobenius_distance < epsilon) {
            break;
        }
        prev_H = cur_H;
    }

    return cur_H;
}






typedef struct datapoints_wrapper {
    double **datapoints;
    int dimension;
} datapoints_wrapper;


datapoints_wrapper* initialize_data(const char *filename) {
    /* Initializes a datapoints wrapper in order to read datapoints from file
    Input: 
        - char[] filename: string filepath to .txt file containing datapoints
    Returns:
        Initialized datapoint wrapper.
    */
    FILE *file = fopen(filename, "r");
    int rows = 0;
    char line[1023];
    datapoints_wrapper *wrapper;

    if (file) {
        while (fgets(line, sizeof(line), file)) {
            if (line[0] != '\n') rows++;
        }
        rewind(file); 
    }
    if (!file || rows == 0) return NULL;

    wrapper = malloc(sizeof(datapoints_wrapper));
    wrapper->datapoints = malloc(rows * sizeof(double*));
    wrapper->dimension = rows;
    fclose(file);
    return wrapper;
}

void populate_data(datapoints_wrapper *wrapper, const char *filename) {
    /* Reads datapoints from file into initialized datapoints wrapper
    Input: 
        - datapoints_wrapper *wrapper: initialized datapoints wrapper.
        - char[] filename: string filepath to .txt file containing datapoints
    */
    FILE *file = fopen(filename, "r");
    char *token;
    char line[1023];
    int row = 0;
    int col;
    int dim = wrapper->dimension;

    if (!file) {
        return;
    }    
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '\n') continue;
        wrapper->datapoints[row] = malloc(dim * sizeof(double));
        token = strtok(line, ",");
        for (col = 0; col < dim && token != NULL; ++col) {
            wrapper->datapoints[row][col] = atof(token);
            token = strtok(NULL, ",");
        }
        row++;
    }
    fclose(file);
    return;
}


void sym(datapoints_wrapper *datapoints) {
    /* Wrapper function to calculate similarity matrix as per project instructions
    Input: 
        - datapoints_wrapper *datapoints: datapoints wrapper.
    */
    double **sym_matrix;
    int n = datapoints->dimension;
    sym_matrix = similarity_matrix(datapoints->datapoints, n, n);
    print_matrix(sym_matrix, n, n);
    free_continuous_matrix(sym_matrix);
}

void ddg(datapoints_wrapper *datapoints) {
    /* Wrapper function to calculate diagonal matrix as per project instructions
    Input: 
        - datapoints_wrapper *datapoints: datapoints wrapper.
    */
    double **sym_matrix;
    double **diag_matrix;
    int n = datapoints->dimension;
    sym_matrix = similarity_matrix(datapoints->datapoints, n, n);
    diag_matrix = diagonal_matrix(sym_matrix, n, n);
    print_matrix(diag_matrix, n, n);
    free_continuous_matrix(sym_matrix);
    free_continuous_matrix(diag_matrix);
}

void norm(datapoints_wrapper *datapoints) {
    /* Wrapper function to calculate norm matrix as per project instructions
    Input: 
        - datapoints_wrapper *datapoints: datapoints wrapper.
    */
    double **sym_matrix;
    double **diag_matrix;
    double **normal_matrix;
    int n = datapoints->dimension;
    sym_matrix = similarity_matrix(datapoints->datapoints, n, n);
    diag_matrix = diagonal_matrix(sym_matrix, n, n);
    normal_matrix = norm_matrix(sym_matrix, diag_matrix, n);
    print_matrix(normal_matrix, n, n);
    free_continuous_matrix(sym_matrix);
    free_continuous_matrix(diag_matrix);
    free_continuous_matrix(normal_matrix);
}


int main(int argc, char **argv) {
    /* Main function for SymNMF in C, handles user input and calling appropriate wrapper functions.
    Input:
        - int argc: number of passed in user arguments (must be exactly 3 in order for program to run)
        - char **argv: user arguments (must be in the form (c_filename, goal, filepath))
    */
    datapoints_wrapper *datapoints;
    
    char *goals[] = {"sym", "ddg", "norm"};
    if (argc != 3) {
        printf("An Error has Occurred\n");
        exit(EXIT_FAILURE);
    }
    datapoints = initialize_data(argv[2]);
    populate_data(datapoints, argv[2]);
    
    if (strcmp(goals[0], argv[1]) == 0) {
        sym(datapoints);
    }
    else if (strcmp(goals[1], argv[1]) == 0) {
        ddg(datapoints);
    }
    else if (strcmp(goals[2], argv[1]) == 0) {
        norm(datapoints);
    }
    else {
        free_matrix(datapoints->datapoints, datapoints->dimension);
        free(datapoints);
        printf("An Error has Occurred\n");
        exit(EXIT_FAILURE);
    }
    free_matrix(datapoints->datapoints, datapoints->dimension);
    free(datapoints);
    exit(EXIT_SUCCESS);
} 



