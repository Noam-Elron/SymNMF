#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "sym.h"
#include "diagonal.h"
#include "norm.h"
#define beta 0.5
#define epsilon 1e-4
#define max_iter 300

void free_update_H_matrices(double **w_h_mult, double **h_t, double **h_h_t_mult, double **h_h_t_h_mult){
    /* Frees up H matrices for convenience. Matrices can be NULL as free_continuous_matrix which is used here handles it. 
    Input: 
        - double prev_H[][]: Matrix 1 we are freeing.
        - double h_t[][]: Matrix 2 we are freeing.
        - double h_h_t_mult[][]: Matrix 3 we are freeing.
        - double h_h_t_h_mult[][]: Matrix 4 we are freeing.
    */
    free_continuous_matrix(w_h_mult);
    free_continuous_matrix(h_t);
    free_continuous_matrix(h_h_t_mult);
    free_continuous_matrix(h_h_t_h_mult);
    return;
}

double **update_H(double **prev_H, double **W, int n, int k) {
    /* Updates H to next iteration as per project instructions. Returns NULL on error.
    Input: 
        - double prev_H[][]: Previous iteration of H we are trying to update.
        - double W[][]: Norm matrix we are using to calculate next iteration of H.
        - int n: Size of norm matrix, number of rows in H.
        - int k: Number of columns in H.
    Returns:
        Next iteration of H. 
    */
    double **w_h_mult, **h_t, **h_h_t_mult, **h_h_t_h_mult, **next_H;
    int i, j;
    w_h_mult = matrix_multiplication(W, prev_H, n, n, k);
    h_t = matrix_transpose(prev_H, n, k);
    if (w_h_mult == NULL || h_t == NULL) {
        free_update_H_matrices(w_h_mult, h_t, NULL, NULL);
        return NULL;
    }
    h_h_t_mult = matrix_multiplication(prev_H, h_t, n, k, n);
    if (h_h_t_mult == NULL) {
        free_update_H_matrices(w_h_mult, h_t, h_h_t_mult, NULL);
        return NULL;
    }
    h_h_t_h_mult = matrix_multiplication(h_h_t_mult, prev_H, n, n, k);
    next_H = continuous_matrix_creation(n, k);
    if (h_h_t_h_mult == NULL || next_H == NULL) {
        free_update_H_matrices(w_h_mult, h_t, h_h_t_mult, h_h_t_h_mult);
        free_continuous_matrix(next_H);
        return NULL;
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < k; j++) {
            next_H[i][j] = prev_H[i][j] * (1 - beta + beta*(w_h_mult[i][j]/h_h_t_h_mult[i][j]));
        }
    }
    free_update_H_matrices(w_h_mult, h_t, h_h_t_mult, h_h_t_h_mult);
    return next_H;
}   



double frobenius_norm_squared(double **matrix, int n, int k) {
    /* Calculates squared frobenius norm of matrix. Returns -1.0 on error.
    Input: 
        - double matrix[][]: Matrix we are trying to calculate frobenius norm of. 
        - int n: Number of rows in matrix.
        - int k: Number of columns in matrix.
    Returns:
        Frobenius norm of given matrix. 
    */
    double **transpose;
    double **m_t_mult;
    double trace;
    transpose = matrix_transpose(matrix, n, k);
    if (transpose == NULL) {
        return -1.0;
    }
    m_t_mult = matrix_multiplication(transpose, matrix, k, n, k);
    if (m_t_mult == NULL) {
        free_continuous_matrix(transpose);
        return -1.0;
    }
    trace = matrix_trace(m_t_mult, k);
    free_continuous_matrix(transpose);
    free_continuous_matrix(m_t_mult);
    return trace;
}


void converge_H_memory_freer(double **prev_H, double **cur_H, double **distance_matrix){
    /* Frees up H matrices for convenience. Matrices can be NULL as free_continuous_matrix which is used here handles it. 
    Input: 
        - double prev_H[][]: Matrix 1 we are freeing.
        - double prev_H[][]: Matrix 2 we are freeing.
        - double distance_matrix: Matrix 3 we are freeing.
    */
    free_continuous_matrix(prev_H);
    free_continuous_matrix(cur_H);
    free_continuous_matrix(distance_matrix);
}

double **converge_H(double **initial_H, double **W, int n, int k) {
    /* Continuously updates H until either convergence or until reaching max iterations, as per project instructions. Returns NULL on error.
    Input: 
        - double Initial_H[][]: Initial H matrix we received from Python.
        - double W[][]: Norm matrix.
        - int n: Size of norm matrix, number of rows in H.
        - int k: Number of columns in H.
    Returns:
        Final iteration of H. 
    */
    double **prev_H = NULL, **cur_H = NULL, **distance_matrix = NULL;
    double frobenius_distance_squared;
    int iteration;
    prev_H = matrix_deep_copy(initial_H, n, k);
    if (prev_H == NULL) {return NULL;}
    for (iteration = 0; iteration < max_iter; iteration++) {
        cur_H = update_H(prev_H, W, n, k);
        if (cur_H == NULL) {
            converge_H_memory_freer(prev_H, cur_H, distance_matrix);
            return NULL;
        }
        distance_matrix = matrix_subtraction(cur_H, prev_H, n, k);
        if (distance_matrix == NULL) {
            converge_H_memory_freer(prev_H, cur_H, distance_matrix);
            return NULL;
        }
        frobenius_distance_squared = frobenius_norm_squared(distance_matrix, n, k);
        free_continuous_matrix(prev_H);
        free_continuous_matrix(distance_matrix);
        if (frobenius_distance_squared == -1) {
            free_continuous_matrix(cur_H);
            return NULL;
        }
        if (frobenius_distance_squared < epsilon) {
            break;
        }
        prev_H = cur_H;
    }
    return cur_H;
}

typedef struct datapoints_wrapper {
    double **datapoints;
    int num_points;
    int dimension;
} datapoints_wrapper;

void datapoints_on_error_handler(datapoints_wrapper *datapoints){
    /* Function to handle deallocating datapoints memory in case of error. Includes print message.
    Input: 
        - datapoints_wrapper *datapoints: datapoints wrapper.
    */
    printf("An Error Has Occurred\n");
    free_matrix(datapoints->datapoints, datapoints->num_points);
    free(datapoints);
}

void invalid_file_read_error_handler(FILE *file, datapoints_wrapper *wrapper, int row) {
    /* Function to handle deallocating datapoints memory and file memory in case of error. Fully handles error's by deallocating memory and exiting.
    Input: 
        - FILE *file: pointer to file we need to close.
        - datapoints_wrapper *datapoints: datapoints wrapper.
        - int row: current row file failed on so datapoints_on_error_handler can properly deallocate datapoints matrix.
    */
    fclose(file);
    wrapper->num_points = row;
    datapoints_on_error_handler(wrapper); 
    exit(EXIT_FAILURE);
}

void invalid_file_read_error_handler_simple_case(FILE *file, datapoints_wrapper *wrapper) {
    /* Function to handle deallocating datapoints memory and file memory in case of error. Fully handles error's by deallocating memory and exiting.
    Input: 
        - FILE *file: pointer to file we need to close.
        - datapoints_wrapper *datapoints: datapoints wrapper.
    */
    fclose(file);
    datapoints_on_error_handler(wrapper); 
    exit(EXIT_FAILURE);
}

datapoints_wrapper* initialize_data(const char *filename) {
    /* Initializes a datapoints wrapper in order to read datapoints from file. Fully handles errors by deallocating memory and exiting program.
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
    if (!file || rows == 0) {
        printf("An Error Has Occurred\n");
        exit(EXIT_FAILURE);
    }

    wrapper = malloc(sizeof(datapoints_wrapper));
    wrapper->datapoints = malloc(rows * sizeof(double*));
    if (wrapper->datapoints == NULL) {
        printf("An Error Has Occurred\n");
        exit(EXIT_FAILURE);
    }
    wrapper->num_points = rows;
    fclose(file);
    return wrapper;
}

void determine_data_dimension(FILE *file, datapoints_wrapper *wrapper) {
    /* Reads the first data line to determine point dimension and stores it in the wrapper. Fully handles error's by deallocating memory and exiting.
    Input:
        - FILE *file: Pointer to the opened input file.
        - datapoints_wrapper *wrapper: Pointer to the wrapper to store dimension in.
    */
    char line[1024];
    char *line_copy = NULL;
    char *token;
    int cols = 0;
    size_t line_len;
    while (fgets(line, sizeof(line), file) != NULL) {
        if (line[0] != '\n') {
            break; 
        }
    }
    if (line[0] == '\n' || feof(file) || ferror(file)) {
        invalid_file_read_error_handler_simple_case(file, wrapper);
    }
    line_len = strlen(line) + 1; 
    line_copy = malloc(line_len); 
    if (line_copy == NULL) { 
        invalid_file_read_error_handler_simple_case(file, wrapper);
    }
    strcpy(line_copy, line); 
    if (line_copy == NULL) {
        invalid_file_read_error_handler_simple_case(file, wrapper);
    }
    token = strtok(line_copy, ",");
    while (token != NULL) {
        cols++;
        token = strtok(NULL, ",");
    }
    free(line_copy); 
    if (cols <= 0) { 
        invalid_file_read_error_handler_simple_case(file, wrapper);
    }
    wrapper->dimension = cols;
}

void read_and_store_datapoints(FILE *file, datapoints_wrapper *wrapper) {
    /* Reads data points from file, allocates rows, and stores values. Exits on error.
    Input:
        - FILE *file: Pointer to the opened input file (positioned correctly).
        - datapoints_wrapper *wrapper: Pointer to the wrapper to populate.
    */
    char line[1024];
    char *token;
    int row = 0;
    int col;
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '\n') continue;
        if (row >= wrapper->num_points) {
            invalid_file_read_error_handler(file, wrapper, row);
        }
        wrapper->datapoints[row] = malloc(wrapper->dimension * sizeof(double));
        if (wrapper->datapoints[row] == NULL) {
            invalid_file_read_error_handler(file, wrapper, row);
        }
        token = strtok(line, ","); 
        for (col = 0; col < wrapper->dimension; ++col) {
            if (token == NULL) { 
                invalid_file_read_error_handler(file, wrapper, row + 1); /* Row + 1 to deallocate the additional row that was allocated whose line reading failed*/
            }
            wrapper->datapoints[row][col] = atof(token);
            token = strtok(NULL, ",");
        }
        row++;
    }
    if (row != wrapper->num_points) {
        invalid_file_read_error_handler(file, wrapper, row);
    }
}

void populate_data(datapoints_wrapper *wrapper, const char *filename) {
    /* Main function to populate datapoints wrapper from file. Orchestrates dimension determination and data reading via helpers. 
    Fully handles error's by deallocating memory and exiting.   
    Input: 
        - datapoints_wrapper *wrapper: initialized datapoints wrapper.
        - char[] filename: string filepath to .txt file containing datapoints
    */
    FILE *file = fopen(filename, "r");

    if (!file) { 
        datapoints_on_error_handler(wrapper); 
        exit(EXIT_FAILURE); 
    }    

    determine_data_dimension(file, wrapper); /* Sets wrapper->dimension */

    rewind(file); /* Reset file pointer for full read */

    read_and_store_datapoints(file, wrapper); /* Allocates rows and reads data */

    fclose(file); 
    return; 
}

void sym(datapoints_wrapper *datapoints) {
    /* Wrapper function to calculate similarity matrix as per project instructions. Fully handles errors by deallocating memory and exiting.
    Input: 
        - datapoints_wrapper *datapoints: datapoints wrapper.
    */
    double **sym_matrix;
    int n = datapoints->num_points;
    int d = datapoints->dimension;
    sym_matrix = similarity_matrix(datapoints->datapoints, n, d);
    if (sym_matrix == NULL) {
        datapoints_on_error_handler(datapoints);
        exit(EXIT_FAILURE);
    }
    print_matrix(sym_matrix, n, n);
    free_continuous_matrix(sym_matrix);
}

void ddg(datapoints_wrapper *datapoints) {
    /* Wrapper function to calculate diagonal matrix as per project instructions. Fully handles errors by deallocating memory and exiting.
    Input: 
        - datapoints_wrapper *datapoints: datapoints wrapper.
    */
    double **sym_matrix;
    double **diag_matrix;
    int n = datapoints->num_points;
    int d = datapoints->dimension;
    sym_matrix = similarity_matrix(datapoints->datapoints, n, d);
    if (sym_matrix == NULL) {
        datapoints_on_error_handler(datapoints);
        exit(EXIT_FAILURE);
    }
    diag_matrix = diagonal_matrix(sym_matrix, n);
    if (diag_matrix == NULL) {
        free_continuous_matrix(sym_matrix);
        datapoints_on_error_handler(datapoints);
        exit(EXIT_FAILURE);
    }
    print_matrix(diag_matrix, n, n);
    free_continuous_matrix(sym_matrix);
    free_continuous_matrix(diag_matrix);
}

void norm(datapoints_wrapper *datapoints) {
    /* Wrapper function to calculate norm matrix as per project instructions. Fully handles errors by deallocating memory and exiting.
    Input: 
        - datapoints_wrapper *datapoints: datapoints wrapper.
    */
    double **sym_matrix;
    double **diag_matrix;
    double **normal_matrix;
    int n = datapoints->num_points;
    int d = datapoints->dimension;
    sym_matrix = similarity_matrix(datapoints->datapoints, n, d);
    if (sym_matrix == NULL) {
        datapoints_on_error_handler(datapoints);
        exit(EXIT_FAILURE);
    }
    diag_matrix = diagonal_matrix(sym_matrix, n);
    if (diag_matrix == NULL) {
        free_continuous_matrix(sym_matrix);
        datapoints_on_error_handler(datapoints);
        exit(EXIT_FAILURE);
    }
    normal_matrix = norm_matrix(sym_matrix, diag_matrix, n);
    if (normal_matrix == NULL) {
        free_continuous_matrix(sym_matrix);
        free_continuous_matrix(diag_matrix);
        datapoints_on_error_handler(datapoints);
        exit(EXIT_FAILURE);
    }
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
        printf("An Error Has Occurred\n");
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
        datapoints_on_error_handler(datapoints);
        exit(EXIT_FAILURE);
    }
    free_matrix(datapoints->datapoints, datapoints->num_points);
    free(datapoints);
    exit(EXIT_SUCCESS);
}
