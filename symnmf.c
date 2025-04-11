#include <math.h>
#include "utils.h"
#define beta 0.5
#define epsilon 0.0001 
#define max_iter 300


double **update_H(double **prev_H, double **W, int n, int k) {
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
        for (j = 0; i < k; j++) {
            next_H[i][j] = prev_H[i][j] * (1 - beta + beta*(w_h_mult[i][j]/h_h_t_h_mult[i][j]));
        }
    }

    free_matrix(w_h_mult, n);
    free_matrix(h_t, k);
    free_matrix(h_h_t_mult, n);
    free_matrix(h_h_t_h_mult, n);
    return next_H;
}   

double frobenius_norm(double **matrix, int m, int n) {
    double **transpose;
    double **m_t_mult;
    double trace;
    transpose = matrix_transpose(matrix, m, n);
    m_t_mult = matrix_multiplication(transpose, matrix, m, n, m);
    trace = matrix_trace(m_t_mult, m);
    return sqrt(trace);
}


double **converge_H(double **initial_H, double **W, int n, int k){
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
        free_matrix(prev_H, n);
        if (frobenius_distance < epsilon) {
            break;
        }
        prev_H = cur_H;
    }

    return cur_H;
}