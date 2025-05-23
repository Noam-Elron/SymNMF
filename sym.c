#include <math.h>
#include <stdlib.h>
#include "utils.h"

double euclidean_distance_squared(double point[], double other_point[], int point_dimension) {
    /* Calculates squared euclidean distance between two datapoints
    Input:
        - double point[]: First point we are comparing euclidean distance with.
        - double other_point[]: Second point we are comparing euclidean distance with.
        - int point_dimension: Number of coordinates in each point.
    Returns:
        Squared euclidean distance between the two points.
    */
    int i;
    double total = 0.0;
    for (i = 0; i < point_dimension; i++) {
        total += (point[i] - other_point[i]) * (point[i] - other_point[i]);
    }
    return total;
}


double **similarity_matrix(double **datapoints, int num_points, int point_dimension) {
    /* Creates similarity matrix as per project instructions. Returns NULL on error.
    Input: 
        - double Datapoints[][]: 2D Array, each element in it is a point who is itself an array of coordinates.datapoints
        - int num_points: Number of points in Datapoints, this is also the size of the similarity matrix as each entry in it corresponds to a point
        - int point_dimension: Number of coordinates in each point, we pass it in as euclidean_distance_squared has need of it in order to calculate distance between two points.
    Returns:
        2D Similarity Matrix
    */

    int i;
    int j;
    double **sym_matrix;
    double distance_squared;

    sym_matrix = continuous_matrix_creation(num_points, num_points);
    if (sym_matrix == NULL) {
        return NULL;
    }
    for (i = 0; i < num_points; i++) {
        for (j = 0; j < num_points; j++) {
            /* Calloc instantiates all elements to zero, therefore no need to set a_ii = 0 manually */
            if (i != j) {
                distance_squared = euclidean_distance_squared(datapoints[i], datapoints[j], point_dimension);
                sym_matrix[i][j] = exp(-(distance_squared / 2.0));
            }

        }
    }
    return sym_matrix;
}


