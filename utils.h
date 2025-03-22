double **continuous_matrix_creation(int m, int n);

double **matrix_deep_copy(double **matrix_to_copy, int m, int n);

double **matrix_multiplication(double **matrix, double **other_matrix, int m, int s, int n);

void free_matrix(double **matrix, int num_rows);

void free_continuous_matrix(double **continuous_matrix);

void print_matrix(double **matrix, int m, int n);
