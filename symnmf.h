
typedef struct datapoints_wrapper datapoints_wrapper;

void free_update_H_matrices(double **w_h_mult, double **h_t, double **h_h_t_mult, double **h_h_t_h_mult);

double **update_H(double **prev_H, double **W, int n, int k);

double frobenius_norm_squared(double **matrix, int m, int n);

void converge_H_memory_freer(double **prev_H, double **distance_matrix);

double **converge_H(double **initial_H, double **W, int n, int k);

void datapoints_on_error_handler(datapoints_wrapper *datapoints);

void invalid_file_read_error_handler(FILE *file, datapoints_wrapper *wrapper, int row);

void invalid_file_read_error_handler_simple_case(FILE *file, datapoints_wrapper *wrapper);

datapoints_wrapper* initialize_data(const char *filename);

void determine_data_dimension(FILE *file, datapoints_wrapper *wrapper);

void read_and_store_datapoints(FILE *file, datapoints_wrapper *wrapper);

void populate_data(datapoints_wrapper *wrapper, const char *filename);

void sym(datapoints_wrapper *datapoints);

void ddg(datapoints_wrapper *datapoints);

void norm(datapoints_wrapper *datapoints);