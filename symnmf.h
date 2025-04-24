
double **update_H(double **prev_H, double **W, int n, int k);

double frobenius_norm(double **matrix, int m, int n);

double **converge_H(double **initial_H, double **W, int n, int k);

typedef struct datapoints_wrapper;

datapoints_wrapper* initialize_data(const char *filename);

void populate_data(datapoints_wrapper *wrapper, const char *filename);

void sym(datapoints_wrapper *datapoints);

void ddg(datapoints_wrapper *datapoints);

void norm(datapoints_wrapper *datapoints);