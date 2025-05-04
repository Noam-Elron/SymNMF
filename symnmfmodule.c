#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "utils.h"
#include "sym.h"
#include "diagonal.h"
#include "norm.h"
#include "symnmf.h"

typedef struct c_matrix_wrapper {
    double **matrix;
    Py_ssize_t rows;
    Py_ssize_t cols;
} c_matrix_wrapper;


void py_matrix_to_c_matrix_error_handler(c_matrix_wrapper *wrapper, int cur_num_rows) {
    /* Function to handle deallocating memory in case of error. 
    Input: 
        - c_matrix_wrapper *wrapper: c matrix wrapper we are freeing.
        - int cur_num_rows: number of currently allocated rows in wrapper->matrix.
    */
    if (cur_num_rows > 0) {
        free_matrix(wrapper->matrix, cur_num_rows);
    }
    free(wrapper);
}


c_matrix_wrapper *py_matrix_to_c_matrix(PyObject *matrix_py_ptr) {
    /* Converts python matrix to c matrix. Returns NULL on error.
    Input: 
        - PyObject *matrix_py_ptr: Python matrix we want to convert to c matrix
    Returns:
        Wrapper for c matrix, includes matrix itself as well as matrix dimensions. */
    PyObject *temp_row_py_ptr, *coord_py_ptr;
    Py_ssize_t i, j;
    c_matrix_wrapper *wrapper;
    wrapper = malloc(sizeof(c_matrix_wrapper));
    if (wrapper == NULL) {return NULL;}
    wrapper->rows = PyList_Size(matrix_py_ptr);
    wrapper->matrix = (double**)calloc(wrapper->rows, sizeof(double*));
    if (wrapper->matrix == NULL) {
        py_matrix_to_c_matrix_error_handler(wrapper, 0);
        return NULL;
    }
    for (i = 0; i < wrapper->rows; i++) {
        temp_row_py_ptr = PyList_GetItem(matrix_py_ptr, i);
        if (!PyList_Check(temp_row_py_ptr)) {
            py_matrix_to_c_matrix_error_handler(wrapper, i);
            return NULL;
        }
        wrapper->cols = PyList_Size(temp_row_py_ptr);
        wrapper->matrix[i] = (double*)calloc(wrapper->cols, sizeof(double));
        if (wrapper->matrix[i] == NULL) {
            py_matrix_to_c_matrix_error_handler(wrapper, i);
            return NULL;
        }
        for (j = 0; j < wrapper->cols; j++) {
            coord_py_ptr = PyList_GetItem(temp_row_py_ptr, j);
            if (Py_IS_TYPE(coord_py_ptr, &PyFloat_Type) == 0) {
                py_matrix_to_c_matrix_error_handler(wrapper, i);
                return NULL;
            }
            wrapper->matrix[i][j] = PyFloat_AsDouble(coord_py_ptr);
        }
    }
    return wrapper;
}


PyObject *c_matrix_to_py_matrix(double **matrix, Py_ssize_t m, Py_ssize_t n) {
    /* Converts c matrix to python matrix
    Input: 
        - double **matrix: C matrix we want to convert to Python matrix
        - Py_ssize_t m: Number of rows in matrix
        - Py_ssize_t n: Number of columns in matrix
    Returns:
        Created python matrix equivalent of given C matrix.
    */
    PyObject *matrix_py, *temp_matrix_row_py, *matrix_entry_py;
    Py_ssize_t i, j;
    matrix_py = PyList_New(m);
    for (i = 0; i < m; i++) {
        temp_matrix_row_py = PyList_New(n);
        PyList_SetItem(matrix_py, i, temp_matrix_row_py);
        for (j = 0; j < n; j++) {
            matrix_entry_py = PyFloat_FromDouble(matrix[i][j]);
            PyList_SetItem(temp_matrix_row_py, j, matrix_entry_py);
        }
    }
    return matrix_py;
}


void wrapper_function_memory_deallocator(double **sim_matrix, double **diag_matrix, double **normal_matrix, double **symnmf_matrix, c_matrix_wrapper *matrix_wrapper_1, c_matrix_wrapper *matrix_wrapper_2) {
    /* Frees up main project matrices and datapoints wrapper for convenience. 
    Input: 
        Items can be NULL as free_continuous_matrix and free which are used here handle it. Includes error print.
        - double sim_matrix[][]: Matrix 1 we are freeing.
        - double diag_matrix[][]: Matrix 2 we are freeing.
        - double normal_matrix[][]: Matrix 3 we are freeing.
        - c_matrix_wrapper *matrix_wrapper_1: Wrapper 1 for matrix we are freeing.
        - c_matrix_wrapper *matrix_wrapper_2: Wrapper 2 for matrix we are freeing.
    */
    free_continuous_matrix(sim_matrix);
    free_continuous_matrix(diag_matrix);
    free_continuous_matrix(normal_matrix);
    free_continuous_matrix(symnmf_matrix);
    if (matrix_wrapper_1 != NULL) {
        free_matrix(matrix_wrapper_1->matrix, matrix_wrapper_1->rows);
        free(matrix_wrapper_1);
    }
    if (matrix_wrapper_2 != NULL) {
    free_matrix(matrix_wrapper_2->matrix, matrix_wrapper_2->rows);
    free(matrix_wrapper_2);
    }
}


void wrapper_function_error_handler(double **sim_matrix, double **diag_matrix, double **normal_matrix, double **symnmf_matrix, c_matrix_wrapper *matrix_wrapper_1, c_matrix_wrapper *matrix_wrapper_2) {
    /* Handles errors in project's main wrapper functions by deallocating memory, printing standard error message and exiting program.
    Input: 
        Items can be NULL as free_continuous_matrix and free which are used here handle it. Includes error print.
        - double sim_matrix[][]: Matrix 1 we are freeing.
        - double diag_matrix[][]: Matrix 2 we are freeing.
        - double normal_matrix[][]: Matrix 3 we are freeing.
        - c_matrix_wrapper *matrix_wrapper_1: Wrapper 1 for matrix we are freeing.
        - c_matrix_wrapper *matrix_wrapper_2: Wrapper 2 for matrix we are freeing.
    */
    printf("An Error Has Occurred\n");
    wrapper_function_memory_deallocator(sim_matrix, diag_matrix, normal_matrix, symnmf_matrix, matrix_wrapper_1, matrix_wrapper_2);
    exit(EXIT_FAILURE);
}


static PyObject* sym_c_wrapper(PyObject *self, PyObject *args) {
    /* Python-C Extension wrapper for calculating similarity matrix in C and returning it to Python program. Fully handles errors by deallocating memory and exiting program.
    Input: 
        - PyObject *self: reference to wrapper.
        - PyObject *args: Python arguments calling c function. 
    Returns:
        Python similarity matrix
    */
    PyObject *datapoints_matrix_py_ptr, *sym_matrix_py_ptr;
    c_matrix_wrapper *datapoints_wrapper;
    double **sim_matrix;
    
    if (!PyArg_ParseTuple(args, "O", &datapoints_matrix_py_ptr)) {
        printf("An Error Has Occurred\n");
        exit(EXIT_FAILURE);
    }
    if (!PyList_Check(datapoints_matrix_py_ptr)) {
        printf("An Error Has Occurred\n");
        exit(EXIT_FAILURE);
    }
    datapoints_wrapper = py_matrix_to_c_matrix(datapoints_matrix_py_ptr);
    if (datapoints_wrapper == NULL) {
        wrapper_function_error_handler(NULL, NULL, NULL, NULL, NULL, NULL);
    }
    sim_matrix = similarity_matrix(datapoints_wrapper->matrix, datapoints_wrapper->rows, datapoints_wrapper->cols);
    if (sim_matrix == NULL) {
        wrapper_function_error_handler(NULL, NULL, NULL, NULL, datapoints_wrapper, NULL);
    }
    sym_matrix_py_ptr = c_matrix_to_py_matrix(sim_matrix, datapoints_wrapper->rows, datapoints_wrapper->rows);
    wrapper_function_memory_deallocator(sim_matrix, NULL, NULL, NULL, datapoints_wrapper, NULL);
    return sym_matrix_py_ptr;
}


static PyObject* diag_c_wrapper(PyObject *self, PyObject *args) {
    /* Python-C Extension wrapper for calculating diagonal matrix in C and returning it to Python program. Fully handles error by deallocating memory and exiting program.
    Input: 
        - PyObject *self: reference to wrapper.
        - PyObject *args: Python arguments calling c function.     
    Returns:
        Python diagonal matrix
    */
    PyObject *datapoints_matrix_py_ptr, *diag_matrix_py_ptr;
    c_matrix_wrapper *datapoints_wrapper;
    double **sim_matrix, **diag_matrix;

    if (!PyArg_ParseTuple(args, "O", &datapoints_matrix_py_ptr)) {
        printf("An Error Has Occurred\n");
        exit(EXIT_FAILURE);
    }
    if (!PyList_Check(datapoints_matrix_py_ptr)) {
        printf("An Error Has Occurred\n");
        exit(EXIT_FAILURE);
    }

    datapoints_wrapper = py_matrix_to_c_matrix(datapoints_matrix_py_ptr);
    if (datapoints_wrapper == NULL) {
        wrapper_function_error_handler(NULL, NULL, NULL, NULL, NULL, NULL);
    }
    sim_matrix = similarity_matrix(datapoints_wrapper->matrix, datapoints_wrapper->rows, datapoints_wrapper->cols);
    if (sim_matrix == NULL) {
        wrapper_function_error_handler(NULL, NULL, NULL, NULL, datapoints_wrapper, NULL);
    }
    diag_matrix = diagonal_matrix(sim_matrix, datapoints_wrapper->rows);
    if (diag_matrix == NULL) {
        wrapper_function_error_handler(sim_matrix, NULL, NULL, NULL, datapoints_wrapper, NULL);
    }
    diag_matrix_py_ptr = c_matrix_to_py_matrix(diag_matrix, datapoints_wrapper->rows, datapoints_wrapper->rows);
    wrapper_function_memory_deallocator(sim_matrix, diag_matrix, NULL, NULL, datapoints_wrapper, NULL);
    return diag_matrix_py_ptr;;
}


static PyObject* norm_c_wrapper(PyObject *self, PyObject *args) {
    /* Python-C Extension wrapper for calculating norm matrix in C and returning it to Python program. Fully handles errors by deallocating memory and exiting program.
    Input: 
        - PyObject *self: reference to wrapper.
        - PyObject *args: Python arguments calling c function. 
    Returns:
        Python norm matrix
    */
    PyObject *datapoints_matrix_py_ptr, *norm_matrix_py_ptr;
    c_matrix_wrapper *datapoints_wrapper;
    double **sim_matrix, **diag_matrix, **nm_matrix;
    if (!PyArg_ParseTuple(args, "O", &datapoints_matrix_py_ptr)) {
        printf("An Error Has Occurred\n");
        exit(EXIT_FAILURE);
    }
    if (!PyList_Check(datapoints_matrix_py_ptr)) {
        printf("An Error Has Occurred\n");
        exit(EXIT_FAILURE);
    }
    datapoints_wrapper = py_matrix_to_c_matrix(datapoints_matrix_py_ptr);
    if (datapoints_wrapper == NULL) {
        wrapper_function_error_handler(NULL, NULL, NULL, NULL, NULL, NULL);
    }
    sim_matrix = similarity_matrix(datapoints_wrapper->matrix, datapoints_wrapper->rows, datapoints_wrapper->cols);
    if (sim_matrix == NULL) {
        wrapper_function_error_handler(NULL, NULL, NULL, NULL, datapoints_wrapper, NULL);
    }
    diag_matrix = diagonal_matrix(sim_matrix, datapoints_wrapper->rows);
    if (diag_matrix == NULL) {
        wrapper_function_error_handler(sim_matrix, NULL, NULL, NULL, datapoints_wrapper, NULL);
    }
    nm_matrix = norm_matrix(sim_matrix, diag_matrix, datapoints_wrapper->rows);
    if (nm_matrix == NULL) {
        wrapper_function_error_handler(sim_matrix, diag_matrix, NULL, NULL, datapoints_wrapper, NULL);
    }
    norm_matrix_py_ptr = c_matrix_to_py_matrix(nm_matrix, datapoints_wrapper->rows, datapoints_wrapper->rows);
    wrapper_function_memory_deallocator(sim_matrix, diag_matrix, nm_matrix, NULL, datapoints_wrapper, NULL);
    return norm_matrix_py_ptr;;
}


static PyObject* symnmf_c_wrapper(PyObject *self, PyObject *args) {
    /* Python-C Extension wrapper for calculating SymNMF matrix in C and returning it to Python program. Fully handles errors by deallocating memory and exiting program.
    Input: 
        - PyObject *self: reference to wrapper.
        - PyObject *args: Python arguments calling c function. 
    Returns:
        Python symnmf matrix
    */
    double **symnmf_matrix;
    c_matrix_wrapper *initial_H_wrapper, *norm_wrapper;
    PyObject *initial_H_py_ptr, *norm_matrix_py_ptr, *symnmf_matrix_py_ptr;
    if (!PyArg_ParseTuple(args, "OO", &initial_H_py_ptr, &norm_matrix_py_ptr)) {
        printf("An Error Has Occurred\n");
        exit(EXIT_FAILURE);
    }
    if (!PyList_Check(initial_H_py_ptr) || !PyList_Check(norm_matrix_py_ptr)) {
        printf("An Error Has Occurred\n");
        exit(EXIT_FAILURE);
    }
    norm_wrapper = py_matrix_to_c_matrix(norm_matrix_py_ptr);
    if (norm_wrapper == NULL) {
        wrapper_function_error_handler(NULL, NULL, NULL, NULL, NULL, NULL);
    }
    initial_H_wrapper = py_matrix_to_c_matrix(initial_H_py_ptr);
    if (norm_wrapper == NULL) {
        wrapper_function_error_handler(NULL, NULL, NULL, NULL, norm_wrapper, NULL);
    }
    symnmf_matrix = converge_H(initial_H_wrapper->matrix, norm_wrapper->matrix, initial_H_wrapper->rows, initial_H_wrapper->cols);
    if (symnmf_matrix == NULL) {
        wrapper_function_error_handler(NULL, NULL, NULL, NULL, norm_wrapper, initial_H_wrapper);
    }
    symnmf_matrix_py_ptr = c_matrix_to_py_matrix(symnmf_matrix, initial_H_wrapper->rows, initial_H_wrapper->cols);
    wrapper_function_memory_deallocator(NULL, NULL, NULL, symnmf_matrix, norm_wrapper, initial_H_wrapper);
    return symnmf_matrix_py_ptr;
}


static PyMethodDef SymNMFMethods[] = {
    {
        "sym", 
        (PyCFunction) sym_c_wrapper,
        METH_VARARGS,
        "sym C Wrapper"
    },
    {
        "diag", 
        (PyCFunction) diag_c_wrapper,
        METH_VARARGS,
        "diag C Wrapper"
    },
    {
        "norm", 
        (PyCFunction) norm_c_wrapper,
        METH_VARARGS,
        "norm C Wrapper"
    },
    {
        "symnmf", 
        (PyCFunction) symnmf_c_wrapper,
        METH_VARARGS,
        "SymNMF C Wrapper"
    },
    {NULL, NULL, 0, NULL}
  };

  
static struct PyModuleDef SymNMFModule = {
    PyModuleDef_HEAD_INIT,
    "symnmf_c",
    "Python interface for the SymNMF C module",
    -1,
    SymNMFMethods
  };
  
  
PyMODINIT_FUNC PyInit_symnmf_c(void) {
    /* Python-C interface module creation, created as shown in class.
    */
      PyObject *module;
      module = PyModule_Create(&SymNMFModule);
      if (!module) {
          return NULL;
      }
      return module;
}