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


c_matrix_wrapper *py_matrix_to_c_matrix(PyObject *matrix_py_ptr) {
    /* Converts python matrix to c matrix
    Input: 
        - PyObject *matrix_py_ptr: Python matrix we want to convert to c matrix
    Returns:
        Wrapper for c matrix, includes matrix itself as well as matrix dimensions.
    */
    PyObject *temp_row_py_ptr;
    PyObject *coord_py_ptr;
    Py_ssize_t i;
    Py_ssize_t j;
    c_matrix_wrapper *wrapper;

    wrapper = malloc(sizeof(c_matrix_wrapper));
    wrapper->rows = PyList_Size(matrix_py_ptr);
    wrapper->matrix = (double**)calloc(wrapper->rows, sizeof(double*));

    for (i = 0; i < wrapper->rows; i++) {
        temp_row_py_ptr = PyList_GetItem(matrix_py_ptr, i);
        if (!PyList_Check(temp_row_py_ptr)) {
            printf("An Error has Occurred\n");
            exit(EXIT_FAILURE);
        }
        wrapper->cols = PyList_Size(temp_row_py_ptr);
        wrapper->matrix[i] = (double*)calloc(wrapper->cols, sizeof(double));
        for (j = 0; j < wrapper->cols; j++) {
            coord_py_ptr = PyList_GetItem(temp_row_py_ptr, j);
            if (Py_IS_TYPE(coord_py_ptr, &PyFloat_Type) == 0) {
                printf("An Error has Occurred\n");
                exit(EXIT_FAILURE);
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
    PyObject *matrix_py;
    PyObject *temp_matrix_row_py;
    PyObject *matrix_entry_py;
    Py_ssize_t i;
    Py_ssize_t j;
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



static PyObject* sym_c_wrapper(PyObject *self, PyObject *args) {
    /* Python-C Extension wrapper for calculating similarity matrix in C and returning it to Python program.
    Input: 
        - PyObject *self: reference to wrapper.
        - PyObject *args: Python arguments calling c function. 
        
    Returns:
        Python similarity matrix
    */
    PyObject *datapoints_matrix_py_ptr;
    c_matrix_wrapper *datapoints_wrapper;
    double **sim_matrix;
    PyObject* sym_matrix_py_ptr;
    
    if (!PyArg_ParseTuple(args, "O", &datapoints_matrix_py_ptr)) {
        printf("An Error has Occurred\n");
        exit(EXIT_FAILURE);
    }

    if (!PyList_Check(datapoints_matrix_py_ptr)) {
        printf("An Error has Occurred\n");
        exit(EXIT_FAILURE);
    }
    datapoints_wrapper = py_matrix_to_c_matrix(datapoints_matrix_py_ptr);
    sim_matrix = similarity_matrix(datapoints_wrapper->matrix, datapoints_wrapper->rows, datapoints_wrapper->cols);
    sym_matrix_py_ptr = c_matrix_to_py_matrix(sim_matrix, datapoints_wrapper->rows, datapoints_wrapper->cols);
    free_matrix(datapoints_wrapper->matrix, datapoints_wrapper->rows);
    free(datapoints_wrapper);
    free_continuous_matrix(sim_matrix);
    return sym_matrix_py_ptr;;
}

static PyObject* diag_c_wrapper(PyObject *self, PyObject *args) {
    /* Python-C Extension wrapper for calculating diagonal matrix in C and returning it to Python program.
    Input: 
        - PyObject *self: reference to wrapper.
        - PyObject *args: Python arguments calling c function. 
        
    Returns:
        Python diagonal matrix
    */
    PyObject *datapoints_matrix_py_ptr;
    c_matrix_wrapper *datapoints_wrapper;
    double **sim_matrix;
    double **diag_matrix;
    PyObject *diag_matrix_py_ptr;

    if (!PyArg_ParseTuple(args, "O", &datapoints_matrix_py_ptr)) {
        printf("An Error has Occurred\n");
        exit(EXIT_FAILURE);
    }

    if (!PyList_Check(datapoints_matrix_py_ptr)) {
        printf("An Error has Occurred\n");
        exit(EXIT_FAILURE);
    }

    datapoints_wrapper = py_matrix_to_c_matrix(datapoints_matrix_py_ptr);
    sim_matrix = similarity_matrix(datapoints_wrapper->matrix, datapoints_wrapper->rows, datapoints_wrapper->cols);
    diag_matrix = diagonal_matrix(sim_matrix, datapoints_wrapper->rows, datapoints_wrapper->rows);
    diag_matrix_py_ptr = c_matrix_to_py_matrix(diag_matrix, datapoints_wrapper->rows, datapoints_wrapper->rows);
    free_matrix(datapoints_wrapper->matrix, datapoints_wrapper->rows);
    free(datapoints_wrapper);
    free_continuous_matrix(sim_matrix);
    free_continuous_matrix(diag_matrix);
    return diag_matrix_py_ptr;;
}

static PyObject* norm_c_wrapper(PyObject *self, PyObject *args) {
    /* Python-C Extension wrapper for calculating norm matrix in C and returning it to Python program.
    Input: 
        - PyObject *self: reference to wrapper.
        - PyObject *args: Python arguments calling c function. 
        
    Returns:
        Python norm matrix
    */
    PyObject *datapoints_matrix_py_ptr;
    c_matrix_wrapper *datapoints_wrapper;
    double **sim_matrix;
    double **diag_matrix;
    double **nm_matrix;
    PyObject* norm_matrix_py_ptr;

    if (!PyArg_ParseTuple(args, "O", &datapoints_matrix_py_ptr)) {
        printf("An Error has Occurred\n");
        exit(EXIT_FAILURE);
    }
    if (!PyList_Check(datapoints_matrix_py_ptr)) {
        printf("An Error has Occurred\n");
        exit(EXIT_FAILURE);
    }
    datapoints_wrapper = py_matrix_to_c_matrix(datapoints_matrix_py_ptr);
    sim_matrix = similarity_matrix(datapoints_wrapper->matrix, datapoints_wrapper->rows, datapoints_wrapper->cols);
    diag_matrix = diagonal_matrix(sim_matrix, datapoints_wrapper->rows, datapoints_wrapper->rows);
    nm_matrix = norm_matrix(sim_matrix, diag_matrix, datapoints_wrapper->rows);
    norm_matrix_py_ptr = c_matrix_to_py_matrix(nm_matrix, datapoints_wrapper->rows, datapoints_wrapper->rows);
    free_matrix(datapoints_wrapper->matrix, datapoints_wrapper->rows);
    free(datapoints_wrapper);
    free_continuous_matrix(sim_matrix);
    free_continuous_matrix(diag_matrix);
    free_continuous_matrix(nm_matrix);
    return norm_matrix_py_ptr;;
}




static PyObject* symnmf_c_wrapper(PyObject *self, PyObject *args) {
    /* Python-C Extension wrapper for calculating SymNMF matrix in C and returning it to Python program.
    Input: 
        - PyObject *self: reference to wrapper.
        - PyObject *args: Python arguments calling c function. 
        
    Returns:
        Python symnmf matrix
    */
    double **symnmf_matrix;
    c_matrix_wrapper *initial_H_wrapper;
    c_matrix_wrapper *norm_wrapper;
    PyObject *initial_H_py_ptr;
    PyObject *norm_matrix_py_ptr;
    PyObject *symnmf_matrix_py_ptr;
    
    if (!PyArg_ParseTuple(args, "OO", &initial_H_py_ptr, &norm_matrix_py_ptr)) {
        printf("An Error has Occurred\n");
        exit(EXIT_FAILURE);
    }

    if (!PyList_Check(initial_H_py_ptr) || !PyList_Check(norm_matrix_py_ptr)) {
        printf("An Error has Occurred\n");
        exit(EXIT_FAILURE);
    }
    norm_wrapper = py_matrix_to_c_matrix(norm_matrix_py_ptr);
    initial_H_wrapper = py_matrix_to_c_matrix(initial_H_py_ptr);
    symnmf_matrix = converge_H(initial_H_wrapper->matrix, norm_wrapper->matrix, initial_H_wrapper->rows, initial_H_wrapper->cols);
    symnmf_matrix_py_ptr = c_matrix_to_py_matrix(symnmf_matrix, initial_H_wrapper->rows, initial_H_wrapper->cols);
    free_matrix(initial_H_wrapper->matrix, initial_H_wrapper->rows);
    free_matrix(norm_wrapper->matrix, norm_wrapper->rows);
    free(initial_H_wrapper);
    free(norm_wrapper);
    return symnmf_matrix_py_ptr;;
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