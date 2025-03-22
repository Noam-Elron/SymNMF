#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "sym.h"
#include "diagonal.h"
#include "norm.h"
#include "symnmf.h"

typedef struct datapoints_wrapper {
    double **datapoints;
    Py_ssize_t num_points;
    Py_ssize_t point_dimension;
} datapoints_wrapper;


datapoints_wrapper *datapoints_to_c_array(PyObject *datapoints_py_ptr) {
    PyObject *temp_datapoint_py_ptr;
    PyObject *coord_py_ptr;
    Py_ssize_t i;
    Py_ssize_t j;
    datapoints_wrapper *wrapper;

    wrapper = malloc(sizeof(datapoints_wrapper));
    wrapper->num_points = PyList_Size(datapoints_py_ptr);
    wrapper->datapoints = (double**)calloc(wrapper->num_points, sizeof(double*));

    for (i = 0; i < wrapper->num_points; i++) {
        temp_datapoint_py_ptr = PyList_GetItem(datapoints_py_ptr, i);
        if (!PyList_Check(temp_datapoint_py_ptr)) {
            printf("An Error has Occurred\n");
            exit(EXIT_FAILURE);
        }
        wrapper->point_dimension = PyList_Size(datapoints_py_ptr);
        wrapper->datapoints[i] = (double*)calloc(wrapper->point_dimension, sizeof(double));
        for (j = 0; j < wrapper->point_dimension; j++) {
            coord_py_ptr = PyList_GetItem(temp_datapoint_py_ptr, j);
            if (Py_IS_TYPE(coord_py_ptr, &PyFloat_Type) == 0) {
                printf("An Error has Occurred\n");
                exit(EXIT_FAILURE);
            }
            wrapper->datapoints[i][j] = PyFloat_AsDouble(coord_py_ptr);
        }
    }
    return wrapper;
}


PyObject *c_matrix_to_py_matrix(double **matrix, Py_ssize_t m, Py_ssize_t n) {
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
    PyObject *datapoints_py_ptr;
    datapoints_wrapper *wrapper;
    double **sim_matrix;
    
    if (!PyArg_ParseTuple(args, "O", &datapoints_py_ptr)) {
        printf("An Error has Occurred\n");
        exit(EXIT_FAILURE);
    }

    if (!PyList_Check(datapoints_py_ptr)) {
        printf("An Error has Occurred\n");
        exit(EXIT_FAILURE);
    }
    wrapper = datapoints_to_c_array(datapoints_py_ptr);
    sim_matrix = similarity_matrix(wrapper->datapoints, wrapper->num_points, wrapper->point_dimension);
    return c_matrix_to_py_matrix(sim_matrix, wrapper->num_points, wrapper->point_dimension);
}

static PyObject* diag_c_wrapper(PyObject *self, PyObject *args) {
    PyObject *datapoints_py_ptr;
    datapoints_wrapper *wrapper;
    double **sim_matrix;
    double **diag_matrix;

    if (!PyArg_ParseTuple(args, "O", &datapoints_py_ptr)) {
        printf("An Error has Occurred\n");
        exit(EXIT_FAILURE);
    }

    if (!PyList_Check(datapoints_py_ptr)) {
        printf("An Error has Occurred\n");
        exit(EXIT_FAILURE);
    }

    wrapper = datapoints_to_c_array(datapoints_py_ptr);
    sim_matrix = similarity_matrix(wrapper->datapoints, wrapper->num_points, wrapper->point_dimension);
    diag_matrix = diagonal_matrix(sim_matrix, wrapper->num_points, wrapper->num_points);
    return c_matrix_to_py_matrix(diag_matrix, wrapper->num_points, wrapper->num_points);
}

static PyObject* norm_c_wrapper(PyObject *self, PyObject *args) {
    PyObject *datapoints_py_ptr;
    datapoints_wrapper *wrapper;
    double **sim_matrix;
    double **diag_matrix;
    double **nm_matrix;

    if (!PyArg_ParseTuple(args, "O", &datapoints_py_ptr)) {
        printf("An Error has Occurred\n");
        exit(EXIT_FAILURE);
    }

    if (!PyList_Check(datapoints_py_ptr)) {
        printf("An Error has Occurred\n");
        exit(EXIT_FAILURE);
    }

    wrapper = datapoints_to_c_array(datapoints_py_ptr);
    sim_matrix = similarity_matrix(wrapper->datapoints, wrapper->num_points, wrapper->point_dimension);
    diag_matrix = diagonal_matrix(sim_matrix, wrapper->num_points, wrapper->num_points);
    nm_matrix = norm_matrix(sim_matrix, diag_matrix, wrapper->num_points);
    return c_matrix_to_py_matrix(nm_matrix, wrapper->num_points, wrapper->num_points);
}

static PyObject* symnmf_c_wrapper(PyObject *self, PyObject *args) {
    PyObject *symnmf_matrix;
    return symnmf_matrix;
}

static PyMethodDef SymNMFModules[] = {
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
    "symnmf",
    "Python interface for the SymNMF C module",
    -1,
    SymNMFModules
  };
  
PyMODINIT_FUNC PyInit_symnmf(void) {
      PyObject *module;
      module = PyModule_Create(&SymNMFModule);
      if (!module) {
          return NULL;
      }
      return module;
}