#include "spkmeans.h"
#define PY_SSIZE_T_CLEAN
#include <Python.h>

static PyObject* fit_capi(PyObject* self, PyObject* args);
static PyObject* fit_capi_pp(PyObject* self, PyObject* args);
PyObject* cToPyObject(double** T, int dimension, int num_of_points, int k);
double** convert_python_to_c(PyObject* data_points_p, int dimension, int num_of_points);
int* convert__array_python_to_c(PyObject* array_p, int num_of_cells);


/*
when calling fit() from python, this function is called. 
getting arguments from python and pass it to the adequate goal function
*/
static PyObject* fit_capi(PyObject* self, PyObject* args)
{
    int k;
    int goal;
    int dimension_p;
    int num_of_points_p;
    double** data_points;
    PyObject* data_points_p;
    double** spk_result;

    if (!(PyArg_ParseTuple(args, "iiiiO", &k, &goal, &dimension_p, &num_of_points_p, &data_points_p)))
    {
        printf("An Error Has Occured");
        exit(0);
    }
    if (!PyList_Check(data_points_p))
    {
        printf("An Error Has Occured");
        exit(0);
    }
    data_points = convert_python_to_c(data_points_p, dimension_p, num_of_points_p);
    
    switch (goal)
    {
    case 1:
        wam_goal(data_points, num_of_points_p, dimension_p);
        exit(0);
    case 2:
        ddg_goal(data_points, num_of_points_p, dimension_p);
        exit(0);
    case 3:
        lnorm_goal(data_points, num_of_points_p, dimension_p);
        exit(0);
    case 4:
        jacobi_goal(data_points, num_of_points_p);
        exit(0);
    case 0:
        if (k == 0)
        {
            k = (int)spk_goal_python(data_points, num_of_points_p, dimension_p, k, 1)[0][0];
        }
        spk_result = spk_goal_python(data_points, num_of_points_p, dimension_p, k, 2);
        return Py_BuildValue("O",
                             cToPyObject(spk_result, k, num_of_points_p, k));
    }
    }
    Py_RETURN_NONE;
}

/*
when calling fit_pp() from python, this function is called. 
getting arguments from python and pass it to kmeans_pp
*/
static PyObject* fit_capi_pp(PyObject* self, PyObject* args)
{
    int k;
    int dimension_p;
    int num_of_points_p;
    PyObject* centroids_locations;
    PyObject* data_points_p;
    int* centroids_locations_c;
    double** data_points_c;

    if (!(PyArg_ParseTuple(args, "iiiOO", &k, &dimension_p, 
                    &num_of_points_p, &centroids_locations, &data_points_p)))
    {
        printf("An Error Has Occured");
        exit(0);
    }
    if (!PyList_Check(centroids_locations) || !PyList_Check(data_points_p))
    {
        printf("An Error Has Occured");
        exit(0);
    }
    data_points_c = convert_python_to_c(data_points_p, dimension_p, num_of_points_p);
    centroids_locations_c = convert__array_python_to_c(centroids_locations, k);
    kmeans_pp(k, dimension_p, num_of_points_p, centroids_locations_c, data_points_c);
    Py_RETURN_NONE;
}

/*
converting all data points from python into data points in C
*/
double** convert_python_to_c(PyObject* data_points_p, int dimension, int num_of_points)
{
    int cnt = 0;
    int i = 0;
    int j = 0;
    double** data_points = (double**)calloc(num_of_points, sizeof(*data_points));
    assert(data_points != NULL && "An Error Has Occured");

    for (i = 0; i < num_of_points; i++)
    {
        data_points[i] = (double *)calloc(dimension, sizeof(*data_points[i]));
        assert(data_points[i] != NULL && "An Error Has Occured");
        for (j = 0; j < dimension; j++)
        {
            data_points[i][j] = PyFloat_AsDouble(PyList_GetItem(data_points_p, cnt));
            cnt++;
        }
    }
    return data_points;
}

/*
converting a list from python into C array
*/
int* convert__array_python_to_c(PyObject* array_p, int num_of_cells)
{
    int i = 0;
    int* result = (int*)calloc(num_of_cells, sizeof(int));
    assert(result != NULL && "An Error Has Occured");

    for (i = 0; i < num_of_cells; i++)
    {
        result[i] = (int) PyLong_AsLong(PyList_GetItem(array_p, i));
    }
    return result;
}

/*
after finishing running kmeans algorithm we want to return the results to python 
converting types from C to python
*/
PyObject* cToPyObject(double** T, int dimension, int num_of_points, int k)
{
    PyObject* points_py;
    int i = 0;
    int j = 0;
    PyObject* value;

    /***adding another cell for k***/
    points_py = PyList_New(num_of_points+1);
    for (i = 0; i < num_of_points; i++)
    {
        PyObject* curr_vector;
        curr_vector = PyList_New(dimension);
        for (j = 0; j < dimension; j++)
        {
            value = Py_BuildValue("d", T[i][j]);
            PyList_SetItem(curr_vector, j, value);
        }
        /*
        adding the PyObject centroid to the PyList clusters
        */
        PyList_SetItem(points_py, i, curr_vector);
    }
    /***getting k as last var***/
    value = Py_BuildValue("i", k);
    PyList_SetItem(points_py, num_of_points, value);
    free_memory(T, num_of_points);
    return points_py;
}

/*
building spkmeans module...
*/
static PyMethodDef spkmeansMethods[] = {
    {"fit",
     (PyCFunction) fit_capi,
     METH_VARARGS,
     PyDoc_STR("spkmeans algorithem")},
    {"fit_pp",
     (PyCFunction) fit_capi_pp,
     METH_VARARGS,
     PyDoc_STR("spkmeans algorithem2")},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef =
    {
        PyModuleDef_HEAD_INIT,
        "spkmeansmodule",
        NULL,
        -1,
        spkmeansMethods};

PyMODINIT_FUNC
PyInit_spkmeansmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m)
    {
        return NULL;
    }
    return m;
}