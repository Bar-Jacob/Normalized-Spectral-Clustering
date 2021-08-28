#include "spkmeans.h"

static PyObject* fit_capi(PyObject* self, PyObject* args);
static PyObject* fit_capi_pp(PyObject* self, PyObject* args);


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
    PyObject *data_points_p;

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
        return Py_BuildValue("O", 
        spk_goal_python(data_points, num_of_points_p, dimension_p, k));
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
    kmeans_pp(k, dimension_p, num_of_points_p, centroids_locations, data_points_p);
    Py_RETURN_NONE;
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