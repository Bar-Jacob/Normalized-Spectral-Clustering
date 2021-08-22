#include "spkmeans.h"

static void fit_capi(PyObject *self, PyObject *args);

/*
when calling fit() from python, this function is called. 
getting arguments from python and pass it to the adequate goal function
*/
static void fit_capi(PyObject *self, PyObject *args)
{
    int k;
    char goal;
    int dimension_p;
    int num_of_points_p;
    double** data_points;
    PyObject *data_points_p;
    printf("in fit\n");

    printf("before first if");
    if (!(PyArg_ParseTuple(args, "isiiO", &k, &goal, &dimension_p, &num_of_points_p, &data_points_p)))
    {
        printf("before error in first if");
        printf("An Error Has Occured");
        exit(0);
    }

    printf("before second if");
    if (!PyList_Check(data_points_p))
    {
        printf("before error in second if");
        printf("An Error Has Occured");
        exit(0);
    }
    printf("before convert");
    data_points = convert_python_to_c(data_points_p, dimension_p, num_of_points_p);
    
    printf("before goals");
    switch (goal)
    {
    case 'w':
        wam_goal(data_points, num_of_points_p, dimension_p);
        break;
    case 'd':
        ddg_goal(data_points, num_of_points_p, dimension_p);
        break;
    case 'l':
        lnorm_goal(data_points, num_of_points_p, dimension_p);
        break;
    case 'j':
        printf("in jacobi goal");
        jacobi_goal(data_points, num_of_points_p, dimension_p);
        break;
    case 's':
        break;
    }

    printf("at the end of fit");
}
/*
building spkmeans module...
*/
static PyMethodDef spkmeansMethods[] = {
    {"fit",
     (PyCFunction) fit_capi,
     METH_VARARGS,
     PyDoc_STR("spkmeans algorithem")},
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