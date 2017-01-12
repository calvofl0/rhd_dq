#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define __LTEBACKGR_H_NO_FUNCTION_PROTOTYPES__
#define PyNone Py_BuildValue("s", NULL)

#define PyArray_FSimpleNewFromData(nd, dims, typenum, data) \
        PyArray_New(&PyArray_Type, nd, dims, typenum, NULL, \
                    data, 0, NPY_ARRAY_FARRAY, NULL)

#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>

#ifdef NAN
#define NAND nan("")
#endif

#define max2(x,y) (((x)>(y))?(x):(y))

/* Docstrings */

static char module_docstring[] =
  "This module provides the 3D level function";
static char level_docstring[] =
  "Compute the l-level of a 3D scalar field (surface where the field is l). " \
	  "The fields needs to be decreasing in the vertical direction.";
  
/* Module function prototypes */

static PyObject* pylevel_level(PyObject *self, PyObject *args, \
		PyObject *kwargs);

static PyMethodDef module_methods[] = {
	{"level", (PyCFunction)pylevel_level, \
		0, level_docstring},
	{NULL, NULL, 0, NULL}
};

static char* level_kwlist[] = {
	"field", "level", "start", NULL
};

PyMODINIT_FUNC initpylevel(void) {
	/* Initialize python module */
	PyObject *m = Py_InitModule3("pylevel", module_methods, \
			module_docstring);
	if (m == NULL) {
		return;
	}

	/* Load `numpy` functionality */
	import_array();

}

static PyObject* pylevel_level(PyObject *self, PyObject *args, \
		PyObject *kwargs)
{
	/* Input arguments */
	float level;
	PyObject *field_obj = PyNone;
	PyObject *start_obj = PyNone;

	/* Local variables */
	PyObject *out = PyNone;
	npy_intp *dims_intp;
	PyArrayObject *field_arr = NULL, *start_arr = NULL;
	float *field = NULL;
	float *depth = NULL;
	float *start = NULL;
	int *start_int = NULL;
	float f0, f1;
	int nx, ny, nz;
	int i, j, k, k1;

	/* Parse input arguments */
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, \
				"Of|O", level_kwlist, \
				&field_obj, &level, &start_obj)) {
		return NULL;
	}

	/* Interprete input objects as numpy arrays */
	field_arr = (PyArrayObject*)PyArray_FROM_OTF(field_obj, \
			NPY_FLOAT, NPY_ARRAY_IN_FARRAY);
	if (start_obj) {
		start_arr = (PyArrayObject*)PyArray_FROM_OTF(start_obj, \
				NPY_FLOAT, NPY_ARRAY_IN_FARRAY);
	}

	/* Throw an exception if that didn't work */
	if (field_arr == NULL) {
		Py_XDECREF(field_arr);
		PyErr_SetString(PyExc_ValueError, \
				"pylevel_level: Failed parsing arguments.");
		return NULL;
	}
	if (start_obj && (start_arr == NULL)) {
		Py_XDECREF(field_arr);
		Py_XDECREF(start_arr);
		PyErr_SetString(PyExc_ValueError, \
				"pylevel_level: Failed parsing arguments.");
		return NULL;
	}

	if (PyArray_NDIM(field_arr) != 3) {
		Py_XDECREF(field_arr);
		PyErr_SetString(PyExc_ValueError, \
				"pylevel_level: Input scalar field must be 3-dimensional.");
		return NULL;
	}

	field = PyArray_DATA(field_arr);

	dims_intp = PyArray_DIMS(field_arr);
	nx = (int)dims_intp[0];
	ny = (int)dims_intp[1];
	nz = (int)dims_intp[2];
	depth = (float*)PyMem_Malloc(nx*ny*sizeof(float));
	if (depth == NULL) {
		Py_XDECREF(field_arr);
		PyErr_SetString(PyExc_MemoryError, \
				"pylevel_level: Could not allocated memory.");
		return NULL;
	}

	start_int = (int*)malloc(nx*ny*sizeof(int));
	if(start_arr) {
		start = PyArray_DATA(start_arr);
		for (i=0; i<nx*ny; i++) start_int[i] = max2(0,(int)ceilf(start[i]-1));
	}
	else for (i=0; i<nx*ny; i++) start_int[i] = 0;

	/* Computations */
	for (i=0; i<nx; i++) {
		for (j=0; j<ny; j++) {
			/* depth[i+dims_intp[0]*j] = field[i+dims_intp[0]*j]; */
			k1 = -1;
			for (k=start_int[i+nx*j]; k<nz; k++) {
				if (field[i+nx*j+nx*ny*k] <= level) {
					k1 = k;
					break;
				}
			}
			if (k1 < 1) depth[i+nx*j] = NAND;
			else {
				f0 = field[i+nx*j+nx*ny*(k1-1)];
				f1 = field[i+nx*j+nx*ny*k1];
				depth[i+nx*j] = (k1-1.) + (level-f0)/(f1-f0);
			}
		}
	}
	free(start_int);

	Py_XDECREF(field_arr);
	if (start_arr) Py_XDECREF(start_arr);
	out = PyArray_FSimpleNewFromData(2, dims_intp, NPY_FLOAT, depth);
	PyArray_ENABLEFLAGS((PyArrayObject*)out, NPY_ARRAY_OWNDATA);

	return out;
}
