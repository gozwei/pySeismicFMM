#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_10_API_VERSION
#include <numpy/arrayobject.h>
#include "3DFMM.h"



struct module_state {
    PyObject *error;
};

#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))

static PyObject* FMM3D(PyObject* self, PyObject* args)
{
	PyArrayObject *model_obj;
	PyArrayObject *time_obj;
	PyArrayObject *accepted_obj;
	PyArrayObject *lat_obj;
	PyArrayObject *lon_obj;
	PyArrayObject *h_obj;
	PyArrayObject *trace_obj;
	int N;

    if (!PyArg_ParseTuple(args, "OOOOOOOi", &model_obj, &time_obj, &accepted_obj, &lat_obj, &lon_obj, &h_obj, &trace_obj, &N))
	{
		Py_INCREF(Py_None);
		return Py_None;
	}

	float *MODEL	= static_cast<float *>(PyArray_DATA(model_obj));
	float *TIME		= static_cast<float *>(PyArray_DATA(time_obj));
	bool *ACCEPTED	= static_cast<bool *>(PyArray_DATA(accepted_obj));
	double *LAT		= static_cast<double *>(PyArray_DATA(lat_obj));
	double *LON		= static_cast<double *>(PyArray_DATA(lon_obj));
	double *H		= static_cast<double *>(PyArray_DATA(h_obj));
	int *TRACE		= static_cast<int *>(PyArray_DATA(trace_obj));

	_FMM3D(MODEL, TIME, ACCEPTED, LAT, LON, H, TRACE, N);

	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject* SetModelSize(PyObject* self, PyObject* args)
{
	int _size_lat, _size_lon, _size_z;
	double _dlat, _dlon, _dz;
	
	if (!PyArg_ParseTuple(args, "iiiddd", &_size_lat, &_size_lon, &_size_z, &_dlat, &_dlon, &_dz))
	{
		Py_INCREF(Py_None);
		return Py_None;
	}
	
	_SetModelSize(_size_lat, _size_lon, _size_z, _dlat, _dlon, _dz);
	
	Py_INCREF(Py_None);
	return Py_None;
	
}

static PyMethodDef FMMMethods[] = {
    {"FMM3D", FMM3D, METH_VARARGS, "..."},
	{"SetModelSize", SetModelSize, METH_VARARGS, "..."},
    {NULL, NULL, 0, NULL}
};

static int FMM_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int FMM_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}


static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "FMM3D",
        NULL,
        sizeof(struct module_state),
        FMMMethods,
        NULL,
        FMM_traverse,
        FMM_clear,
        NULL
};

extern "C" PyObject * PyInit_FMM3D(void)
{
	PyObject *module = PyModule_Create(&moduledef);
	if (module == NULL)
        return NULL;
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException("FMM3D.Error", NULL, NULL);
    if (st->error == NULL) 
    {
        Py_DECREF(module);
        return NULL;
    }
	import_array();
	Py_INCREF(module);
    return module;

}