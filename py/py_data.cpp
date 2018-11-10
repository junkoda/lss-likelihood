#include <cassert>
#include "py_util.h"
#include "py_data.h"


static void py_data_free(PyObject *obj)
{
  // Delete the data object
  // Called automatically by Python
  Data* const data= (Data*) PyCapsule_GetPointer(obj, "_Data");
  assert(data);

  delete data;
}


PyObject* py_data_alloc(PyObject* self, PyObject* args)
{
  PyObject *py_P0, *py_P2, *py_P4;
  if(!PyArg_ParseTuple(args, "OOO",
		       &py_P0, &py_P2, &py_P4))
    return NULL;

  Data* const data = new Data();

  py_util_array_as_vector("P0", py_P0, data->P0, 0);
  const int len= static_cast<int>(data->P0.size());
  
  py_util_array_as_vector("P2", py_P2, data->P2, len);
  py_util_array_as_vector("P4", py_P4, data->P4, len);

  return PyCapsule_New(data, "_Data", py_data_free);
}
