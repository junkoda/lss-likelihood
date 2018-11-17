#include "py_util.h"
#include "py_error.h"

using namespace std;

static void decode_array(const char name[],
			 PyObject* py_obj, Py_buffer* buf,
			 const Py_ssize_t len=0, const Py_ssize_t ncol=1,
			 const bool read_only=true)
{
  // name: name of the array for error message
  // py_obs: array object
  // buf: resulting buffer object
  // len: expected length; raise error if the length of the array is not len
  
  char msg[128];

  int flag = 0;
  if(read_only)
    flag= PyBUF_FULL_RO;
  else
    flag= PyBUF_FULL;
  
  if(PyObject_GetBuffer(py_obj, buf, PyBUF_FORMAT | flag) == -1)
      throw TypeError();

  if(ncol == 1) {
    if(buf->ndim != 1) {
      sprintf(msg, "Expected a 1-dimensional array for %s", name);
      PyErr_SetString(PyExc_TypeError, msg);
      throw TypeError();
    }
  }
  else if(ncol > 1) {
    if(buf->ndim != 2) {
      sprintf(msg, "Expected a 2-dimensional array for %s", name);
      PyErr_SetString(PyExc_TypeError, msg);
      throw TypeError();
    }
  }


  if(strcmp(buf->format, "d") != 0) {
    sprintf(msg, "Expected an array of double for %s: %s", name, buf->format);
    PyErr_SetString(PyExc_TypeError, msg);
    throw TypeError();
  }

  if(len > 0 && buf->shape[0] != len) {
    sprintf(msg, "Expected the length arrays of %d for %s: %d",
	    (int) len, name, (int) buf->shape[0]);
    PyErr_SetString(PyExc_TypeError, msg);
    throw TypeError();
  }
  
  if(ncol > 0 && buf->shape[1] != ncol) {
    sprintf(msg, "Expected number of columns %d for %s: %d",
	    (int) ncol, name, (int) buf->shape[2]);
    PyErr_SetString(PyExc_TypeError, msg);
    throw TypeError();
  }

}


void py_util_array_as_vector(const char name[],
			     PyObject* py_obj,
			     vector<double>& v,
			     const Py_ssize_t len_expect)
{
  // 1-dimensional array as a vector
  Py_buffer buf;

  decode_array(name, py_obj, &buf, len_expect, 0, true);

  double const * x= (double const *) buf.buf;
  const size_t n= buf.shape[0];
  const size_t stride= buf.strides[0];

  v.reserve(n);
  
  for(size_t i=0; i<n; ++i) {
    v.push_back(*x);
    x = (double const *) ((char const *) x + stride);
  }

  PyBuffer_Release(&buf);
}

size_t py_util_array2_as_vector(const char name[],
				PyObject* py_array,
				vector<double>& v,
				const Py_ssize_t len_expect,
				const Py_ssize_t ncol_expect)
{
  // 2-dimensional array as a vector
  Py_buffer buf;

  decode_array(name, py_array, &buf, len_expect, ncol_expect, true);

  if(buf.ndim != 2) {
    char msg[128];
    sprintf(msg, "Expected a 2-dimensional array for %s: %d",
	    name, (int) buf.ndim);
    PyErr_SetString(PyExc_TypeError, msg);
    throw TypeError();
  }


  double const * x_row= (double const *) buf.buf;
  const size_t n= buf.shape[0];
  const size_t ncol= buf.shape[1];
  const size_t stride_row= buf.strides[0];
  const size_t stride_col= buf.strides[1];

  v.reserve(n*ncol);
  
  for(size_t i=0; i<n; ++i) {
    double const *x = x;
    for(size_t j=0; j<ncol; ++j) {
      v.push_back(*x);

      x = (double const *) ((char const *) x + stride_col); // next column
    }
    x_row = (double const *) ((char const *) x_row + stride_row); // next row
  }

  PyBuffer_Release(&buf);

  return ncol;
}

void py_util_vector_as_array(const char name[], vector<double>& v,
			     PyObject* py_array)
{
  // copy vector v content to array py_array
  // The length of the array must be the same as that of the vector
  Py_buffer buf;

  decode_array(name, py_array, &buf, 0, 0, false);

  double * x= (double *) buf.buf;
  const size_t n= buf.shape[0];
  const size_t stride= buf.strides[0];

  char msg[128];
  if(v.size() != (size_t) buf.shape[0]) {
    sprintf(msg, "Expected the length of arrays of %d for %s: %d",
	    (int) v.size(), name, (int) buf.shape[0]);
    PyErr_SetString(PyExc_TypeError, msg);
    throw TypeError();
  }
  
  for(size_t i=0; i<n; ++i) {
    *x = v[i];
    x = (double *) ((char const *) x + stride);
  }

  PyBuffer_Release(&buf);
}

void py_util_sequence_as_vector(const char name[], PyObject* py_list,
				vector<double>& v)
{
  char msg[128];
  if(PySequence_Check(py_list) == 0) {
    sprintf(msg, "%s is not is not a sequence", name);
    PyErr_SetString(PyExc_TypeError, msg);
    throw TypeError();
  }

  const int len= (int) PySequence_Length(py_list);
  v.reserve(len);
  
  for(int i=0; i<len; ++i) {
    PyObject* py_float= PySequence_GetItem(py_list, i);
    if(!PyFloat_Check(py_float)) {
      sprintf(msg, "%s[%d] is float", name, i);
      PyErr_SetString(PyExc_TypeError, msg);
      throw TypeError();
    }

    v.push_back(PyFloat_AsDouble(py_float));
    Py_DECREF(py_float);
  }
}
