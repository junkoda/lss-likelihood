#include <vector>
#include <cmath>
#include <cassert>

#include "py_model.h"
#include "py_error.h"
#include "py_assert.h"

using namespace std;

//
//
//
class Model {
 public:
  virtual ~Model();
  virtual void evaluate(const double param[]) const = 0;

};

class Kaiser : public Model {
 public:
  Kaiser(PyObject* py_k, PyObject* py_P);
  virtual void evaluate(const double param[]) const;
private:
  vector<double> v_k, v_P;
};

//
// static function
//
static void decode_array(const char name[],
			 PyObject* py_obj, Py_buffer* buf,
			 Py_ssize_t len=0)
{
  // name: name of the array for error message
  // py_obs: array object
  // buf: resulting buffer object
  // len: expected length; raise error if the length of the array is not len
  
  char msg[128];
  
  if(PyObject_GetBuffer(py_obj, buf, PyBUF_FORMAT | PyBUF_FULL_RO) == -1)
    throw TypeError();
  
  if(buf->ndim != 1) {
    sprintf(msg, "Expected a 1-dimensional array for %s", name);
    PyErr_SetString(PyExc_TypeError, msg);
    throw TypeError();
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
}

static void array_as_vector(const char name[],
			    PyObject* py_obj,
			    vector<double>& v,
			    const Py_ssize_t len=0)
{
  Py_buffer buf;

  decode_array(name, py_obj, &buf);

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


static void py_model_free(PyObject *obj)
{
  // Delete the power spectrum object
  // Called automatically by Python
  Model* const model= (Model*) PyCapsule_GetPointer(obj, "_Model");
  assert(model);

  delete model;
}

//
// Model implementation
//
Model::~Model()
{

}

// !!! deprecated !!!
PyObject* py_model_alloc(PyObject* self, PyObject* args)
{
  return NULL;
}

//
// Kaiser
//
Kaiser::Kaiser(PyObject* py_k, PyObject* py_P)
{
  array_as_vector("k", py_k, v_k);
  array_as_vector("P", py_P, v_P, v_k.size());
}

void Kaiser::evaluate(const double param[]) const
{

}

PyObject* py_model_kaiser_alloc(PyObject* self, PyObject* args)
{
  // _model_kaiser_alloc(k, P)
  PyObject *py_k, *py_P;
    
  if(!PyArg_ParseTuple(args, "OO", &py_k, &py_P))
    return NULL;

  Py_buffer k, P;
  try {
    return PyCapsule_New(new Kaiser(py_k, py_P), "_Model", py_model_free);
  }
  catch(TypeError) {
    return NULL;
  }

  return NULL;
}


