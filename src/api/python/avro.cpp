#include <stdio.h>

#include <python.h>

void
py_avro( char* args )
{
  printf("hello from avro!\n");
}

static PyObject*
py_avro_wrapper( PyObject* self , PyObject* args )
{
  char * input;
  char * result;
  PyObject * ret;

  // parse arguments
  if (!PyArg_ParseTuple(args, "s", &input)) {
    return NULL;
  }

  // run the actual function
  result = hello(input);

  // build the resulting string into a Python object.
  ret = PyString_FromString(result);
  free(result);

  return ret;
}

static PyMethodDef avroMethods[] = { {"hello",py_avro_wrapper,METH_VARARGS,"run avro"} , {NULL,NULL,0,NULL} };

DL_EXPORT(void) initavro(void)
{
  Py_InitModule("hello",avroMethods);
}
