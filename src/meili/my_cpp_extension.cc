#include <boost/python.hpp>
#include <assert.h>
#include <iostream>

using namespace boost::python;

class MyCPPException : public std::exception
{
private:
  std::string message;
  std::string extraData;
public:
  MyCPPException(std::string message, std::string extraData)
  {
    this->message = message;
    this->extraData = extraData;
  }
  const char *what() const throw()
  {
    return this->message.c_str();
  }
  ~MyCPPException() throw()
  {
  }
  std::string getMessage()
  {
    return this->message;
  }
  std::string getExtraData()
  {
    return this->extraData;
  }
};

void my_cpp_function(bool throwException)
{
  std::cout << "Called a C++ function." << std::endl;
  if (throwException)
    {
      throw MyCPPException("Throwing an exception as requested.",
               "This is the extra data.");
    }
}

PyObject *myCPPExceptionType = NULL;

void translateMyCPPException(MyCPPException const &e)
{
  assert(myCPPExceptionType != NULL);
  boost::python::object pythonExceptionInstance(e);
  PyErr_SetObject(myCPPExceptionType, pythonExceptionInstance.ptr());
}

BOOST_PYTHON_MODULE(my_cpp_extension)
{
  boost::python::class_<MyCPPException> myCPPExceptionClass("MyCPPException",
            boost::python::init<std::string, std::string>());
  myCPPExceptionClass.add_property("message", &MyCPPException::getMessage)
    .add_property("extra_data", &MyCPPException::getExtraData);
  myCPPExceptionType = myCPPExceptionClass.ptr();
  boost::python::register_exception_translator<MyCPPException> (&translateMyCPPException);
  boost::python::def("my_cpp_function", &my_cpp_function);
}
