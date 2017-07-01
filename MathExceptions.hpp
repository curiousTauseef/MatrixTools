#ifndef MathExceptions_hpp__
#define MathExceptions_hpp__

#include "MathDefines.h"
#include <exception>
#include <cassert>

MATH_START
MATH_BASIC_START

class GeneralException : std::exception
{
public :
  enum class Type : char { InvalidParams = 0 };
  GeneralException(GeneralException::Type CrtType) : std::exception(ToString(CrtType)) 
  {
  }
private :
  const char* ToString(GeneralException::Type CrtType)
  {
    switch (CrtType)
    {
    case Type::InvalidParams :
      return "InvalidParams";
    default :
      {
        assert(false && L"Not covered...");
        return "GeneralException unknown";
      }
    }
  }
};

class StorageException : std::exception
{
public :
  enum class Type : char { WrongMatrixType = 0, InvalidDimsForOperation, InvalidIndex, InvalidMemoryAcces };
  StorageException(StorageException::Type CrtType) : std::exception(ToString(CrtType)) 
  {
  }
private :
  const char* ToString(StorageException::Type CrtType)
  {
    switch (CrtType)
    {
    case Type::WrongMatrixType :
      return "WrongMatrixType";
    case Type::InvalidDimsForOperation :
      return "InvalidDimsForOperation";
    case Type::InvalidIndex :
      return "InvalidIndex";
    case Type::InvalidMemoryAcces :
      return "InvalidMemoryAcces";
    default :
      {
        assert(false && L"Not covered...");
        return "StorageException unknown";
      }
    }
  }
};

class MatrixException : std::exception
{
public :
  enum class Type : char { NonSquareMatrix = 0, SingularMatrix};
  MatrixException(MatrixException::Type CrtType) : std::exception(ToString(CrtType)) 
  {
  }
private :
  const char* ToString(MatrixException::Type CrtType)
  {
    switch (CrtType)
    {
    case Type::NonSquareMatrix :
      return "NonSquareMatrix";
    case Type::SingularMatrix :
      return "SingularMatrix";
    default :
      {
        assert(false && L"Not covered...");
        return "MatrixException unknown";
      }
    }
  }
};



MATH_BASIC_END
MATH_END

#endif // MathExceptions_hpp__
