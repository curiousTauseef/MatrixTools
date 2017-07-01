#ifndef MathIterators_hpp__
#define MathIterators_hpp__

// Math Header
#include "MathDefines.h"

MATH_START
MATH_MISC_START

// - Internal trait for adding compile time const types
template<typename Type>
struct AddConst
{
  typedef const Type TypeWithCst;
};

template<typename Type>
struct AddConst<Type const>
{
  typedef Type TypeWithCst;
};

template <typename Type>
class RawIterator : public std::iterator<std::random_access_iterator_tag, Type, std::size_t>
{
public :
  /// Ctor
  RawIterator(Type* Ptr);
  /// Copy ctor
  RawIterator(const RawIterator<Type>& RHS);
  /// Copy ctor from a const object of the current type (const objects will be created from non const objects through this mechanism but not vice versa)
  RawIterator(const RawIterator<typename AddConst<Type>::TypeWithCst>& RHS);
  /// Dtor
  ~RawIterator();

  // - Operators

  /// Assignment from another iterator
  RawIterator<Type>& operator = (const RawIterator<Type>& RHS);
  /// Assignment from another iterator that wraps a const object (const objects will be created from non const objects through this mechanism)
  RawIterator<Type>& operator = (const RawIterator<typename AddConst<Type>::TypeWithCst>& RHS);
  /// Assignment from a pointer to the actual type
  RawIterator<Type>& operator = (Type* TypePtr);
  /// Assignment from a pointer to a const type
  RawIterator<Type>& operator = (typename AddConst<Type>::TypeWithCst* TypePtr);
  /// Test if the crt object wraps a valid pointer
  operator bool() const;
  /// Comparison
  bool operator == (const RawIterator<Type>& Iterator) const;
  /// Complemented comparison
  bool operator != (const RawIterator<Type>& Iterator) const;
  /// Add an offset
  const RawIterator<Type> operator + (const ptrdiff_t& Delta) const;
  /// Subtract an offset
  const RawIterator<Type> operator - (const ptrdiff_t& Delta) const;
  /// Obtain the offset distance from the 2 objects
  ptrdiff_t operator - (const RawIterator<Type>& Iterator) const;
  /// Add an offset and assign
  RawIterator<Type>& operator += (const ptrdiff_t& Delta);
  /// Subtract an offset and assign
  RawIterator<Type>& operator -= (const ptrdiff_t& Delta);
  /// Preincrement
  RawIterator<Type>& operator ++ ();
  /// Predecrement
  RawIterator<Type>& operator -- ();
  /// Postincrement
  RawIterator<Type> operator ++ (int);
  /// Postdecrement
  RawIterator<Type> operator -- (int);
  /// Dereference (internal ptr doesn't get modified) (if the type is const the operation will be read only)
  Type& operator * () const; 
  /// Indirection  (if the type is const the operation will be read only)
  Type* operator -> () const;
  /// Dereference with offset (if the type is const the operation will be read only)
  Type& operator [] (ptrdiff_t Delta) const;

  // - Accesors

  /// Get the currently wrapped pointer 
  Type* GetWrappedPtr();

protected :
  Type* m_CrtPtr;
};


template <typename Type>
class RevRawIterator : public RawIterator<Type>
{
public :
  /// Ctor
  RevRawIterator(Type* Ptr = nullptr);
  /// Copy ctor
  RevRawIterator(const RevRawIterator<Type>& RHS);
  /// Copy ctor from a const object of the current type (const objects will be created from non const objects through this mechanism but not vice versa)
  RevRawIterator(const RevRawIterator<typename AddConst<Type>::TypeWithCst>& RHS);
  /// Dtor
  ~RevRawIterator();

  // - Operators

  /// Assignment from another iterator
  RevRawIterator<Type>& operator = (const RevRawIterator<Type>& RHS);
  /// Assignment from another iterator that wraps a const object (const objects will be created from non const objects through this mechanism)
  RevRawIterator<Type>& operator = (const RevRawIterator<typename AddConst<Type>::TypeWithCst>& RHS);
  /// Assignment from a pointer to the actual type
  RevRawIterator<Type>& operator = (Type* TypePtr);
  /// Assignment from a pointer to a const type
  RevRawIterator<Type>& operator = (typename AddConst<Type>::TypeWithCst* TypePtr);
  /// Comparison
  bool operator == (const RevRawIterator<Type>& Iterator) const;
  /// Complemented comparison
  bool operator != (const RevRawIterator<Type>& Iterator) const;
  /// Add an offset
  const RevRawIterator<Type> operator + (const ptrdiff_t& Delta) const;
  /// Subtract an offset
  const RevRawIterator<Type> operator - (const ptrdiff_t& Delta) const;
  /// Obtain the offset distance from the 2 objects
  ptrdiff_t operator - (const RevRawIterator<Type>& Iterator) const;
  /// Add an offset and assign
  RevRawIterator<Type>& operator += (const ptrdiff_t& Delta);
  /// Subtract an offset and assign
  RevRawIterator<Type>& operator -= (const ptrdiff_t& Delta);
  /// Preincrement
  RevRawIterator<Type>& operator ++ ();
  /// Predecrement
  RevRawIterator<Type>& operator -- ();
  /// Postincrement
  RevRawIterator<Type> operator ++ (int);
  /// Postdecrement
  RevRawIterator<Type> operator -- (int);
  /// Dereference with offset (if the type is const the operation will be read only)
  Type& operator [] (ptrdiff_t Delta) const;
};

#include "MathIterators.inl"

MATH_MISC_END
MATH_END

#endif // MathIterators_hpp__
