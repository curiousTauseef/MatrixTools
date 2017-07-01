#ifndef MathMatrix_hpp__
#define MathMatrix_hpp__

#include "MathExceptions.hpp"
#include "MathBasic.hpp"
#include "MathIterators.hpp"
#include "MathHelpers.hpp"
#include <cstddef>
#include <utility>

MATH_START
MATH_MATRIX_START

typedef std::pair<std::size_t, std::size_t> Range;      ///< Used to define ranges between row,i to row,j (or to columns)
typedef std::pair<std::size_t, std::size_t> IdxPair;    ///< Storage access pair for rows and columns (useful for identifying element (i,j))

template <typename T, std::size_t Rows, std::size_t Cols>
class DS_2DBase
{
public :
  typedef T Type;                                  ///< The current type used to instantiate the generic class
  typedef DS_2DBase<T, Rows, Cols> CrtSpec;        ///< The actual type of the current object
protected :
  /// Default ctor
  DS_2DBase();
  /// Copy ctor
  DS_2DBase(const DS_2DBase<T, Rows, Cols>& RHS);
  /// Move ctor
  DS_2DBase(DS_2DBase<T, Rows, Cols>&& RHS);
public :
  /// Dtor
  virtual ~DS_2DBase();

protected :
  /// Get the address of the element stored at the specified offset
  T* GetAddress(std::size_t Offset);

  // - Operators

protected :
  /// Copy assignment
  DS_2DBase<T, Rows, Cols>& operator = (const DS_2DBase<T, Rows, Cols>& RHS_Object);
  /// Move assignment
  DS_2DBase<T, Rows, Cols>& operator = (DS_2DBase<T, Rows, Cols>&& RHS_Object);
public :
  /// Access row and line (like standard 2D contiguous matrices mat[i][j])
  T* operator [] (std::size_t nIdx);
  /// Access row and line (like standard 2D contiguous matrices mat[i][j])
  const T* operator [] (std::size_t nIdx) const;
  /// Access row and line (the individual elements)
  T& operator () (std::size_t nIdxRow, std::size_t nIdxCol);
  /// Access row and line (the individual elements)
  const T& operator () (std::size_t nIdxRow, std::size_t nIdxCol) const;
  /// Linearly access element (the individual elements)
  T& operator () (std::size_t nIdx);
  /// Linearly access element(the individual elements)
  const T& operator () (std::size_t nIdx) const;
  /// Insert element into matrix (circular buffer type)
  DS_2DBase<T, Rows, Cols>& operator << (const T& Value);

  // - Mutators
  /// Swap the contents
  void SwapContents(const DS_2DBase<T, Rows, Cols>& Source);
  /// Copy data into the specified row
  void CopyOnRow(std::size_t nRow, T Array[], std::size_t nArraySize);
  /// Copy data into the specified row (the number of columns must be the same)   
  void CopyOnRow(const DS_2DBase<T, Misc::ValueIfRtDynamic<Cols, Misc::DynamicSize, 1>::Value, Misc::ValueIfRtDynamic<Cols, Misc::DynamicSize, Cols>::Value>& RowLine, std::size_t nRow);
  /// Copy data into the specified column
  void CopyOnCol(std::size_t nCol, T Array[], std::size_t nArraySize);
  /// Copy data into the specified column (the number of rows must be the same)
  void CopyOnCol(const DS_2DBase<T, Misc::ValueIfRtDynamic<Rows, Misc::DynamicSize, Rows>::Value, Misc::ValueIfRtDynamic<Rows, Misc::DynamicSize, 1>::Value>& ColLine, std::size_t nCol);
  /// Copy data from the specified massive
  void CopyFromMassive(T Array[], std::size_t nArraySize);
  /// Swap 2 rows
  void SwapRows(std::size_t nRowI, std::size_t nRowJ);
  /// Swap 2 columns
  void SwapCols(std::size_t nColI, std::size_t nColJ);
  /// Swap 2 cells
  void SwapCells(const Math::Matrix::IdxPair& Cell_1, const Math::Matrix::IdxPair& Cell_2);
  /// Apply a row permutation to the current storage space
  void ApplyRowPermutation(const Math::Misc::Permutation<Rows>& RowPermutation);
  /// Apply a column permutation to the current storage space
  void ApplyColPermutation(const Math::Misc::Permutation<Cols>& ColPermutation);
  /// Apply a simple function object, function pointer or a lambda to each element in the storage (prototype should be void func(Type& Val))
  template <typename F>
  void ProcessEachElement(F Function);
  /// Apply a simple function object, function pointer or a lambda to each element in the storage (prototype should be void func(Type& Val))
  template <typename F>
  void ProcessRow(F Function, std::size_t nIdxRow);
  /// Apply a simple function object, function pointer or a lambda to each element in the storage (prototype should be void func(Type& Val))
  template <typename F>
  void ProcessCol(F Function, std::size_t nIdxCol);

  // - Accessors
  /// Get the number of rows
  virtual std::size_t GetNoRows() const = 0;
  /// Get the number of columns
  virtual std::size_t GetNoCols() const = 0;
  /// Get the linear size of the internal massive
  virtual std::size_t GetLinearSize() const = 0;
  /// Serialize data into an ostream
  std::ostream& Serialize(std::ostream& Stream) const;
  /// Copy data to the specified massive (internal allocation)
  void CopyToMassive(T*& Array, std::size_t& nSize) const;

  // - Internal members
protected :
  Type* m_pDataMassive;   ///< actual data massive 
  std::size_t m_nCrtIdx;  ///< current position (used when pushing elements inside with operator <<)
};

//////////////////////////////////////////////////////////////////////////
/// @ brief 2D Dense storage with compile time number of rows and columns
///
/// By not specifying the number of rows and columns the partial specialization will be used
/// Example of usage Math::Storage::DenseStorage<double, 5, 5> 2DTable;
//////////////////////////////////////////////////////////////////////////
template <typename T, std::size_t Rows = Misc::DynamicSize, std::size_t Cols = Misc::DynamicSize>
class DenseStorage : public DS_2DBase<T, Rows, Cols>
{
protected :
  /// Ctor (no need to provide these arguments in the case of a static element, only here to have a common interface)
  DenseStorage(std::size_t nRows = Misc::DynamicSize, std::size_t nCols = Misc::DynamicSize);
  /// Copy ctor
  DenseStorage(const DenseStorage<T, Rows, Cols>& RHS);
  /// Move ctor
  DenseStorage(DenseStorage<T, Rows, Cols>&& RHS);
public :
  /// Dtor
  virtual ~DenseStorage();

public :
  /// Create a deep copy of the current object
  DenseStorage<T, Rows, Cols>* MakeCopy() const;

  // - Operators

protected :
  /// Copy assignment
  DenseStorage<T, Rows, Cols>& operator = (const DenseStorage<T, Rows, Cols>& RHS_Object);
  /// Move assignment
  DenseStorage<T, Rows, Cols>& operator = (DenseStorage<T, Rows, Cols>&& RHS_Object);
public :

  // - Accessors
  /// Get the number of rows
  std::size_t GetNoRows() const;
  /// Get the number of columns
  std::size_t GetNoCols() const;
  /// Get the linear size of the internal massive
  std::size_t GetLinearSize() const;
  /// Create a dynamic (with reshape capabilities) storage object from the existing static matrix
  DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize> CreateRTDynamicStorage() const;
};

//////////////////////////////////////////////////////////////////////////
/// @brief 2D Dense storage with runtime time number of rows and columns (partial specialization)
///
/// The storage is dynamic which basically means that it can be reshaped at runtime
/// Example of usage Math::Storage::DenseStorage<double> 2DTable(5, 5);
//////////////////////////////////////////////////////////////////////////
template <typename T>
class DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize> : public DS_2DBase<T, Misc::DynamicSize, Misc::DynamicSize>
{
protected :
  /// Ctor for square matrices
  explicit DenseStorage(std::size_t nSize);
  /// Ctor for rectangular matrices
  DenseStorage(std::size_t nRows, std::size_t nCols);
  /// Copy ctor
  DenseStorage(const DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>& RHS);
  /// Move ctor
  DenseStorage(DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>&& RHS);
public :
  /// Dtor
  virtual ~DenseStorage();

public :
  /// Create a deep copy of the current object
  DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>* MakeCopy() const;

  // - Operators

public :
  /// Copy assignment
  DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>& operator = (const DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>& RHS_Object);
  /// Move assignment
  DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>& operator = (DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>&& RHS_Object);

  // - Mutators
  /// Set the size for the current matrix (all info is kept, internal counters are reset)
  void Reshape(std::size_t nRowsNo, std::size_t nColsNo);
  /// Append a row to the current table (reshape is called and a deep copy is performed) (the added line has to be dynamic)
  void AppendRow(const DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>& Row);
  /// Append a row to the current table (reshape is called and a deep copy is performed) (the added line has to be static)
  template <std::size_t Cols>
  void AppendRow(const DenseStorage<T, 1, Cols>& Row);
  /// Append a column to the current table (reshape is called) (the added line has to be dynamic)
  void AppendCol(const DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>& Col);
  /// Append a column to the current table (reshape is called and a deep copy is performed) (the added line has to be static)
  template <std::size_t Rows>
  void AppendCol(const DenseStorage<T, Rows, 1>& Col);
  
  // - Accessors
  /// Get the number of rows
  std::size_t GetNoRows() const;
  /// Get the number of columns
  std::size_t GetNoCols() const;
  /// Get the size of the massive
  std::size_t GetLinearSize() const;

  // - Internal members
protected :
  std::size_t m_nNoRows;  ///< Number of rows
  std::size_t m_nNoCols;  ///< Number of cols
};

//////////////////////////////////////////////////////////////////////////
/// @brief 1D Dense storage with compile time number of columns (partial specialization)
///
/// The storage is actually unidimensional in the sense that there's a single row and multiple columns (row vector)
/// Example of usage Math::Storage::DenseStorage<double, 1, 5> Vector();
//////////////////////////////////////////////////////////////////////////
template <typename T, std::size_t Cols>
class DenseStorage<T, 1, Cols> : public DS_2DBase<T, 1, Cols>
{
protected :
  /// Ctor (no need to provide these arguments in the case of a static element, only here to have a common interface)
  DenseStorage(std::size_t nRows = 1, std::size_t nCols = Cols);
  /// Copy ctor
  DenseStorage(const DenseStorage<T, 1, Cols>& RHS);
  /// Move ctor
  DenseStorage(DenseStorage<T, 1, Cols>&& RHS);
public :
  /// Dtor
  virtual ~DenseStorage();

protected :
  /// Create a deep copy of the current object
  DenseStorage<T, 1, Cols>* MakeCopy() const;

  // - Operators

protected :
  /// Copy assignment
  DenseStorage<T, 1, Cols>& operator = (const DenseStorage<T, 1, Cols>& RHS_Object);
  /// Move assignment
  DenseStorage<T, 1, Cols>& operator = (DenseStorage<T, 1, Cols>&& RHS_Object);
public :
  /// Access column (like standard 1D contiguous arrays mat[i])
  T& operator [] (std::size_t nIdx);
  /// Access column (like standard 1D contiguous arrays mat[i])
  const T& operator [] (std::size_t nIdx) const;
  /// Access column (the individual elements)
  T& operator () (std::size_t nIdxCol);
  /// Access column (the individual elements)
  const T& operator () (std::size_t nIdxCol) const;

  // - Accessors
  /// Convert to column vector (transposed form)
  DenseStorage<T, Cols, 1> GetColumnVector();
  /// Get the number of rows
  std::size_t GetNoRows() const;
  /// Get the number of columns
  std::size_t GetNoCols() const;
  /// Get the size of the massive
  std::size_t GetLinearSize() const;
  /// Create a dynamic (with reshape capabilities) storage object from the existing static matrix
  DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize> CreateRTDynamicStorage() const;
};

//////////////////////////////////////////////////////////////////////////
/// @brief 1D Dense storage with compile time number of rows (partial specialization)
///
/// The storage is actually unidimensional in the sense that there's a single column and multiple rows (column vector)
/// Example of usage Math::Storage::DenseStorage<double, 5, 1> Vector();
//////////////////////////////////////////////////////////////////////////
template <typename T, std::size_t Rows>
class DenseStorage<T, Rows, 1> : public DS_2DBase<T, Rows, 1>
{
protected :
  /// Ctor (no need to provide these arguments in the case of a static element, only here to have a common interface)
  DenseStorage(std::size_t nRows = Rows, std::size_t nCols = 1);
  /// Copy ctor
  DenseStorage(const DenseStorage<T, Rows, 1>& RHS);
  /// Move ctor
  DenseStorage(DenseStorage<T, Rows, 1>&& RHS);
public :
  /// Dtor
  virtual ~DenseStorage();

protected :
  /// Create a deep copy of the current object
  DenseStorage<T, Rows, 1>* MakeCopy() const;

  // - Operators

protected :
  /// Copy assignment
  DenseStorage<T, Rows, 1>& operator = (const DenseStorage<T, Rows, 1>& RHS_Object);
  /// Move assignment
  DenseStorage<T, Rows, 1>& operator = (DenseStorage<T, Rows, 1>&& RHS_Object);
public :
  /// Access column (like standard 1D contiguous arrays mat[i])
  T& operator [] (std::size_t nIdx);
  /// Access column (like standard 1D contiguous arrays mat[i])
  const T& operator [] (std::size_t nIdx) const;
  /// Access column (the individual elements)
  T& operator () (std::size_t nIdxRow);
  /// Access column (the individual elements)
  const T& operator () (std::size_t nIdxRow) const;

  // - Accessors
  /// Convert to row vector (transposed form)
  DenseStorage<T, 1, Rows> GetRowVector();
  /// Get the number of rows
  std::size_t GetNoRows() const;
  /// Get the number of columns
  std::size_t GetNoCols() const;
  /// Get the size of the massive
  std::size_t GetLinearSize() const;
  /// Create a dynamic (with reshape capabilities) storage object from the existing static matrix
  DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize> CreateRTDynamicStorage() const;
};

// - For converting 2 indices into a linear offset
#define ACC(i, j, COLS_NO) i * COLS_NO + j

#include "MathDenseStorage.inl"

MATH_MATRIX_END
MATH_END


#endif // MathMatrix_h__
