#ifndef MathDenseMatrix_hpp__
#define MathDenseMatrix_hpp__

#include "MathDenseStorage.hpp"

MATH_START
MATH_MATRIX_START

// - Predefined functors that can be used by methods that require a position criteria predicate (CopyFromMatrixWithPositionCriteria)
struct MainDiagonal
{
  bool operator()(std::size_t Row, std::size_t Col) const
  {
    return Row == Col;
  }
};

struct AboveMainDiagonal
{
  bool operator()(std::size_t Row, std::size_t Col) const
  {
    return Row < Col;
  }
};

struct BelowMainDiagonal
{
  bool operator()(std::size_t Row, std::size_t Col) const
  {
    return Row > Col;
  }
};

class SecondaryDiagonal
{
  std::size_t m_nSize;  // the square matrix's size
public :
  SecondaryDiagonal(std::size_t nSize) : m_nSize(nSize)
  {
  }
  bool operator()(std::size_t Row, std::size_t Col) const
  {
    return (Row + Col) == m_nSize;
  }
};

template <typename Scalar, std::size_t Rows = Misc::DynamicSize, std::size_t Cols = Misc::DynamicSize>
class DenseMatrix : public DenseStorage<Scalar, Rows, Cols>
{
public :
  typedef DenseMatrix<Scalar, Rows, Cols> CrtSpec;
  enum class FilterType : unsigned char { Max = 0, Min, MaxAbs, MinAbs };

  // - typedefs for use with STL standard algorithms

  typedef Misc::RawIterator<Scalar>          iterator;
  typedef Misc::RawIterator<const Scalar>    const_iterator;
  typedef Misc::RevRawIterator<Scalar>       reverse_iterator;
  typedef Misc::RevRawIterator<const Scalar> const_reverse_iterator;
  typedef Scalar        value_type;
  typedef std::size_t   size_type;
  typedef std::size_t   difference_type;
  typedef Scalar*       pointer;
  typedef Scalar&       reference;
  typedef const Scalar* const_pointer;
  typedef const Scalar& const_reference;
  
public :
  /// Ctor (dimension params have to be specified only for dynamic matrices, for static matrices they are ignored)
  explicit DenseMatrix(const typename DenseMatrix<Scalar, Rows, Cols>::Type& Tolerance, 
                       std::size_t nRows = Rows, std::size_t nCols = Cols);
  /// Copy ctor
  DenseMatrix(const DenseMatrix<Scalar, Rows, Cols>& RHS);
  /// Move ctor
  DenseMatrix(DenseMatrix<Scalar, Rows, Cols>&& RHS);
  /// Dtor
  ~DenseMatrix();

  /// Create a deep copy of the current object
  DenseMatrix<Scalar, Rows, Cols>* MakeCopy() const;

  // - Operators

  /// Copy assignment
  DenseMatrix<Scalar, Rows, Cols>& operator = (const DenseMatrix<Scalar, Rows, Cols>& RHS_Object);
  /// Move assignment
  DenseMatrix<Scalar, Rows, Cols>& operator = (DenseMatrix<Scalar, Rows, Cols>&& RHS_Object);
  /// Test if the 2 matrices are equal (must be of the same type)
  bool operator == (const DenseMatrix<Scalar, Rows, Cols>& RHS_Object) const;
  /// Test if the 2 matrices aren't equal (must be of the same type)
  bool operator != (const DenseMatrix<Scalar, Rows, Cols>& RHS_Object) const;
  /// Addition 
  const DenseMatrix<Scalar, Rows, Cols> operator + (const DenseMatrix<Scalar, Rows, Cols>& RHS_Object) const;
  /// Subtraction
  const DenseMatrix<Scalar, Rows, Cols> operator - (const DenseMatrix<Scalar, Rows, Cols>& RHS_Object) const;
  /// Scalar multiplication on the right
  const DenseMatrix<Scalar, Rows, Cols> operator * (const typename DenseMatrix<Scalar, Rows, Cols>::Type& ScalarValue) const;
  /// Matrix multiplication
  template<std::size_t RowsB, std::size_t ColsB>
  const DenseMatrix<Scalar, Rows, ColsB> operator * (const DenseMatrix<Scalar, RowsB, ColsB>& RHS_Object) const;
  /// Multiply by -1
  const DenseMatrix<Scalar, Rows, Cols> operator - () const;

  // - Accessors

  /// Get an iterator to the beginning (STL compliant)
  typename DenseMatrix<Scalar, Rows, Cols>::iterator begin();
  /// Get a const iterator to the beginning (STL compliant)
  typename DenseMatrix<Scalar, Rows, Cols>::const_iterator cbegin() const;
  /// Get an iterator to the end (STL compliant)
  typename DenseMatrix<Scalar, Rows, Cols>::iterator end();
  /// Get a const iterator to the end (STL compliant)
  typename DenseMatrix<Scalar, Rows, Cols>::const_iterator cend() const;

  /// Get a reverse iterator to the beginning (actually the end) (STL compliant)
  typename DenseMatrix<Scalar, Rows, Cols>::reverse_iterator rbegin();
  /// Get a const reverse iterator to the beginning (actually the end) (STL compliant)
  typename DenseMatrix<Scalar, Rows, Cols>::const_reverse_iterator crbegin() const;
  /// Get a reverse iterator to the end (actually the beginning) (STL compliant)
  typename DenseMatrix<Scalar, Rows, Cols>::reverse_iterator rend();
  /// Get a const reverse iterator to the end (actually the beginning) (STL compliant)
  typename DenseMatrix<Scalar, Rows, Cols>::const_reverse_iterator crend() const;

  /// Get the stored tolerance
  Scalar GetTolerance() const;
  /// Get a dynamic matrix from the current static matrix
  DenseMatrix<Scalar, Misc::DynamicSize, Misc::DynamicSize> CreateDynamicMatrix() const;
  /// Extract a row vector
  DenseMatrix<Scalar, Math::Misc::ValueIfRtDynamic<Rows, Misc::DynamicSize, 1>::Value, Math::Misc::ValueIfRtDynamic<Cols, Misc::DynamicSize, Cols>::Value>
  ExtractRowVector(std::size_t nRowIdx) const;
  /// Extract a column vector (not exactly cache friendly)
  DenseMatrix<Scalar, Math::Misc::ValueIfRtDynamic<Rows, Misc::DynamicSize, Rows>::Value, Math::Misc::ValueIfRtDynamic<Cols, Misc::DynamicSize, 1>::Value>
  ExtractColumnVector(std::size_t nColumnIdx) const;
  /// Copy a submatrix from the current matrix to the specified matrix (either static or dynamic) at the specified location (not cache friendly)
  template <std::size_t RowsSubMatrix, std::size_t ColsSubMatrix>
  void ExtractSubMatrix(DenseMatrix<Scalar, RowsSubMatrix, ColsSubMatrix>& SubMatrix, const Math::Matrix::IdxPair& DestUprLeftCorner,
                        const Math::Matrix::IdxPair& SrcUprLeftCorner, const Math::Matrix::IdxPair& SrcLwrRightCorner) const;

  // - Special matrix operations

  /// Test if the current matrix is symmetric A = t(A), (valid only for square matrices, will not compile for non square static matrices)
  typename Misc::CompileIf<Rows == Cols, bool, void>::RetValue IsSymmetric() const;
  /// Get the transpose
  DenseMatrix<Scalar, Cols, Rows> GetTransposed() const;

  // - Mutators

  /// Copy a set of elements that obey the specific rule from a source matrix to the crt matrix (the no. rows and cols must be identical, predicate can be a functor or a lambda of the form bool F(row, col))
  template<typename Pred>
  DenseMatrix<Scalar, Rows, Cols>& CopyFromMatrixWithPositionCriteria(DenseMatrix<Scalar, Rows, Cols>& Source, Pred Fctor);
  /// Set the current matrix as identity (only for square matrices, will not compile for compile time matrices that aren't square)
  typename Misc::CompileIf<Rows == Cols, DenseMatrix<Scalar, Rows, Cols>&, void>::RetValue SetIdentity();
  /// Set the current matrix equal to a row permutation matrix
  DenseMatrix<Scalar, Rows, Cols>& SetRowPermutationMatrix(const Misc::Permutation<Rows>& RowPermutation);
  /// Set the current matrix equal to a row permutation matrix
  DenseMatrix<Scalar, Rows, Cols>& SetColumnPermutationMatrix(const Misc::Permutation<Cols>& ColPermutation);
  /// Zero out the entire matrix
  DenseMatrix<Scalar, Rows, Cols>& ZeroOut();
  /// Zero out the specified row (cache friendly)
  DenseMatrix<Scalar, Rows, Cols>& ZeroOutRow(std::size_t nRow);
  /// Zero out the specified column (not cache friendly)
  DenseMatrix<Scalar, Rows, Cols>& ZeroOutCol(std::size_t nCol);
  /// Multiply a row by the factor specified
  DenseMatrix<Scalar, Rows, Cols>& MultiplyRow(std::size_t nIdxRow, const typename DenseMatrix<Scalar, Rows, Cols>::Type& Value);
  /// Multiply a column by the factor specified
  DenseMatrix<Scalar, Rows, Cols>& MultiplyCol(std::size_t nIdxCol, const typename DenseMatrix<Scalar, Rows, Cols>::Type& Value);
  /// Add a specified row (multiplied by a scaling coefficient) to another row
  DenseMatrix<Scalar, Rows, Cols>& AddRowJ_toRowI(std::size_t nRowI, std::size_t nRowJ, typename DenseMatrix<Scalar, Rows, Cols>::Type Scale = 1);
  /// Add a specified column (multiplied by a scaling coefficient) to another column
  DenseMatrix<Scalar, Rows, Cols>& AddColJ_toColI(std::size_t nColI, std::size_t nColJ, typename DenseMatrix<Scalar, Rows, Cols>::Type Scale = 1);

public :
  /// Pivot operation on the specified row, determines the index of the row with the maximum value from [nRowPiv, NoRows-1] considering the column nColPiv, returns the row to be exchanged
  std::size_t PivotizeRow(std::size_t nRowPiv, std::size_t nColPiv) const;
  /// Return the abs max value from a matrix column, considering a specific permutation of the same size as the no. of rows in the matrix (in the specified range)
  std::size_t MaxPermutationIndex(const Misc::Permutation<Math::Misc::ValueIfRtDynamic<Rows, Misc::DynamicSize, Rows>::Value>& Permutation, 
                                  const Math::Matrix::Range& PermutationIdxRange, std::size_t nColPiv) const;

protected :
  CrtSpec& RefTable;
  Type m_Tolerance;
};

#include "MathDenseMatrix.inl"

MATH_MATRIX_END
MATH_END

#endif // MathDenseMatrix_hpp__
