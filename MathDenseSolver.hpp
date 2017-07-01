#ifndef MathDenseSolver_hpp__
#define MathDenseSolver_hpp__

#include "MathDenseMatrix.hpp"

MATH_START
MATH_MATRIX_START

enum class ProblemType : unsigned short { Determinant, LinearEqSolve, Decomposition, Inverse };

#define KEEP_MATRIX   0x01        ///< Internally set
#define PARTIAL_PIVOT 0x02        ///< User set

template <typename Scalar, std::size_t Rows = Misc::DynamicSize, std::size_t Cols = Misc::DynamicSize>
class SolverManager
{
public :
  typedef DenseMatrix<Scalar, Rows, Cols> CrtMatrixType;
protected :
  /// Ctor
  SolverManager(DenseMatrix<Scalar, Rows, Cols>& CrtObject, bool bKeepMatrix = true);
public :
  /// Dtor
  ~SolverManager();

  // - Accessors

  /// Test if the specified option is enabled
  bool OptionEnabled(unsigned short Flag) const;
  
  // - Mutators
  
  /// Enable the specified flag (applies for all the flags except KEEP_MATRIX)
  SolverManager<Scalar, Rows, Cols>& EnableOption(unsigned short Flag, bool bEnable = true);
public :
  /// Upload a matrix (either a copy or a reference to the original matrix)
  void UploadMatrix(DenseMatrix<Scalar, Rows, Cols>& CrtObject, bool bKeepMatrix = true);
protected :
  /// Release matrix
  SolverManager<Scalar, Rows, Cols>& ReleaseMatrix();
  
protected :
  DenseMatrix<Scalar, Rows, Cols>* m_pCrtMatrix;
  unsigned short m_nOptions;

private :
  /// Copy ctor (blocked copy construction)
  SolverManager(const SolverManager<Scalar, Rows, Cols>& RHS);
  /// Copy assignment (blocked copy assignment)
  SolverManager<Scalar, Rows, Cols>& operator = (const SolverManager<Scalar, Rows, Cols>& RHS);
  /// Move assignment (blocked move assignment)
  SolverManager<Scalar, Rows, Cols>& operator = (SolverManager<Scalar, Rows, Cols>&& RHS);
};

/// Blank specialization, specializations will carry out the entire work
template <ProblemType RunPurpose, typename Scalar, std::size_t Rows = Misc::DynamicSize, std::size_t Cols = Misc::DynamicSize>
struct Problem : public SolverManager<Scalar, Rows, Cols>
{
};

/// Partial specialization for the base functor that performs a Lower Upper decomposition for a square matrix (either using Crout's or Doolittle's approach)
template<typename Scalar, std::size_t Rows, std::size_t Cols>
struct Problem<ProblemType::Decomposition, Scalar, Rows, Cols> : public SolverManager<Scalar, Rows, Cols>
{
  /// Ctor
  explicit Problem(DenseMatrix<Scalar, Rows, Cols>& CrtObject, bool bKeepMatrix = true);

  enum class LU_Strategy : unsigned char { Crout, Doolittle } eStrategy;

  /// Perform the actual calculation (returns true if the calculation can be carried out successfully)
  typename Misc::CompileIf<Rows == Cols, bool, void>::RetValue Solve();

  // - Results

  std::size_t                    PermutationsPerformed;   ///< The number of permutations performed
  Misc::Permutation<Rows>        RowPermutation;          ///< The row permutation performed
  DenseMatrix<Scalar, Rows, Cols> LUDecMatrix;            ///< the actual decomposition (For Crout the lower and upper matrices are stored in LUDecMatrix below/above the main diagonal, note that the Upr matrix has the main diagonal filled with 1s in theory), returns false if the matrix is singular)
                                                          ///< For Doolittle the lower and upper matrices are stored in LUDecMatrix below/above the main diagonal, note that the Lwr matrix has the main diagonal filled with 1s in theory
};

/// Partial specialization that computes the determinant of a square matrix (only for square matrices static (compile time dim checks) or dynamic (runtime dim checks))
template<typename Scalar, std::size_t Rows, std::size_t Cols>
struct Problem<ProblemType::Determinant, Scalar, Rows, Cols> : public SolverManager<Scalar, Rows, Cols>
{
  /// Ctor
  explicit Problem(DenseMatrix<Scalar, Rows, Cols>& CrtObject, bool bKeepMatrix = true);

  enum class Det_Strategy : unsigned char { Gauss, LU_Crout, LU_Doolittle } eStrategy;

  /// Perform the actual calculation (returns true if the calculation can be carried out successfully)
  typename Misc::CompileIf<Rows == Cols, bool, void>::RetValue Solve();

  // - Results

  std::size_t             PermutationsPerformed;       ///< The number of permutations performed
  Misc::Permutation<Rows> RowPermutation;              ///< The row permutation performed
  Scalar                  Determinant;                 ///< The determinant value
};

/// Partial specialization that computes the inverse of a square matrix
template<typename Scalar, std::size_t Rows, std::size_t Cols>
struct Problem<ProblemType::Inverse, Scalar, Rows, Cols> : public SolverManager<Scalar, Rows, Cols>
{
  /// Ctor
  explicit Problem(DenseMatrix<Scalar, Rows, Cols>& CrtObject, bool bKeepMatrix = true);

  enum class Inv_Strategy : unsigned char { Gauss_Jordan } eStrategy;

  /// Perform the actual calculation (returns true if the calculation can be carried out successfully)
  typename Misc::CompileIf<Rows == Cols, bool, void>::RetValue Solve();

  // - Results

  std::size_t             PermutationsPerformed;       ///< The number of permutations performed
  Misc::Permutation<Rows> RowPermutation;              ///< The row permutation performed
  DenseMatrix<Scalar, Rows, Cols> Inverse;             ///< The inverted matrix
};

/// Partial specialization that determines the solutions of a linear system of equations
template<typename Scalar, std::size_t Rows, std::size_t Cols>
struct Problem<ProblemType::LinearEqSolve, Scalar, Rows, Cols> : public SolverManager<Scalar, Rows, Cols>
{
  /// Ctor
  explicit Problem(DenseMatrix<Scalar, Rows, Cols>& CrtObject, bool bKeepMatrix = true);

  enum class Res_Strategy : unsigned char { Gauss } eStrategy;
 
  /// Perform the actual calculation (returns true if the calculation can be carried out successfully)
  typename Misc::CompileIf<(Rows + 1) == Cols || (Rows == Misc::DynamicSize && Cols == Misc::DynamicSize), bool, void>::RetValue Solve();

  // - Results

  std::size_t             PermutationsPerformed;                                                                                                               ///< The number of permutations performed
  Misc::Permutation<Rows> RowPermutation;                                                                                                                      ///< The row permutation performed
  DenseMatrix<Scalar, Misc::ValueIfRtDynamic<Rows, Misc::DynamicSize, Rows>::Value, Misc::ValueIfRtDynamic<Cols, Misc::DynamicSize, 1>::Value> Solutions;      ///< The actual solutions vector
};

#include "MathDenseSolver.inl"

MATH_MATRIX_END
MATH_END

#endif // MathDenseSolver_h__
