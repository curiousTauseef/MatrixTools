// ConsoleApplication1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "MathDenseMatrix.hpp"
#include "MathDenseSolver.hpp"
#include <algorithm>

bool ExtractSubMatrix()
{
  using Math::Matrix::DenseMatrix;

  Math::Matrix::IdxPair SrcLeft(2, 2), SrcRight(3, 3), DestLeft(2, 2);

  // - Test out matrix extraction from a static matrix
  DenseMatrix<double, 4, 4> MatrixSourceStatic(1e-4);
  DenseMatrix<double, 4, 4> SubMatrixDest(1e-4), SubMatrixDestRef(1e-4);

  SubMatrixDest.SetIdentity();
  MatrixSourceStatic << 0 << 1 << 2 << 3
                     << 4 << 5 << 6 << 7
                     << 8 << 9 << 10 << 11
                     << 12 << 13 << 14 << 15;

  SubMatrixDestRef << 1 << 0 << 0 << 0
                   << 0 << 1 << 0 << 0
                   << 0 << 0 << 10 << 11
                   << 0 << 0 << 14 << 15;

  MatrixSourceStatic.ExtractSubMatrix(SubMatrixDest, DestLeft, SrcLeft, SrcRight);

  if (SubMatrixDest != SubMatrixDestRef)
    return false;

  // - Test out matrix extraction from a dynamic matrix into a static matrix (allowed)
  DenseMatrix<double> MatrixDynamic(1e-4, 4, 4);

  MatrixDynamic << 0 << 1 << 2 << 3
                << 4 << 5 << 6 << 7
                << 8 << 9 << 10 << 11
                << 12 << 13 << 14 << 15;
  MatrixDynamic.ExtractSubMatrix(SubMatrixDest, DestLeft, SrcLeft, SrcRight);

  if (SubMatrixDest != SubMatrixDestRef)
    return false;

  return true;
};

bool PerformStaticMatrixOperations()
{
  using Math::Matrix::DenseMatrix;
  DenseMatrix<double, 4, 4> Matrix(1e-4), MatrixRef(1e-4);
  DenseMatrix<double, 1, 4> RowVector(1e-4);
  DenseMatrix<double, 4, 1> ColVector(1e-4);

  MatrixRef << 2 << 2 << 0 << 1
            << 2 << 6 << 0 << 1
            << 0 << 0 << 1 << 1
            << 0 << 0 << 0 << 1;

  Matrix.ZeroOut();
  Matrix.SetIdentity();
  Matrix.ZeroOutCol(0);
  Matrix.ZeroOutRow(1);
  Matrix.SetIdentity();
  Matrix.MultiplyRow(0, 2.);
  Matrix.MultiplyCol(1, 4.);
  Matrix.AddRowJ_toRowI(1, 0, 1.);
  Matrix.AddColJ_toColI(1, 0, 1.);

  RowVector.ProcessEachElement([=](double& x) mutable { x = 1; });
  ColVector = RowVector.GetTransposed();

  Matrix.CopyOnCol(ColVector, 3);

  if (Matrix != MatrixRef)
    return false;

  return true;
}

bool PerformDynamicMatrixOperations()
{
  using Math::Matrix::DenseMatrix;

  DenseMatrix<double> Matrix(1e-7, 4, 4), MatrixRef(1e-7, 4, 4);

  MatrixRef << 2 << 2 << 0 << 0
            << 2 << 6 << 0 << 0
            << 0 << 0 << 1 << 0
            << 0 << 0 << 0 << 1;

  Matrix.ZeroOut();
  Matrix.SetIdentity();
  Matrix.ZeroOutCol(0);
  Matrix.ZeroOutRow(1);
  Matrix.SetIdentity();
  Matrix.MultiplyRow(0, 2.);
  Matrix.MultiplyCol(1, 4.);
  Matrix.AddRowJ_toRowI(1, 0, 1.);
  Matrix.AddColJ_toColI(1, 0, 1.);

  if (Matrix != MatrixRef)
    return false;

  return true;
}

bool PerformIteratorOperations()
{
  using Math::Matrix::DenseMatrix;
  DenseMatrix<double, 4, 4> MatrixSource(1e-7);

  MatrixSource << 0 << 1 << 2 << 3
               << 4 << 5 << 6 << 7
               << 8 << 9 << 10 << 11
               << 12 << 13 << 14 << 15;

  double* Ptr = nullptr;
  ptrdiff_t Diff = 0;
  int i = 0;
  DenseMatrix<double, 4, 4>::iterator itCrt(MatrixSource.begin());
  DenseMatrix<double, 4, 4>::iterator itEnd(MatrixSource.end());

  Ptr = itCrt.GetWrappedPtr();

  if (Ptr != &MatrixSource(0))
    return false;

  Diff = itEnd - itCrt;

  if (Diff != 16)
    return false;

  for (itCrt = MatrixSource.begin(); itCrt != MatrixSource.end(); ++itCrt)
    *itCrt = 1.;

  for (std::size_t j = 0; j < MatrixSource.GetLinearSize(); ++j)
  {
    if (Math::Basic::ValueCmp(MatrixSource(j), 1., 1e-4) != 0)
      return false;
  }

  itCrt = MatrixSource.begin();
  for (std::size_t i = 0; i < MatrixSource.GetLinearSize(); ++i)
    itCrt[i] = 2;

  for (std::size_t j = 0; j < MatrixSource.GetLinearSize(); ++j)
  {
    if (Math::Basic::ValueCmp(MatrixSource(j), 2., 1e-4) != 0)
      return false;
  }

  itCrt = MatrixSource.end();
  DenseMatrix<double, 4, 4>::iterator itNew = itCrt - 3;
  itNew -= 1;
  itNew += 1;
  itNew++;
  itNew--;

  *itNew = 10.;
  itNew = Ptr;    // test a pointer assignment
  *itNew = 5.;

  if (Math::Basic::ValueCmp(MatrixSource(13), 10., 1e-4) != 0 || Math::Basic::ValueCmp(MatrixSource(0), 5., 1e-4) != 0)
    return false;

  DenseMatrix<double, 4, 4>::iterator       itCrt2 = (itCrt - 3) + 1;
  DenseMatrix<double, 4, 4>::const_iterator itCst(MatrixSource.begin());

  // DenseMatrix<double, 4, 4>::iterator itNew2(itCst); // try and create an iterator from a const iterator...compilation should fail (intended)
  // itNew = itCst;                                     // try and assign to an iterator from a const iterator...compilation should fail (intended)

  double dfValue = 6;
  DenseMatrix<double, 4, 4>::iterator itCrtFind = std::find_if(MatrixSource.begin(), MatrixSource.end(), [=](double x) -> bool { return x == dfValue; });

  if (!itCrtFind)
    return false;

  std::size_t k = 0;
  std::for_each(MatrixSource.begin(), MatrixSource.end(), [=](double& x) mutable { x = (k++) % 2 == 0 ? 0. : 1.; });

  DenseMatrix<double, 4, 4> MatRef(1e-5);

  MatRef << 0 << 1 << 0 << 1
         << 0 << 1 << 0 << 1
         << 0 << 1 << 0 << 1
         << 0 << 1 << 0 << 1;

  if (MatRef != MatrixSource)
    return false;

  return true;
}

void PerformConstIteratorOperations()
{
  using Math::Matrix::DenseMatrix;
  DenseMatrix<double, 4, 4> MatrixSource(1e-7);

  MatrixSource << 0 << 1 << 2 << 3
               << 4 << 5 << 6 << 7
               << 8 << 9 << 10 << 11
               << 12 << 13 << 14 << 15;

  double* Ptr = nullptr;
  const double* CPtr = nullptr;
  ptrdiff_t Diff = 0;
  int i = 0;
  DenseMatrix<double, 4, 4>::const_iterator itCrt(MatrixSource.begin());  // test the construction of a const_iterator from an iterator
  DenseMatrix<double, 4, 4>::const_iterator itEnd(MatrixSource.end());

  // Ptr = itCrt.GetWrappedPtr();   // this should fail to compile (intended)
  CPtr = itCrt.GetWrappedPtr();

  Diff = itEnd - itCrt;

  double dfVal = 0.;
  for (itCrt = MatrixSource.begin(); itCrt != MatrixSource.end(); ++itCrt)
  {
    // *itCrt = 1.;                  // this should fail to compile (intended)
    dfVal = *itCrt;
  }

  itCrt = MatrixSource.begin();
  //for (std::size_t i = 0; i < Matrix4s.GetLinearSize(); ++i)
  // itCrt[i] = 2;                  // this should fail to compile (intended)     

  itCrt = MatrixSource.end();
  DenseMatrix<double, 4, 4>::const_iterator itNew = itCrt - 3;
  itNew -= 1;
  itNew += 1;
  itNew++;
  itNew--;

  itNew = CPtr;    // test a pointer assignment

                   //DenseMatrix<double, 4, 4>::iterator itCrt2 = (itCrt - 3) + 1;      // this should fail to compile since you can't construct/assign an iterator from a const iterator
                   //DenseMatrix<double, 4, 4>::const_iterator itCst(Matrix4s.begin());

                   // DenseMatrix<double, 4, 4>::iterator itNew2(itCst);                // try and create an iterator from a const iterator...compilation should fail (intended)
                   // itNew = itCst;                                                    // try and assign to an iterator from a const iterator...compilation should fail (intended)*/
}

void PerformReverseIteratorOperations()
{
  using Math::Matrix::DenseMatrix;
  DenseMatrix<double, 4, 4> MatrixSource(1e-7);

  MatrixSource << 0 << 1 << 2 << 3
               << 4 << 5 << 6 << 7
               << 8 << 9 << 10 << 11
               << 12 << 13 << 14 << 15;

  double* Ptr = nullptr;
  ptrdiff_t Diff = 0;
  int i = 0;
  DenseMatrix<double, 4, 4>::reverse_iterator itCrt(MatrixSource.rbegin());
  DenseMatrix<double, 4, 4>::reverse_iterator itEnd(MatrixSource.rend());

  Ptr = itCrt.GetWrappedPtr();

  Diff = itEnd - itCrt;

  for (itCrt = MatrixSource.rbegin(); itCrt != MatrixSource.rend(); ++itCrt)
    *itCrt = 1.;

  itCrt = MatrixSource.rbegin();
  for (std::size_t i = 0; i < MatrixSource.GetLinearSize(); ++i)
    itCrt[i] = 2;

  itCrt = MatrixSource.rend();
  DenseMatrix<double, 4, 4>::reverse_iterator itNew = itCrt + 3;
  itNew += 1;
  itNew -= 1;
  itNew--;
  itNew++;

  itNew = Ptr;                                                  // test a pointer assignment

  DenseMatrix<double, 4, 4>::reverse_iterator       itCrt2 = (itCrt - 3) + 1;
  DenseMatrix<double, 4, 4>::const_reverse_iterator itCst(MatrixSource.rbegin());

  // DenseMatrix<double, 4, 4>::reverse_iterator itNew2(itCst); // try and create an iterator from a const iterator...compilation should fail (intended)
  // itNew = itCst;                                             // try and assign to an iterator from a const iterator...compilation should fail (intended)
}

bool SolveGauss22()
{
  using namespace Math::Matrix;

  // - Using an extended static matrix and copying the free terms into the coefficients matrix
  {
    DenseMatrix<double, 2, 3> MatrixCoeffs(1e-7);
    DenseMatrix<double, 2, 1> MatrixFreeTerms(1e-7);
    DenseMatrix<double, 2, 1> MatrixRef(1e-7);

    MatrixCoeffs << 2 << 3 << 0
      << 4 << 9 << 0;
    MatrixFreeTerms << 6
      << 15;

    MatrixRef << 1.5
      << 1;

    MatrixCoeffs.CopyOnCol(MatrixFreeTerms, 2);


    Problem<ProblemType::LinearEqSolve, double, 2, 3> SolverInstance(MatrixCoeffs, true);
    SolverInstance.Solve();

    if (SolverInstance.Solutions != MatrixRef)
      return false;
  }

  // - Using dynamic matrices by appending the free terms column
  {
    DenseMatrix<double> MatrixCoeffs(1e-7, 2, 2);
    DenseMatrix<double> MatrixFreeTerms(1e-7, 2, 1);
    DenseMatrix<double> MatrixSols(1e-7, 2, 1);
    DenseMatrix<double> MatrixRef(1e-7, 2, 1);

    MatrixCoeffs << 2 << 3
      << 4 << 9;
    MatrixFreeTerms << 6
      << 15;

    MatrixRef << 1.5
      << 1;

    MatrixCoeffs.AppendCol(MatrixFreeTerms);

    Problem<ProblemType::LinearEqSolve, double> SolverInstance(MatrixCoeffs);
    SolverInstance.Solve();

    if (SolverInstance.Solutions != MatrixRef)
      return false;
  }

  // - Using an extended dynamic matrix and copying the free terms into the coefficients matrix
  {
    DenseMatrix<double> MatrixCoeffs(1e-7, 2, 3);
    DenseMatrix<double> MatrixFreeTerms(1e-7, 2, 1);
    DenseMatrix<double> MatrixRef(1e-7, 2, 1);

    MatrixCoeffs << 2 << 3 << 0
      << 4 << 9 << 0;

    MatrixFreeTerms << 6
      << 15;

    MatrixRef << 1.5
      << 1;

    MatrixCoeffs.CopyOnCol(MatrixFreeTerms, 2);

    Problem<ProblemType::LinearEqSolve, double> SolverInstance(MatrixCoeffs);
    SolverInstance.Solve();

    if (SolverInstance.Solutions != MatrixRef)
      return false;
  }

  return true;
}

bool SolveGauss55()
{
  using namespace Math::Matrix;

  DenseMatrix<double> MatrixCoeffs(1e-7, 5, 6), MatrixSols(1e-7, 5, 1), MatrixRef(1e-7, 5, 1);

  MatrixCoeffs << 2 << 0 << 4 << 5 << 6 << -1
               << 4 << 9 << 8 << 7 << 2 << -3
               << 2 << 0 << 2 << 6 << 9 << -7
               << 1 << 0 << 6 << 5 << 3 << -9
               << 3 << 0 << 7 << 2 << 1 << -14;

  MatrixRef << 315.666666666666 << -154.74074074074048 << -190.66666666666632 << 301.66666666666612 << -229.66666666666623;

  Problem<ProblemType::LinearEqSolve, double> SolverInstance(MatrixCoeffs);
  SolverInstance.Solve();

  if (SolverInstance.Solutions != MatrixRef)
    return false;

  return true;
}

bool ReshapeAndTransposeDynamicMatrix()
{
  using Math::Matrix::DenseMatrix;

  DenseMatrix<double> Matrix(1e-7, 4, 4);
  DenseMatrix<double> RowVector(1e-7, 1, 4);
  DenseMatrix<double> ColVector(1e-7, 5, 1);

  double dfMultiplicationFactor = 3;

  Matrix.SetIdentity();
  Matrix.ProcessRow([=](double& x) { x *= dfMultiplicationFactor; }, 0);
  Matrix.ProcessCol([=](double& x) { x += dfMultiplicationFactor; }, 0);

  RowVector.ProcessEachElement([=](double& x) mutable { x = dfMultiplicationFactor; dfMultiplicationFactor += 1.; });
  ColVector.ProcessEachElement([=](double& x) mutable { x = 2 * dfMultiplicationFactor; dfMultiplicationFactor += 1.; });

  Matrix.AppendCol(RowVector.GetTransposed());
  Matrix.AppendRow(ColVector.GetTransposed());

  Matrix.Reshape(2, Matrix.GetNoCols() - 1);

  DenseMatrix<double> MatrixTr(1e-7, Matrix.GetNoCols(), Matrix.GetNoRows());
  MatrixTr = Matrix.GetTransposed();

  return true;
}

bool TestEquality()
{
  using Math::Matrix::DenseMatrix;

  DenseMatrix<double, 4, 4> Matrix(1e-7, 4, 4);
  DenseMatrix<double, 4, 4> MatrixTr(1e-7, 4, 4);

  Matrix.SetIdentity();
  MatrixTr = Matrix.GetTransposed();

  return Matrix == MatrixTr;
}

bool InvertMatrix22_1()
{
  using namespace Math::Matrix;

  DenseMatrix<double> MatrixCoeffs(1e-7, 2, 2);

  MatrixCoeffs << 2 << 3
               << 4 << 9;

  {
    Problem<ProblemType::Determinant, double> SolverInstance(MatrixCoeffs);
    SolverInstance.Solve();

    if (Math::Basic::ValueCmp(SolverInstance.Determinant, 6., 1e-4) != 0)
      return false;
  }

  {
    DenseMatrix<double> MatrixInvRef(1e-7, 2, 2);
    Problem<ProblemType::Inverse, double> SolverInstance(MatrixCoeffs);
    SolverInstance.Solve();

    MatrixInvRef << 1.5 << -0.5
                 << -0.666666666666 << 0.333333333333;

    if (SolverInstance.Inverse != MatrixInvRef)
      return false;
  }

  return true;
}

bool InvertMatrix22_2()
{
  using namespace Math::Matrix;

  DenseMatrix<double> MatrixCoeffs(1e-7, 2, 2);

  MatrixCoeffs << 2 << 3
               << 0 << 9;

  {
    Problem<ProblemType::Determinant, double> SolverInstance(MatrixCoeffs);
    SolverInstance.Solve();

    if (Math::Basic::ValueCmp(SolverInstance.Determinant, 18., 1e-4) != 0)
      return false;
  }

  {
    DenseMatrix<double> MatrixInvRef(1e-7, 2, 2);
    Problem<ProblemType::Inverse, double> SolverInstance(MatrixCoeffs);
    SolverInstance.Solve();

    MatrixInvRef << 0.5 << -0.1666666666666
      << 0. << 0.111111111111;

    if (SolverInstance.Inverse != MatrixInvRef)
      return false;
  }

  return true;
}

bool PermuteColsRows()
{
  using Math::Matrix::DenseMatrix;
  using Math::Matrix::IdxPair;
  DenseMatrix<double, 4, 4> MatrixSource(1e-7), MatrixRef(1e-7);

  MatrixSource << 0 << 1 << 2 << 3
               << 4 << 5 << 6 << 7
               << 8 << 9 << 10 << 11
               << 12 << 13 << 14 << 15;

  MatrixRef << 15 << 13 << 0 << 12
            << 7 << 5 << 6 << 4
            << 11 << 9 << 10 << 8
            << 3 << 1 << 2 << 14;

  IdxPair Cell_1(0, 0), Cell_2(3, 2);
  MatrixSource.SwapCells(Cell_1, Cell_2);
  MatrixSource.SwapRows(0, 3);
  MatrixSource.SwapCols(0, 3);

  if (MatrixSource != MatrixRef)
    return false;

  return true;
}

bool ComputeDeterminant55()
{
  using namespace Math::Matrix;

  DenseMatrix<double> MatrixCoeffs(1e-7, 5, 5);

  MatrixCoeffs << 2 << 0 << 4 << 5 << 6
               << 4 << 9 << 8 << 7 << 2
               << 2 << 0 << 2 << 6 << 9
               << 1 << 0 << 6 << 5 << 3
               << 3 << 0 << 7 << 2 << 1;

  {
    Problem<ProblemType::Determinant, double> SolverInstance(MatrixCoeffs);
    SolverInstance.Solve();                                                   // use the default Gauss method

    if (!Math::Basic::ValueCmp(SolverInstance.Determinant, -27., 1e-8) == 0)  // verified (Scilab)
      return false;

    SolverInstance.eStrategy = Problem<ProblemType::Determinant, double>::Det_Strategy::LU_Crout;
    SolverInstance.UploadMatrix(MatrixCoeffs);
    SolverInstance.Solve();

    if (!Math::Basic::ValueCmp(SolverInstance.Determinant, -27., 1e-8) == 0)  // verified (Scilab)
      return false;

    SolverInstance.eStrategy = Problem<ProblemType::Determinant, double>::Det_Strategy::LU_Doolittle;
    SolverInstance.UploadMatrix(MatrixCoeffs);
    SolverInstance.Solve();

    if (!Math::Basic::ValueCmp(SolverInstance.Determinant, -27., 1e-8) == 0)  // verified (Scilab)
      return false;
  }

  return true;
}

bool InvertMatrix55()
{
  using namespace Math::Matrix;
  DenseMatrix<double> MatrixCoeffs(1e-4, 5, 5), RefInv(1e-4, 5, 5);

  MatrixCoeffs << 2 << 0 << 4 << 5 << 6
               << 4 << 9 << 8 << 7 << 2
               << 2 << 0 << 2 << 6 << 9
               << 1 << 0 << 6 << 5 << 3
               << 3 << 0 << 7 << 2 << 1;
  RefInv << 39.666667 << 0. << -22.333333 << -9.6666667 << -8.
         << -19.518519 << 0.1111111 << 11.037037 << 4.5185185 << 4.
         << -23.666667 << -1.776E-15 << 13.333333 << 5.6666667 << 5.
         << 37.666667 << 0. << -21.333333 << -8.6666667 << -8.
         << -28.666667 << 0. << 16.333333 << 6.6666667 << 6.;   // verified (Scilab)

  Problem<ProblemType::Inverse, double> SolverInstance(MatrixCoeffs);
  SolverInstance.Solve();

  return SolverInstance.Inverse == RefInv;
}

bool Serialize()
{
  {
    using Math::Matrix::DenseMatrix;
    DenseMatrix<double> MatrixCoeffs(1e-4, 5, 5);

    MatrixCoeffs << 2 << 0 << 4 << 5 << 6
                 << 4 << 9 << 8 << 7 << 2
                 << 2 << 0 << 2 << 6 << 9
                 << 1 << 0 << 6 << 5 << 3
                 << 3 << 0 << 7 << 2 << 1;

    std::ostringstream strBuffer;
    std::string        strString;

    MatrixCoeffs.Serialize(strBuffer);
    strString = strBuffer.str();

    if (strString != "2 0 4 5 6 4 9 8 7 2 2 0 2 6 9 1 0 6 5 3 3 0 7 2 1 ")
      return false;
  }

  {
    using Math::Matrix::DenseMatrix;
    DenseMatrix<double, 1, 5> MatrixCoeffs(1e-4);

    MatrixCoeffs << 2 << 0 << 4 << 5 << 6;

    std::ostringstream strBuffer;
    std::string        strString;

    MatrixCoeffs.Serialize(strBuffer);
    strString = strBuffer.str();

    if (strString != "2 0 4 5 6 ")
      return false;
  }

  return true;
};

bool LUDecompCrout()
{
  using namespace Math::Matrix;

  DenseMatrix<double> MatrixTerms(1e-4, 5, 5), MatrixTermsPermuted(1e-4, 5, 5), PermMatrix(1e-4, 5, 5);
  DenseMatrix<double> LMatrix(1e-4, 5, 5), UMatrix(1e-4, 5, 5), Res(1e-4, 5, 5);

  MatrixTerms << 2 << 0 << 4 << 5 << 6
              << 4 << 9 << 8 << 7 << 2
              << 2 << 0 << 2 << 6 << 9
              << 1 << 0 << 6 << 5 << 3
              << 3 << 0 << 7 << 2 << 1;

  MatrixTermsPermuted = MatrixTerms;

  Problem<ProblemType::Decomposition, double> SolverInstance(MatrixTerms);
  SolverInstance.eStrategy = Problem<ProblemType::Decomposition, double>::LU_Strategy::Crout;
  SolverInstance.Solve();

  PermMatrix.SetRowPermutationMatrix(SolverInstance.RowPermutation);

  LMatrix.ZeroOut();
  UMatrix.SetIdentity();

  UMatrix.CopyFromMatrixWithPositionCriteria(SolverInstance.LUDecMatrix, Math::Matrix::AboveMainDiagonal());   // copy the elements above the main diagonal (the actual diagonal is not copied)
  LMatrix.CopyFromMatrixWithPositionCriteria(SolverInstance.LUDecMatrix, [](std::size_t i, std::size_t j) -> bool { return i >= j; });  // copy the elements below the main diagonal (the actual diagonal is also copied)

  MatrixTermsPermuted = PermMatrix * MatrixTerms;     // multiply by the permutation matrix to test it out...or we could apply the row permutation directly
  Res = LMatrix * UMatrix;                            // we should obtain the actual matrix...if we use pivoting we will obtain the actual matrix, but pivoted

  if (Res != MatrixTermsPermuted)
    return false;

  return true;
}

bool LUDecompDoolittle()
{
  using namespace Math::Matrix;

  DenseMatrix<double> MatrixTerms(1e-4, 5, 5), MatrixTermsPermuted(1e-4, 5, 5), PermMatrix(1e-4, 5, 5);
  DenseMatrix<double> LMatrix(1e-4, 5, 5), UMatrix(1e-4, 5, 5), Res(1e-4, 5, 5);

  MatrixTerms << 2 << 0 << 4 << 5 << 6
              << 4 << 9 << 8 << 7 << 2
              << 2 << 0 << 2 << 6 << 9
              << 1 << 0 << 6 << 5 << 3
              << 3 << 0 << 7 << 2 << 1;

  MatrixTermsPermuted = MatrixTerms;

  Problem<ProblemType::Decomposition, double> SolverInstance(MatrixTerms);
  SolverInstance.eStrategy = Problem<ProblemType::Decomposition, double>::LU_Strategy::Doolittle;
  SolverInstance.Solve();

  PermMatrix.SetRowPermutationMatrix(SolverInstance.RowPermutation);

  LMatrix.SetIdentity();
  UMatrix.ZeroOut();

  UMatrix.CopyFromMatrixWithPositionCriteria(SolverInstance.LUDecMatrix, [](std::size_t i, std::size_t j) -> bool { return j >= i; });   // copy the elements above the main diagonal (the actual diagonal is also copied)
  LMatrix.CopyFromMatrixWithPositionCriteria(SolverInstance.LUDecMatrix, Math::Matrix::BelowMainDiagonal());    // copy the elements below the main diagonal (the actual diagonal is not copied)

  MatrixTermsPermuted = PermMatrix * MatrixTerms;     // multiply by the permutation matrix to test it out...or we could apply the row permutation directly
  Res = LMatrix * UMatrix;                            // we should obtain the actual matrix...if we use pivoting we will obtain the actual matrix, but pivoted

  if (Res != MatrixTermsPermuted)
    return false;

  return true;
}

bool PermutationTests()
{
  using Math::Misc::Permutation;
  std::map<std::size_t, std::size_t> mapPairs;        // This is actually 1, 4, 3, 2, 0 - (official wiki example)

  mapPairs.insert(std::map<std::size_t, std::size_t>::value_type(0, 1));
  mapPairs.insert(std::map<std::size_t, std::size_t>::value_type(1, 4));
  mapPairs.insert(std::map<std::size_t, std::size_t>::value_type(2, 3));
  mapPairs.insert(std::map<std::size_t, std::size_t>::value_type(3, 2));
  mapPairs.insert(std::map<std::size_t, std::size_t>::value_type(4, 0));

  Permutation<5> CrtPerm(mapPairs);

  std::list<Permutation<5> > listCycles;
  std::list<Permutation<5> > listExpectedCycles;

  {
    std::map<std::size_t, std::size_t> mapExpPairs;   // (0, 1, 4)
    mapExpPairs.insert(std::map<std::size_t, std::size_t>::value_type(0, 1));
    mapExpPairs.insert(std::map<std::size_t, std::size_t>::value_type(1, 4));
    mapExpPairs.insert(std::map<std::size_t, std::size_t>::value_type(2, 2));
    mapExpPairs.insert(std::map<std::size_t, std::size_t>::value_type(3, 3));
    mapExpPairs.insert(std::map<std::size_t, std::size_t>::value_type(4, 0));
    Permutation<5> CrtExpCycle(mapExpPairs);
    listExpectedCycles.push_back(CrtExpCycle);
  }

  {
    std::map<std::size_t, std::size_t> mapExpPairs;   // (2, 3)
    mapExpPairs.insert(std::map<std::size_t, std::size_t>::value_type(0, 0));
    mapExpPairs.insert(std::map<std::size_t, std::size_t>::value_type(1, 1));
    mapExpPairs.insert(std::map<std::size_t, std::size_t>::value_type(2, 3));
    mapExpPairs.insert(std::map<std::size_t, std::size_t>::value_type(3, 2));
    mapExpPairs.insert(std::map<std::size_t, std::size_t>::value_type(4, 4));
    Permutation<5> CrtExpCycle(mapExpPairs);
    listExpectedCycles.push_back(CrtExpCycle);
  }

  CrtPerm.DisjointCycles(listCycles);                 // result should be (1, 4, 3, 2, 0) = (0, 1, 4)*(2, 3)

  std::list<Permutation<5> >::const_iterator itCrt(listCycles.begin()), itEnd(listCycles.end());
  std::list<Permutation<5> >::const_iterator itCrtExp(listExpectedCycles.begin()), itEndExp(listExpectedCycles.end());

  while (itCrt != itEnd && itCrtExp != itEndExp)
  {
    if (*itCrt != *itCrtExp)
      return false;
    else
    {
      ++itCrt;
      ++itCrtExp;
    }
  }

  return true;
}

int main(int argc, char* argv[])
{
  assert(ExtractSubMatrix() && "Failed");
  assert(PerformStaticMatrixOperations() && "Failed");
  assert(PerformDynamicMatrixOperations() && "Failed");
  assert(PerformIteratorOperations() && "Failed");
  PerformConstIteratorOperations();
  PerformReverseIteratorOperations();
  assert(SolveGauss22() && "Failed");
  assert(SolveGauss55() && "Failed");
  ReshapeAndTransposeDynamicMatrix();
  assert(TestEquality() && "Failed");
  assert(InvertMatrix22_1() && "Failed");
  assert(InvertMatrix22_2() && "Failed");
  assert(PermuteColsRows() && "Failed");
  assert(ComputeDeterminant55() && "Failed");
  assert(InvertMatrix55() && "Failed");
  assert(Serialize() && "Failed");
  assert(LUDecompCrout() && "Failed");
  assert(LUDecompDoolittle() && "Failed");
  assert(PermutationTests() && "Failed");

  return 0;
}

