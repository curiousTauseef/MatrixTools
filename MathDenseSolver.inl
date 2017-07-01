template <typename Scalar, std::size_t Rows, std::size_t Cols>
SolverManager<Scalar, Rows, Cols>::SolverManager(DenseMatrix<Scalar, Rows, Cols>& CrtObject, bool bKeepMatrix /* = true */) : m_pCrtMatrix(nullptr), m_nOptions(bKeepMatrix ? KEEP_MATRIX : 0)
{
  EnableOption(PARTIAL_PIVOT, true);
  UploadMatrix(CrtObject, bKeepMatrix);
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
SolverManager<Scalar, Rows, Cols>::~SolverManager()
{
  ReleaseMatrix();
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
bool SolverManager<Scalar, Rows, Cols>::OptionEnabled(unsigned short Flag) const
{
  unsigned short nRes(0);
  
  if (nRes = (m_nOptions & Flag))       
  {
    if (nRes != Flag)                // in case there are multiple flags sent for validation and only some are enabled I will return false 
      return false;
    return true;
  }
  return false;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
inline SolverManager<Scalar, Rows, Cols>& SolverManager<Scalar, Rows, Cols>::EnableOption(unsigned short Flag, bool bEnable = true)
{
  Flag &= ~KEEP_MATRIX;            // make sure that we will set all bits except the Keep matrix flag
  m_nOptions = bEnable ? m_nOptions |= Flag : m_nOptions &= ~Flag;
  return *this;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
inline void SolverManager<Scalar, Rows, Cols>::UploadMatrix(DenseMatrix<Scalar, Rows, Cols>& CrtObject, bool bKeepMatrix /* = true */)
{ 
  ReleaseMatrix();
  bKeepMatrix ? m_nOptions |= KEEP_MATRIX : m_nOptions &= ~KEEP_MATRIX; 
  m_pCrtMatrix = bKeepMatrix ? new DenseMatrix<Scalar, Rows, Cols>(CrtObject) : &CrtObject;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
inline SolverManager<Scalar, Rows, Cols>& SolverManager<Scalar, Rows, Cols>::ReleaseMatrix()
{
  if (OptionEnabled(KEEP_MATRIX) && m_pCrtMatrix)
    delete m_pCrtMatrix;
  return *this;
}




//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
template<typename Scalar, std::size_t Rows, std::size_t Cols>
Problem<ProblemType::Decomposition, Scalar, Rows, Cols>::Problem(DenseMatrix<Scalar, Rows, Cols>& CrtObject, bool bKeepMatrix /* = true */) : eStrategy(LU_Strategy::Crout), SolverManager<Scalar, Rows, Cols>(CrtObject, bKeepMatrix),
                                                                                                                                              RowPermutation(Misc::Permutation<Rows>(m_pCrtMatrix->GetNoRows())),
                                                                                                                                              LUDecMatrix(DenseMatrix<Scalar, Rows, Cols>(m_pCrtMatrix->GetTolerance(), m_pCrtMatrix->GetNoRows(), m_pCrtMatrix->GetNoCols())),
                                                                                                                                              PermutationsPerformed(0)
{
}

template<typename Scalar, std::size_t Rows, std::size_t Cols>
typename Misc::CompileIf<Rows == Cols, bool, void>::RetValue 
Problem<ProblemType::Decomposition, Scalar, Rows, Cols>::Solve()
{
  std::size_t nLSTerms(m_pCrtMatrix->GetLinearSize());
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!nLSTerms)
    throw Basic::StorageException(Basic::StorageException::Type::InvalidDimsForOperation);
  if (m_pCrtMatrix->GetNoCols() != (m_pCrtMatrix->GetNoRows()))
    throw Basic::MatrixException(Basic::MatrixException::Type::NonSquareMatrix);
#endif

  PermutationsPerformed = 0;
  std::size_t nRow = 0, nCol = 0, k = 0, nDiagonalElemIdx = 0, MaxRowIdx = 0, MaxRowPermIdx = 0;
  LUDecMatrix = *m_pCrtMatrix;
  Range PermRange;

  // - Create the LU Decomposition
  for (nDiagonalElemIdx = 0; nDiagonalElemIdx < LUDecMatrix.GetNoRows(); ++nDiagonalElemIdx)
  {
    if (OptionEnabled(PARTIAL_PIVOT))
    {
      PermRange.first  = nDiagonalElemIdx;
      PermRange.second = RowPermutation.GetSize() - 1;

      MaxRowPermIdx = LUDecMatrix.MaxPermutationIndex(RowPermutation, PermRange, nDiagonalElemIdx);   // the permutation position that contains the row of the matrix that has the maximum abs value
      MaxRowIdx     = RowPermutation(MaxRowPermIdx);                                                  // the actual row matrix that has the maximum abs value in the specified range

      if (Basic::ValueCmpZero(LUDecMatrix(MaxRowIdx * LUDecMatrix.GetNoCols() + nDiagonalElemIdx), LUDecMatrix.GetTolerance()) == 0)
        return false;                                        // singular matrix (not an isomorphism)
      else if (nDiagonalElemIdx != MaxRowPermIdx && nDiagonalElemIdx != LUDecMatrix.GetNoRows() - 1)
      {
        RowPermutation.Permute(nDiagonalElemIdx, MaxRowPermIdx);
        ++PermutationsPerformed;
      }
    }

    if (eStrategy == LU_Strategy::Crout)
    {
      for (nCol = nDiagonalElemIdx + 1; nCol < LUDecMatrix.GetNoCols(); ++nCol)
        LUDecMatrix(RowPermutation(nDiagonalElemIdx) * LUDecMatrix.GetNoCols() + nCol) = LUDecMatrix(RowPermutation(nDiagonalElemIdx) * LUDecMatrix.GetNoCols() + nCol) / 
                                                                                         LUDecMatrix(RowPermutation(nDiagonalElemIdx) * LUDecMatrix.GetNoCols() + nDiagonalElemIdx);
    }
    else    // Doolittle
    {
      for (nRow = nDiagonalElemIdx + 1; nRow < LUDecMatrix.GetNoRows(); ++nRow)
        LUDecMatrix(RowPermutation(nRow) * LUDecMatrix.GetNoCols() + nDiagonalElemIdx) = LUDecMatrix(RowPermutation(nRow) * LUDecMatrix.GetNoCols() + nDiagonalElemIdx) / 
                                                                                         LUDecMatrix(RowPermutation(nDiagonalElemIdx) * LUDecMatrix.GetNoCols() + nDiagonalElemIdx);
    }

    for (nRow = nDiagonalElemIdx + 1; nRow < LUDecMatrix.GetNoRows(); ++nRow)
    {
      for (nCol = nDiagonalElemIdx + 1; nCol < LUDecMatrix.GetNoRows(); ++nCol)
        LUDecMatrix(RowPermutation(nRow) * LUDecMatrix.GetNoCols() + nCol) = LUDecMatrix(RowPermutation(nRow) * LUDecMatrix.GetNoCols() + nCol) - 
                                                                             LUDecMatrix(RowPermutation(nRow) * LUDecMatrix.GetNoCols() + nDiagonalElemIdx) * 
                                                                             LUDecMatrix(RowPermutation(nDiagonalElemIdx) * LUDecMatrix.GetNoCols() + nCol);
    }
  }

  if (OptionEnabled(PARTIAL_PIVOT))
    LUDecMatrix.ApplyRowPermutation(RowPermutation);
   
  return true;
}

  


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
template <typename Scalar, std::size_t Rows, std::size_t Cols>
Problem<ProblemType::Determinant, Scalar, Rows, Cols>::Problem(DenseMatrix<Scalar, Rows, Cols>& CrtObject, bool bKeepMatrix /* = true */) : Determinant(0), eStrategy(Det_Strategy::Gauss), 
                                                                                                                                            SolverManager<Scalar, Rows, Cols>(CrtObject, bKeepMatrix),
                                                                                                                                            RowPermutation(Misc::Permutation<Rows>(m_pCrtMatrix->GetNoRows())),
                                                                                                                                            PermutationsPerformed(0)
{
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
typename Misc::CompileIf<Rows == Cols, bool, void>::RetValue 
Problem<ProblemType::Determinant, Scalar, Rows, Cols>::Solve()
{
  std::size_t nLSCoeffs(m_pCrtMatrix->GetLinearSize());
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!nLSCoeffs)
    throw Basic::StorageException(Basic::StorageException::Type::InvalidDimsForOperation);
  if (m_pCrtMatrix->GetNoCols() != m_pCrtMatrix->GetNoRows())
    throw Basic::MatrixException(Basic::MatrixException::Type::NonSquareMatrix);
#endif

  Determinant = 1;
  PermutationsPerformed = 0;

  if (eStrategy != Det_Strategy::Gauss)
  {  
    Problem<ProblemType::Decomposition, Scalar, Rows, Cols> Decomposer(*m_pCrtMatrix, false);      // Don't keep the current matrix that's stored inside the Determinant solver

    switch (eStrategy)
    {
    case Det_Strategy::LU_Crout :
      Decomposer.eStrategy = Problem<ProblemType::Decomposition, Scalar, Rows, Cols>::LU_Strategy::Crout;
      break;
    case Det_Strategy::LU_Doolittle :
      Decomposer.eStrategy = Problem<ProblemType::Decomposition, Scalar, Rows, Cols>::LU_Strategy::Doolittle;
      break;
    }

    Decomposer.EnableOption(PARTIAL_PIVOT, OptionEnabled(PARTIAL_PIVOT));
    if (Decomposer.Solve())
    {
      for (std::size_t i = 0; i < Decomposer.LUDecMatrix.GetNoRows(); ++i)
        Determinant *= Decomposer.LUDecMatrix(i * Decomposer.LUDecMatrix.GetNoCols() + i);
      if (Decomposer.PermutationsPerformed & 1)
        Determinant *= (-1);
      return true;
    }
    return false;
  }

  std::size_t nRow = 0, k = 0, nDiagonalElemIdx = 0, MaxRowIdx = 0, MaxRowPermIdx = 0;
  Scalar Coefficient;
  Range PermRange;

  CrtMatrixType& CrtMatrix(*m_pCrtMatrix);                                                                  // Depending on the option...I'll be converting the initial matrix to reduced echelon form... 

  // - Create the reduced echelon form of the matrix
  for (nDiagonalElemIdx = 0; nDiagonalElemIdx < CrtMatrix.GetNoRows(); ++nDiagonalElemIdx)
  {
    if (OptionEnabled(PARTIAL_PIVOT))
    {
      PermRange.first  = nDiagonalElemIdx;
      PermRange.second = RowPermutation.GetSize() - 1;

      MaxRowPermIdx = CrtMatrix.MaxPermutationIndex(RowPermutation, PermRange, nDiagonalElemIdx);     // the permutation position that contains the row of the matrix that has the maximum abs value
      MaxRowIdx     = RowPermutation(MaxRowPermIdx);                                                  // the actual row matrix that has the maximum abs value in the specified range

      if (Basic::ValueCmpZero(CrtMatrix(MaxRowIdx * CrtMatrix.GetNoCols() + nDiagonalElemIdx), CrtMatrix.GetTolerance()) == 0)
      {
        if (nDiagonalElemIdx == 0)
        {
          Determinant = 0;
          return true;                                                                                 // don't need to determine other stuff
        }
        else
        {
          assert(false && "Failed to compute the determinant. Intermediate null pivot found...to review !");
          return false;
        }
      }
      else if (nDiagonalElemIdx != MaxRowPermIdx && nDiagonalElemIdx != CrtMatrix.GetNoRows() - 1)
      {
        RowPermutation.Permute(nDiagonalElemIdx, MaxRowPermIdx);
        ++PermutationsPerformed;                                                                     // each time we permute 2 rows the determinant gets multiplied by -1                                                         
      }
    }

    for (nRow = nDiagonalElemIdx + 1; nRow < CrtMatrix.GetNoRows(); ++nRow)
    {
      Coefficient = CrtMatrix(RowPermutation(nRow) * CrtMatrix.GetNoCols() + nDiagonalElemIdx) / CrtMatrix(RowPermutation(nDiagonalElemIdx) * CrtMatrix.GetNoCols() + nDiagonalElemIdx);
      for (k = nDiagonalElemIdx + 1; k < CrtMatrix.GetNoCols(); ++k)
        CrtMatrix(RowPermutation(nRow) * CrtMatrix.GetNoCols() + k) = CrtMatrix(RowPermutation(nRow) * CrtMatrix.GetNoCols() + k) - Coefficient * CrtMatrix(RowPermutation(nDiagonalElemIdx) * CrtMatrix.GetNoCols() + k); 
      CrtMatrix(RowPermutation(nRow) * CrtMatrix.GetNoCols() + nDiagonalElemIdx) = 0.;               // zero out the eliminated term
    }
  }

  for (k = 0; k < CrtMatrix.GetNoRows(); ++k)                                                        // fairly easy to compute the determinant since the matrix is diagonal
  {  
    Determinant = Determinant * CrtMatrix(RowPermutation(k) * CrtMatrix.GetNoCols() + k);
    if (Basic::ValueCmpZero(Determinant, CrtMatrix.GetTolerance()) == 0)
      break;
  }

  if (PermutationsPerformed & 1)
    Determinant *= (-1);

  return true;
}




//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
template<typename Scalar, std::size_t Rows, std::size_t Cols>
Problem<ProblemType::Inverse, Scalar, Rows, Cols>::Problem(DenseMatrix<Scalar, Rows, Cols>& CrtObject, bool bKeepMatrix /* = true */) : eStrategy(Inv_Strategy::Gauss_Jordan), 
                                                                                                                                        SolverManager<Scalar, Rows, Cols>(CrtObject, bKeepMatrix),
                                                                                                                                        RowPermutation(Misc::Permutation<Rows>(m_pCrtMatrix->GetNoRows())),
                                                                                                                                        Inverse(DenseMatrix<Scalar, Rows, Cols>(m_pCrtMatrix->GetTolerance(), m_pCrtMatrix->GetNoRows(), m_pCrtMatrix->GetNoCols())),
                                                                                                                                        PermutationsPerformed(0)
{
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
typename Misc::CompileIf<Rows == Cols, bool, void>::RetValue 
Problem<ProblemType::Inverse, Scalar, Rows, Cols>::Solve()
{
  std::size_t nLSCoeffs(m_pCrtMatrix->GetLinearSize());
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!nLSCoeffs)
    throw Basic::StorageException(Basic::StorageException::Type::InvalidDimsForOperation);
  if (m_pCrtMatrix->GetNoCols() != m_pCrtMatrix->GetNoRows())
    throw Basic::MatrixException(Basic::MatrixException::Type::NonSquareMatrix);
#endif

  std::size_t nRow = 0, k = 0, nDiagonalElemIdx = 0, MaxRowIdx = 0, MaxRowPermIdx = 0;
  Scalar Coefficient;
  Range PermRange;

  CrtMatrixType& CrtMatrix(*m_pCrtMatrix);                                                           

  PermutationsPerformed = 0;
  Inverse.SetIdentity();

  // - Create the reduced row echelon form of the matrix
  for (nDiagonalElemIdx = 0; nDiagonalElemIdx < CrtMatrix.GetNoRows(); ++nDiagonalElemIdx)
  {
    if (OptionEnabled(PARTIAL_PIVOT))
    {
      PermRange.first  = nDiagonalElemIdx;
      PermRange.second = RowPermutation.GetSize() - 1;

      MaxRowPermIdx = CrtMatrix.MaxPermutationIndex(RowPermutation, PermRange, nDiagonalElemIdx);     // the permutation position that contains the row of the matrix that has the maximum abs value
      MaxRowIdx     = RowPermutation(MaxRowPermIdx);                                                  // the actual row matrix that has the maximum abs value in the specified range

      if (Basic::ValueCmpZero(CrtMatrix(MaxRowIdx * CrtMatrix.GetNoCols() + nDiagonalElemIdx), CrtMatrix.GetTolerance()) == 0)
        return false;                                        // singular matrix (not an isomorphism)
      else if (nDiagonalElemIdx != MaxRowPermIdx && nDiagonalElemIdx != CrtMatrix.GetNoRows() - 1)
      {
        RowPermutation.Permute(nDiagonalElemIdx, MaxRowPermIdx);
        ++PermutationsPerformed;
      }
    }

    Coefficient = 1 / CrtMatrix(RowPermutation(nDiagonalElemIdx) * CrtMatrix.GetNoCols() + nDiagonalElemIdx);
    CrtMatrix.MultiplyRow(RowPermutation(nDiagonalElemIdx), Coefficient);                          // divide the crt. row by the crt. pivot's value to produce a 1
    Inverse.MultiplyRow(RowPermutation(nDiagonalElemIdx), Coefficient);

    for (nRow = 0; nRow < CrtMatrix.GetNoRows(); ++nRow)
    {
      if (RowPermutation(nRow) != RowPermutation(nDiagonalElemIdx))
      {
        Coefficient = -1 * CrtMatrix(RowPermutation(nRow) * CrtMatrix.GetNoCols() + nDiagonalElemIdx);
        CrtMatrix.AddRowJ_toRowI(RowPermutation(nRow), RowPermutation(nDiagonalElemIdx), Coefficient);
        Inverse.AddRowJ_toRowI(RowPermutation(nRow), RowPermutation(nDiagonalElemIdx), Coefficient);
      }
    }
  }

  Inverse.ApplyRowPermutation(RowPermutation);
  CrtMatrix.ApplyRowPermutation(RowPermutation);

  // Need to check if the crt reduced row echelon form is equal to the identity matrix
  DenseMatrix<Scalar, Rows, Cols> Identity(CrtMatrix.GetTolerance(), CrtMatrix.GetNoRows(), CrtMatrix.GetNoCols());
  Identity.SetIdentity();
  if (CrtMatrix != Identity)
    return false;

  return true;
}




//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
template<typename Scalar, std::size_t Rows, std::size_t Cols>
Problem<ProblemType::LinearEqSolve, Scalar, Rows, Cols>::Problem(DenseMatrix<Scalar, Rows, Cols>& CrtObject, bool bKeepMatrix /* = true */) : eStrategy(Res_Strategy::Gauss), SolverManager<Scalar, Rows, Cols>(CrtObject, bKeepMatrix),
                                                                                                                                              RowPermutation(Misc::Permutation<Rows>(m_pCrtMatrix->GetNoRows())),
                                                                                                                                              Solutions(DenseMatrix<Scalar, Misc::ValueIfRtDynamic<Rows, Misc::DynamicSize, Rows>::Value, Misc::ValueIfRtDynamic<Cols, Misc::DynamicSize, 1>::Value>(m_pCrtMatrix->GetTolerance(), m_pCrtMatrix->GetNoRows(), 1)),
                                                                                                                                              PermutationsPerformed(0)
{
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
typename Misc::CompileIf<(Rows + 1) == Cols || (Rows == Misc::DynamicSize && Cols == Misc::DynamicSize), bool, void>::RetValue 
Problem<ProblemType::LinearEqSolve, Scalar, Rows, Cols>::Solve()
{
  std::size_t nLSCoeffs(m_pCrtMatrix->GetLinearSize());
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!nLSCoeffs)
    throw Basic::StorageException(Basic::StorageException::Type::InvalidDimsForOperation);
  if (m_pCrtMatrix->GetNoCols() != (m_pCrtMatrix->GetNoRows() + 1))
    throw Basic::MatrixException(Basic::MatrixException::Type::NonSquareMatrix);      // non square coefficients matrix
#endif

  PermutationsPerformed = 0;
  std::size_t nRow = 0, k = 0, nDiagonalElemIdx = 0, MaxRowIdx = 0, MaxRowPermIdx = 0;
  Scalar Coefficient, Aux;

  CrtMatrixType& CrtMatrix(*m_pCrtMatrix);  
  Range PermRange;

  // - Create the reduced echelon form of the matrix
  for (nDiagonalElemIdx = 0; nDiagonalElemIdx < CrtMatrix.GetNoRows(); ++nDiagonalElemIdx)
  {
    if (OptionEnabled(PARTIAL_PIVOT))
    {
      PermRange.first  = nDiagonalElemIdx;
      PermRange.second = RowPermutation.GetSize() - 1;

      MaxRowPermIdx = CrtMatrix.MaxPermutationIndex(RowPermutation, PermRange, nDiagonalElemIdx);     // the permutation position that contains the row of the matrix that has the maximum abs value
      MaxRowIdx     = RowPermutation(MaxRowPermIdx);                                                  // the actual row matrix that has the maximum abs value in the specified range

      if (Basic::ValueCmpZero(CrtMatrix(MaxRowIdx * CrtMatrix.GetNoCols() + nDiagonalElemIdx), CrtMatrix.GetTolerance()) == 0)
        return false;
      else if (nDiagonalElemIdx != MaxRowPermIdx && nDiagonalElemIdx != CrtMatrix.GetNoRows() - 1)
      {
        RowPermutation.Permute(nDiagonalElemIdx, MaxRowPermIdx);
        ++PermutationsPerformed;
      }
    }

    for (nRow = nDiagonalElemIdx + 1; nRow < CrtMatrix.GetNoRows(); ++nRow)
    {
      Coefficient = CrtMatrix(RowPermutation(nRow) * CrtMatrix.GetNoCols() + nDiagonalElemIdx) / CrtMatrix(RowPermutation(nDiagonalElemIdx) * CrtMatrix.GetNoCols() + nDiagonalElemIdx);
      for (k = nDiagonalElemIdx + 1; k < CrtMatrix.GetNoCols(); ++k)
        CrtMatrix(RowPermutation(nRow) * CrtMatrix.GetNoCols() + k) = CrtMatrix(RowPermutation(nRow) * CrtMatrix.GetNoCols() + k) - Coefficient * CrtMatrix(RowPermutation(nDiagonalElemIdx) * CrtMatrix.GetNoCols() + k); 
      CrtMatrix(RowPermutation(nRow) * CrtMatrix.GetNoCols() + nDiagonalElemIdx) = 0.;               // zero out the eliminated term
    }
  }

  for (nRow = CrtMatrix.GetNoRows() - 1; nRow >= 0; --nRow)
  {
    Aux = CrtMatrix(RowPermutation(nRow) * CrtMatrix.GetNoCols() + (CrtMatrix.GetNoCols() - 1));
    for (k = CrtMatrix.GetNoRows() - 1; k > nRow; --k)
      Aux = Aux - Solutions(k) * CrtMatrix(RowPermutation(nRow) * CrtMatrix.GetNoCols() + k);
    Solutions(nRow) = Aux / CrtMatrix(RowPermutation(nRow) * CrtMatrix.GetNoCols() + nRow);

    if (nRow == 0)                                                    // just so that I don't plunge into infinity...
      break;
  }

  return true;
}