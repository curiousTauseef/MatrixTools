template <typename Scalar, std::size_t Rows, std::size_t Cols>
DenseMatrix<Scalar, Rows, Cols>::DenseMatrix(const typename DenseMatrix<Scalar, Rows, Cols>::Type& Tolerance, 
                                             std::size_t nRows /* = Rows */, std::size_t nCols /* = Cols*/) : DenseStorage<Scalar, Rows, Cols>(nRows, nCols), 
                                                                                                              m_Tolerance(Tolerance), RefTable(*this)
{

}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
DenseMatrix<Scalar, Rows, Cols>::DenseMatrix(const DenseMatrix<Scalar, Rows, Cols>& RHS) : DenseStorage<Scalar, Rows, Cols>(RHS), 
                                                                                           m_Tolerance(RHS.m_Tolerance), RefTable(*this)
{
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
DenseMatrix<Scalar, Rows, Cols>::DenseMatrix(DenseMatrix<Scalar, Rows, Cols>&& RHS) : DenseStorage<Scalar, Rows, Cols>(std::move(RHS)), 
                                                                                      m_Tolerance(RHS.m_Tolerance), RefTable(*this)
{
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
DenseMatrix<Scalar, Rows, Cols>::~DenseMatrix()
{
  // Dense storage will take care of the actual destructions
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
inline DenseMatrix<Scalar, Rows, Cols>* DenseMatrix<Scalar, Rows, Cols>::MakeCopy() const
{
  return new DenseMatrix<Scalar, Rows, Cols>(*this);
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
DenseMatrix<Scalar, Rows, Cols>& DenseMatrix<Scalar, Rows, Cols>::operator = (const DenseMatrix<Scalar, Rows, Cols>& RHS_Object)
{
  if (this != &RHS_Object)
  {
    DenseStorage<Scalar, Rows, Cols>::operator = (RHS_Object);
    m_Tolerance = RHS_Object.m_Tolerance;
  }
  return *this;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
DenseMatrix<Scalar, Rows, Cols>& DenseMatrix<Scalar, Rows, Cols>::operator = (DenseMatrix<Scalar, Rows, Cols>&& RHS_Object)
{
  if (this != &RHS_Object)
  {
    DenseStorage<Scalar, Rows, Cols>::operator = (std::move(RHS_Object));
    m_Tolerance = RHS_Object.m_Tolerance;
  }
  return *this;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
bool DenseMatrix<Scalar, Rows, Cols>::operator == (const DenseMatrix<Scalar, Rows, Cols>& RHS_Object) const
{
  if ((GetNoRows() != RHS_Object.GetNoRows()) || (GetNoCols() != RHS_Object.GetNoCols()))
    return false;

  Scalar Tolerance = min(m_Tolerance, RHS_Object.m_Tolerance);

  for (std::size_t k = 0; k < GetLinearSize(); ++k)
  {
    if (Basic::ValueCmp(RefTable(k), RHS_Object(k), Tolerance) != 0)
      return false;
  }
    
  return true;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
inline bool DenseMatrix<Scalar, Rows, Cols>::operator != (const DenseMatrix<Scalar, Rows, Cols>& RHS_Object) const
{
  return !(operator == (RHS_Object));
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
const DenseMatrix<Scalar, Rows, Cols> DenseMatrix<Scalar, Rows, Cols>::operator + (const DenseMatrix<Scalar, Rows, Cols>& RHS_Object) const
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (RHS_Object.GetNoCols() != GetNoCols() || RHS_Object.GetNoRows() != GetNoRows()
    || !(RHS_Object.GetNoCols() && GetNoCols() && RHS_Object.GetNoRows() && GetNoRows()))
    throw Basic::StorageException(Basic::StorageException::Type::InvalidDimsForOperation);
#endif

  CrtSpec CrtObject(m_Tolerance, GetNoRows(), GetNoCols());

  std::size_t nSize(GetLinearSize());
  for (std::size_t k = 0; k < nSize; ++k)
    CrtObject(k) = RefTable(k) + RHS_Object(k);

  return CrtObject;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
const DenseMatrix<Scalar, Rows, Cols> DenseMatrix<Scalar, Rows, Cols>::operator - (const DenseMatrix<Scalar, Rows, Cols>& RHS_Object) const
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (RHS_Object.GetNoCols() != GetNoCols() || RHS_Object.GetNoRows() != GetNoRows()
    || !(RHS_Object.GetNoCols() && GetNoCols() && RHS_Object.GetNoRows() && GetNoRows()))
    throw Basic::StorageException(Basic::StorageException::Type::InvalidDimsForOperation);
#endif

  CrtSpec CrtObject(m_Tolerance, GetNoRows(), GetNoCols());

  std::size_t nSize(GetLinearSize());
  for (std::size_t k = 0; k < nSize; ++k)
    CrtObject(k) = RefTable(k) - RHS_Object(k);

  return CrtObject;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
const DenseMatrix<Scalar, Rows, Cols> DenseMatrix<Scalar, Rows, Cols>::operator * (const typename DenseMatrix<Scalar, Rows, Cols>::Type& ScalarValue) const
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!(GetNoCols() && GetNoRows() && m_pDataMassive))
    throw Basic::StorageException(Basic::StorageException::Type::InvalidDimsForOperation);
#endif

  CrtSpec CrtObject(m_Tolerance, GetNoRows(), GetNoCols());

  std::size_t nSize(GetLinearSize());
  for (std::size_t k = 0; k < nSize; ++k)
    CrtObject(k) = RefTable(k) * ScalarValue;

  return CrtObject;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
inline const DenseMatrix<Scalar, Rows, Cols> operator * (const typename DenseMatrix<Scalar, Rows, Cols>::Type& ScalarValue, const DenseMatrix<Scalar, Rows, Cols>& CrtMatrix)
{
  return CrtMatrix * ScalarValue;         // commutative product 
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
template <std::size_t RowsB, std::size_t ColsB>
const DenseMatrix<Scalar, Rows, ColsB> DenseMatrix<Scalar, Rows, Cols>::operator * (const DenseMatrix<Scalar, RowsB, ColsB>& RHS_Object) const
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (RHS_Object.GetNoRows() != GetNoCols() || RHS_Object.GetNoRows() != GetNoRows()
    || !(RHS_Object.GetNoCols() && GetNoCols() && RHS_Object.GetNoRows() && GetNoRows()))
    throw Basic::StorageException(Basic::StorageException::Type::InvalidDimsForOperation);
#endif

  DenseMatrix<Scalar, Rows, ColsB> CrtObject(m_Tolerance, GetNoRows(), RHS_Object.GetNoCols());
  CrtSpec::Type CrtSum = 0;

  PARALLEL_BLOCK
  {
    for (std::size_t i = 0; i < GetNoRows(); ++i)
    {
      for (std::size_t j = 0; j < RHS_Object.GetNoCols(); ++j)
      {
        CrtSum = 0;
        for (std::size_t k = 0; k < RHS_Object.GetNoRows(); ++k)
          CrtSum = CrtSum + (RefTable(i * GetNoCols() + k) * RHS_Object(k * RHS_Object.GetNoCols() + j));

        CrtObject(i * CrtObject.GetNoCols() + j) = CrtSum;
      }
    }
  }

  return CrtObject;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
const DenseMatrix<Scalar, Rows, Cols> DenseMatrix<Scalar, Rows, Cols>::operator - () const
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!(GetNoCols() && GetNoRows() && m_pDataMassive))
    throw Basic::StorageException(Basic::StorageException::Type::InvalidDimsForOperation);
#endif

  CrtSpec CrtObject(m_Tolerance, GetNoRows(), GetNoCols());

  std::size_t nSize(GetLinearSize());
  for (std::size_t k = 0; k < nSize; ++k)
    CrtObject(k) = RefTable(k) * (-1);

  return CrtObject;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
inline typename DenseMatrix<Scalar, Rows, Cols>::iterator 
DenseMatrix<Scalar, Rows, Cols>::begin()
{
  return iterator(GetAddress(0));
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
inline typename DenseMatrix<Scalar, Rows, Cols>::const_iterator 
DenseMatrix<Scalar, Rows, Cols>::cbegin() const
{
  return const_iterator(const_cast<DenseMatrix<Scalar, Rows, Cols>*>(this)->GetAddress(0));
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
inline typename DenseMatrix<Scalar, Rows, Cols>::iterator 
DenseMatrix<Scalar, Rows, Cols>::end()
{
  return iterator(GetAddress(GetLinearSize()));
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
inline typename DenseMatrix<Scalar, Rows, Cols>::const_iterator 
DenseMatrix<Scalar, Rows, Cols>::cend() const
{
  return const_iterator(const_cast<DenseMatrix<Scalar, Rows, Cols>*>(this)->GetAddress(GetLinearSize()));
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
inline typename DenseMatrix<Scalar, Rows, Cols>::reverse_iterator 
DenseMatrix<Scalar, Rows, Cols>::rbegin()
{
  return reverse_iterator(GetAddress(GetLinearSize() - 1));
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
inline typename DenseMatrix<Scalar, Rows, Cols>::const_reverse_iterator 
DenseMatrix<Scalar, Rows, Cols>::crbegin() const
{
  return const_reverse_iterator(const_cast<DenseMatrix<Scalar, Rows, Cols>*>(this)->GetAddress(GetLinearSize() - 1));
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
inline typename DenseMatrix<Scalar, Rows, Cols>::reverse_iterator 
DenseMatrix<Scalar, Rows, Cols>::rend()
{
  return reverse_iterator(GetAddress(0) - 1);
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
inline typename DenseMatrix<Scalar, Rows, Cols>::const_reverse_iterator 
DenseMatrix<Scalar, Rows, Cols>::crend() const
{
  return const_reverse_iterator(const_cast<DenseMatrix<Scalar, Rows, Cols>*>(this)->GetAddress(0) - 1);
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
inline Scalar DenseMatrix<Scalar, Rows, Cols>::GetTolerance() const
{
  return m_Tolerance;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
DenseMatrix<Scalar, Misc::DynamicSize, Misc::DynamicSize> DenseMatrix<Scalar, Rows, Cols>::CreateDynamicMatrix() const
{
  static_assert(Rows != Misc::DynamicSize && Cols != Misc::DynamicSize, "Trying to create a dynamic matrix from a current dynamic object...(should in theory be used for static matrices)");
  return DenseMatrix<Scalar, Misc::DynamicSize, Misc::DynamicSize>(m_Tolerance, GetNoRows(), GetNoCols());
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
DenseMatrix<Scalar, Math::Misc::ValueIfRtDynamic<Rows, Misc::DynamicSize, 1>::Value, Math::Misc::ValueIfRtDynamic<Cols, Misc::DynamicSize, Cols>::Value>
DenseMatrix<Scalar, Rows, Cols>::ExtractRowVector(std::size_t nRowIdx) const
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (nRowIdx >= GetNoRows())
    throw Basic::StorageException(Basic::StorageException::Type::InvalidIndex);
#endif

  DenseMatrix<Scalar, Math::Misc::ValueIfRtDynamic<Rows, Misc::DynamicSize, 1>::Value, Math::Misc::ValueIfRtDynamic<Cols, Misc::DynamicSize, Cols>::Value>
  CrtRow(m_Tolerance, 1, GetNoCols());

  std::size_t nSt(nRowIdx * GetNoCols()), nEnd(nSt + GetNoCols()), j(0);

  for (std::size_t k = nSt; k < nEnd; ++k)
    CrtRow(j++) = RefTable(k);

  return CrtRow;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
DenseMatrix<Scalar, Math::Misc::ValueIfRtDynamic<Rows, Misc::DynamicSize, Rows>::Value, Math::Misc::ValueIfRtDynamic<Cols, Misc::DynamicSize, 1>::Value>
  DenseMatrix<Scalar, Rows, Cols>::ExtractColumnVector(std::size_t nColumnIdx) const
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (nColumnIdx >= GetNoCols())
    throw Basic::StorageException(Basic::StorageException::Type::InvalidIndex);
#endif

  DenseMatrix<Scalar, Math::Misc::ValueIfRtDynamic<Rows, Misc::DynamicSize, Rows>::Value, Math::Misc::ValueIfRtDynamic<Cols, Misc::DynamicSize, 1>::Value>
  CrtCol(m_Tolerance, GetNoRows(), 1);

  std::size_t j(0);

  for (std::size_t k = nColumnIdx; k < GetLinearSize(); k += GetNoCols())
    CrtCol(j++) = RefTable(k);

  return CrtCol;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
template <std::size_t RowsSubMatrix, std::size_t ColsSubMatrix>
void DenseMatrix<Scalar, Rows, Cols>::ExtractSubMatrix(DenseMatrix<Scalar, RowsSubMatrix, ColsSubMatrix>& SubMatrix, const Math::Matrix::IdxPair& DestUprLeftCorner,
                                                       const Math::Matrix::IdxPair& SrcUprLeftCorner, const Math::Matrix::IdxPair& SrcLwrRightCorner) const
{
  std::size_t nLinearSize(GetLinearSize());
  std::size_t nNoRows = SrcLwrRightCorner.first  - SrcUprLeftCorner.first  + 1;           // No. of rows to copy
  std::size_t nNoCols = SrcLwrRightCorner.second - SrcUprLeftCorner.second + 1;           // No. of cols to copy

#ifndef MATH_SAFETY_CHECKS_OFF
  if (!nLinearSize || SrcUprLeftCorner.first >= GetNoRows() || SrcUprLeftCorner.second >= GetNoCols() 
      || SrcLwrRightCorner.first >= GetNoRows() || SrcLwrRightCorner.second >= GetNoCols())
    throw Basic::StorageException(Basic::StorageException::Type::InvalidMemoryAcces);
  if (ACC(SrcUprLeftCorner.first, SrcUprLeftCorner.second, SubMatrix.GetNoCols()) > ACC(SrcLwrRightCorner.first, SrcLwrRightCorner.second, SubMatrix.GetNoCols()) 
      || DestUprLeftCorner.first >= SubMatrix.GetNoRows() || DestUprLeftCorner.second >= SubMatrix.GetNoCols())
    throw Basic::GeneralException(Basic::GeneralException::Type::InvalidParams);
  if (nNoRows > SubMatrix.GetNoRows() || nNoCols > SubMatrix.GetNoCols())
    throw Basic::StorageException(Basic::StorageException::Type::InvalidDimsForOperation);
#endif

  std::size_t iSrc  = SrcUprLeftCorner.first,  iDest = DestUprLeftCorner.first;
  std::size_t jSrc  = SrcUprLeftCorner.second, jDest = DestUprLeftCorner.second;

  while (iSrc <= SrcLwrRightCorner.first && iDest < SubMatrix.GetNoRows())
  {
    jSrc  = SrcUprLeftCorner.second;
    jDest = DestUprLeftCorner.second;
    while (jSrc <= SrcLwrRightCorner.second && jDest < SubMatrix.GetNoCols())
    {
      SubMatrix(ACC(iDest, jDest, SubMatrix.GetNoCols())) = RefTable(ACC(iSrc, jSrc, GetNoCols()));
      ++jSrc, ++jDest;
    }
    ++iSrc, ++iDest;
  }
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
typename Misc::CompileIf<Rows == Cols, bool, void>::RetValue
DenseMatrix<Scalar, Rows, Cols>::IsSymmetric() const
{
  static_assert(Rows == Cols, "The matrix must be square to be able to tell if it is symmetric!");       // for compile time dimensions

#ifndef MATH_SAFETY_CHECKS_OFF
  if (GetNoCols() != GetNoRows())
    throw MatrixException(MatrixException::Type::NonSquareMatrix);
#endif

  for (std::size_t i = 0; j < GetNoRows(); ++j)
  {
    for (std::size_t j = 0; j < GetNoCols(); ++j)
    {
      if (i == j)
        continue;
      if (RefTable(i * GetNoCols() + j) != RefTable(j * GetNoCols() + i))
        return false;
    }
  }

  return true;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
DenseMatrix<Scalar, Cols, Rows> DenseMatrix<Scalar, Rows, Cols>::GetTransposed() const
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!GetLinearSize())
    throw Basic::StorageException(Basic::StorageException::Type::InvalidDimsForOperation);
#endif

  DenseMatrix<Scalar, Cols, Rows> TransposedMat(m_Tolerance, GetNoCols(), GetNoRows());

  for (std::size_t i = 0; i < TransposedMat.GetNoRows(); ++i)
  {
    for (std::size_t j = 0; j < TransposedMat.GetNoCols(); ++j)
      TransposedMat(i * TransposedMat.GetNoCols() + j) = RefTable(j * GetNoCols() + i);
  }

  return TransposedMat;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
template <typename Pred>
DenseMatrix<Scalar, Rows, Cols>& DenseMatrix<Scalar, Rows, Cols>::CopyFromMatrixWithPositionCriteria(DenseMatrix<Scalar, Rows, Cols>& Source, Pred Fctor)
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!GetLinearSize())
    throw Basic::StorageException(Basic::StorageException::Type::InvalidMemoryAcces);
  if (GetNoRows() != Source.GetNoRows() || GetNoCols() != Source.GetNoCols())
    throw Basic::StorageException(Basic::StorageException::Type::InvalidDimsForOperation);
#endif

  if (this != &Source)
  {
    for (std::size_t i = 0; i < GetNoRows(); ++i)
    {
      for (std::size_t j = 0; j < GetNoCols(); ++j)
      {
        if (Fctor(i, j))
          RefTable(i * GetNoCols() + j) = Source(i * GetNoCols() + j);
      }
    }
  }
  return *this;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
typename Misc::CompileIf<Rows == Cols, DenseMatrix<Scalar, Rows, Cols>&, void>::RetValue DenseMatrix<Scalar, Rows, Cols>::SetIdentity()
{
  static_assert(Rows == Cols, "The matrix must be square to be able to set it equal to I !");       // for compile time dimensions

  std::size_t nLinearSize(GetLinearSize());
#ifndef MATH_SAFETY_CHECKS_OFF
  if (GetNoRows() != GetNoCols() || !nLinearSize)
    throw Basic::StorageException(Basic::StorageException::Type::InvalidDimsForOperation);
#endif

  std::size_t k = 0;
  
  ProcessEachElement([=](Type& x) { x = 0; });

  while (k < nLinearSize)
  {
    RefTable(k) = 1;
    k += GetNoCols() + 1;
  }
  return *this;   
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
DenseMatrix<Scalar, Rows, Cols>& DenseMatrix<Scalar, Rows, Cols>::SetRowPermutationMatrix(const Misc::Permutation<Rows>& RowPermutation)
{
  std::size_t nLinearSize(GetLinearSize());
#ifndef MATH_SAFETY_CHECKS_OFF
  if (RowPermutation.GetSize() != GetNoRows())
    throw Basic::GeneralException(Basic::GeneralException::Type::InvalidParams);
#endif

  SetIdentity();
  ApplyRowPermutation(RowPermutation);
  return *this;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
DenseMatrix<Scalar, Rows, Cols>& DenseMatrix<Scalar, Rows, Cols>::SetColumnPermutationMatrix(const Misc::Permutation<Cols>& ColPermutation)
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (RowPermutation.GetSize() != GetNoCols())
    throw Basic::GeneralException(Basic::GeneralException::Type::InvalidParams);
#endif

  SetIdentity();
  ApplyColPermutation(RowPermutation);
  return *this;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
DenseMatrix<Scalar, Rows, Cols>& DenseMatrix<Scalar, Rows, Cols>::ZeroOut()
{
  std::size_t nLinearSize(GetLinearSize());
  if (!nLinearSize)
    return *this;

  std::size_t i = 0;
  while (i < nLinearSize)
    RefTable(i++) = 0;

  return *this;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
DenseMatrix<Scalar, Rows, Cols>& DenseMatrix<Scalar, Rows, Cols>::ZeroOutRow(std::size_t nRow)
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!GetLinearSize() || nRow >= GetNoRows())
    throw Basic::StorageException(Basic::StorageException::Type::InvalidMemoryAcces);
#endif
  
  std::size_t i = nRow * GetNoCols(); // memory block starting point
  std::size_t j = i + GetNoCols();
  for (; i < j; ++i)
    RefTable(i) = 0;

  return *this;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
DenseMatrix<Scalar, Rows, Cols>& DenseMatrix<Scalar, Rows, Cols>::ZeroOutCol(std::size_t nCol)
{
  std::size_t nLinearSize(GetLinearSize());
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!nLinearSize || nCol >= GetNoCols())
    throw Basic::StorageException(Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  std::size_t i = nCol;
  while (i < nLinearSize)
  {
    RefTable(i) = 0;
    i += GetNoCols();
  }

  return *this;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
inline DenseMatrix<Scalar, Rows, Cols>& DenseMatrix<Scalar, Rows, Cols>::MultiplyRow(std::size_t nIdxRow, const typename DenseMatrix<Scalar, Rows, Cols>::Type& Value)
{
  ProcessRow([&Value](CrtSpec::Type& ElemValue){ ElemValue = ElemValue * Value; }, nIdxRow);
  return *this;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
inline DenseMatrix<Scalar, Rows, Cols>& DenseMatrix<Scalar, Rows, Cols>::MultiplyCol(std::size_t nIdxCol, const typename DenseMatrix<Scalar, Rows, Cols>::Type& Value)
{
  ProcessCol([&Value](CrtSpec::Type& ElemValue){ ElemValue = ElemValue * Value; }, nIdxCol);
  return *this;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
DenseMatrix<Scalar, Rows, Cols>& DenseMatrix<Scalar, Rows, Cols>::AddRowJ_toRowI(std::size_t nRowI, std::size_t nRowJ, 
                                                                                 typename DenseMatrix<Scalar, Rows, Cols>::Type Scale)
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!GetLinearSize() || nRowI >= GetNoRows() || nRowJ >= GetNoRows())
    throw Basic::StorageException(Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  std::size_t i(nRowI * GetNoCols()), j(nRowJ * GetNoCols());  // memory block starting points
  std::size_t nEnd = i + GetNoCols();                          // memory block end point

  while (i < nEnd)
  {
    RefTable(i) = RefTable(i) + RefTable(j++) * Scale;
    ++i;
  }

  return *this;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
DenseMatrix<Scalar, Rows, Cols>& DenseMatrix<Scalar, Rows, Cols>::AddColJ_toColI(std::size_t nColI, std::size_t nColJ,
                                                                                 typename DenseMatrix<Scalar, Rows, Cols>::Type Scale)
{
  std::size_t nLinearSize(GetLinearSize());
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!nLinearSize || nColI >= GetNoCols() || nColJ >= GetNoCols())
    throw Basic::StorageException(Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  std::size_t i = nColI;
  std::size_t j = nColJ;

  while (i < nLinearSize && j < nLinearSize)
  {
    RefTable(i) = RefTable(i) + RefTable(j) * Scale;
    i += GetNoCols();
    j += GetNoCols();
  }

  return *this;
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
std::size_t DenseMatrix<Scalar, Rows, Cols>::PivotizeRow(std::size_t nRowPiv, std::size_t nColPiv) const
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!GetLinearSize())
    throw Basic::StorageException(Basic::StorageException::Type::InvalidDimsForOperation);
  if (nRowPiv >= GetNoRows() || nColPiv >= GetNoCols())
    throw Basic::StorageException(Basic::StorageException::Type::InvalidIndex);
#endif

  if (nRowPiv == GetNoRows() - 1)                // Nothing to pivotize
    return (GetNoRows() - 1);

  Type        MaxElem = RefTable(nRowPiv * GetNoCols() + nColPiv);
  std::size_t MaxRowIdx = nRowPiv;

  for (std::size_t nRow = nRowPiv + 1; nRow < RefTable.GetNoRows(); ++nRow)
  {
    if (Basic::Max<Scalar, true>(MaxElem, RefTable(nRow * GetNoCols() + nColPiv)) == RefTable(nRow * GetNoCols() + nColPiv))
    {
      MaxElem   = RefTable(nRow * GetNoCols() + nColPiv);
      MaxRowIdx = nRow;
    }
  }

  return MaxRowIdx;
}

template<typename Scalar, std::size_t Rows, std::size_t Cols>
std::size_t DenseMatrix<Scalar, Rows, Cols>::MaxPermutationIndex(const Misc::Permutation<Math::Misc::ValueIfRtDynamic<Rows, Misc::DynamicSize, Rows>::Value>& Permutation, 
                                                                 const Math::Matrix::Range& PermutationIdxRange, std::size_t nColPiv) const
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!GetLinearSize())
    throw Basic::StorageException(Basic::StorageException::Type::InvalidDimsForOperation);
  if (nColPiv >= GetNoCols())
    throw Basic::StorageException(Basic::StorageException::Type::InvalidIndex);
  if (PermutationIdxRange.first > PermutationIdxRange.second || PermutationIdxRange.second >= Permutation.GetSize() || Permutation.GetSize() != GetNoRows())
    throw Basic::GeneralException(Basic::GeneralException::Type::InvalidParams);
#endif

  Type        MaxElem    = RefTable(Permutation(PermutationIdxRange.first) * GetNoCols() + nColPiv);
  std::size_t MaxPermIdx = PermutationIdxRange.first;

  for (std::size_t i(PermutationIdxRange.first + 1); i <= PermutationIdxRange.second; ++i)
  {
    if (Basic::Max<Scalar, true>(MaxElem, RefTable(Permutation(i) * GetNoCols() + nColPiv)) == RefTable(Permutation(i) * GetNoCols() + nColPiv))
    {
      MaxElem    = RefTable(Permutation(i) * GetNoCols() + nColPiv);
      MaxPermIdx = i;
    }
  }

  return MaxPermIdx;
}

