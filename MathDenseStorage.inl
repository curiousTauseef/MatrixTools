template <typename T, std::size_t Rows, std::size_t Cols>
DS_2DBase<T, Rows, Cols>::DS_2DBase() : m_nCrtIdx(0), m_pDataMassive(nullptr)
{
  // shall not reference any of the defined virtual functions in here...since the vtable hasn't been constructed yet
}

template <typename T, std::size_t Rows, std::size_t Cols>
DS_2DBase<T, Rows, Cols>::DS_2DBase(const DS_2DBase<T, Rows, Cols>& RHS) : m_nCrtIdx(RHS.m_nCrtIdx)
{
  // shall not reference any of the defined virtual functions in here...since the vtable hasn't been constructed yet
}

template <typename T, std::size_t Rows, std::size_t Cols>
DS_2DBase<T, Rows, Cols>::DS_2DBase(DS_2DBase<T, Rows, Cols>&& RHS) : m_nCrtIdx(RHS.m_nCrtIdx), m_pDataMassive(RHS.m_pDataMassive)
{
  RHS.m_pDataMassive = nullptr;
  RHS.m_nCrtIdx = 0;
}

template <typename T, std::size_t Rows, std::size_t Cols>
DS_2DBase<T, Rows, Cols>::~DS_2DBase()
{
}

template <typename T, std::size_t Rows, std::size_t Cols>
inline T* DS_2DBase<T, Rows, Cols>::GetAddress(std::size_t Offset)
{
  return m_pDataMassive + Offset;
}

template <typename T, std::size_t Rows, std::size_t Cols>
DS_2DBase<T, Rows, Cols>& DS_2DBase<T, Rows, Cols>::operator = (const DS_2DBase<T, Rows, Cols>& RHS_Object)
{
  if (this != &RHS_Object)
  {
    std::size_t nSize = GetLinearSize();
    if (!nSize)
    {
      m_nCrtIdx = 0;
      m_pDataMassive = nullptr;
      return *this;
    }

    for (std::size_t k = 0; k < nSize; ++k)             // pseudo destructor call in case
      m_pDataMassive[k].~T();                           // some crazy objects are inserted other than scalars

    m_nCrtIdx = RHS_Object.m_nCrtIdx;

    for (std::size_t k = 0; k < GetLinearSize(); ++k)
      m_pDataMassive[k] = RHS_Object.m_pDataMassive[k];
  }
  return *this;
}

template <typename T, std::size_t Rows, std::size_t Cols>
DS_2DBase<T, Rows, Cols>& DS_2DBase<T, Rows, Cols>::operator = (DS_2DBase<T, Rows, Cols>&& RHS_Object)
{
  if (this != &RHS_Object)
  {
    std::size_t nSize = GetLinearSize();
    if (!nSize)
    {
      m_nCrtIdx = 0;
      m_pDataMassive = nullptr;
      return *this;
    }

    for (std::size_t k = 0; k < nSize; ++k)      // pseudo destructor call in case
      m_pDataMassive[k].~T();                           // some crazy objects are inserted other than scalars

    m_nCrtIdx = RHS_Object.m_nCrtIdx;
    m_pDataMassive = RHS_Object.m_pDataMassive;

    RHS_Object.m_pDataMassive = nullptr;                // prevent destructor actions on this memory address
    RHS_Object.m_nCrtIdx = 0;
  }
  return *this;
}

template <typename T, std::size_t Rows, std::size_t Cols>
inline typename DS_2DBase<T, Rows, Cols>::Type* DS_2DBase<T, Rows, Cols>::operator [] (std::size_t nIdx) 
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (m_pDataMassive == nullptr || nIdx >= GetNoRows())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  return m_pDataMassive + GetNoCols() * nIdx;
}

template <typename T, std::size_t Rows, std::size_t Cols>
inline const typename DS_2DBase<T, Rows, Cols>::Type* DS_2DBase<T, Rows, Cols>::operator [] (std::size_t nIdx) const
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (m_pDataMassive == nullptr || nIdx >= GetNoRows())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  return m_pDataMassive + GetNoCols() * nIdx;
}

template <typename T, std::size_t Rows, std::size_t Cols>
inline typename DS_2DBase<T, Rows, Cols>::Type& DS_2DBase<T, Rows, Cols>::operator () (std::size_t nIdxRow, std::size_t nIdxCol)
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (m_pDataMassive == nullptr || nIdxRow >= GetNoRows() || nIdxCol >= GetNoCols())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  return *(m_pDataMassive + GetNoCols() * nIdxRow + nIdxCol);
}

template <typename T, std::size_t Rows, std::size_t Cols>
inline const typename DS_2DBase<T, Rows, Cols>::Type& DS_2DBase<T, Rows, Cols>::operator () (std::size_t nIdxRow, std::size_t nIdxCol) const
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (m_pDataMassive == nullptr || nIdxRow >= GetNoRows() || nIdxCol >= GetNoCols())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  return *(m_pDataMassive + GetNoCols() * nIdxRow + nIdxCol);
}

template <typename T, std::size_t Rows, std::size_t Cols>
inline typename DS_2DBase<T, Rows, Cols>::Type& DS_2DBase<T, Rows, Cols>::operator () (std::size_t nIdx)
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (m_pDataMassive == nullptr || nIdx >= GetLinearSize())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  return *(m_pDataMassive + nIdx);
}

template <typename T, std::size_t Rows, std::size_t Cols>
inline const typename DS_2DBase<T, Rows, Cols>::Type& DS_2DBase<T, Rows, Cols>::operator () (std::size_t nIdx) const
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (m_pDataMassive == nullptr || nIdx >= GetLinearSize())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  return *(m_pDataMassive + nIdx);
}

template <typename T, std::size_t Rows, std::size_t Cols>
DS_2DBase<T, Rows, Cols>& DS_2DBase<T, Rows, Cols>::operator << (const typename DS_2DBase<T, Rows, Cols>::Type& Value)
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (m_pDataMassive == nullptr)
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  m_pDataMassive[m_nCrtIdx] = Value;
  m_nCrtIdx = (m_nCrtIdx + 1) % GetLinearSize();

  return *this;
}

template <typename T, std::size_t Rows, std::size_t Cols>
std::ostream& operator << (std::ostream& Stream, const DS_2DBase<T, Rows, Cols>& CrtObject)
{
  std::size_t nRows(CrtObject.GetNoRows()), nCols(CrtObject.GetNoCols());

  for (std::size_t i = 0; i < nRows; ++i)
  {
    for (std::size_t j = 0; j < nCols; ++j)
      Stream << CrtObject(i,j);
  }

  return Stream;
}

template <typename T, std::size_t Rows, std::size_t Cols>
inline void DS_2DBase<T, Rows, Cols>::SwapContents(const DS_2DBase<T, Rows, Cols>& Source)
{
  if (this != &Source)
  {
    std::swap(m_pDataMassive, Source.m_pDataMassive);  // exchanging memory addresses
    std::swap(m_nCrtIdx, Source.Source);
  }
}

template <typename T, std::size_t Rows, std::size_t Cols>
void DS_2DBase<T, Rows, Cols>::CopyOnRow(std::size_t nRow, typename DS_2DBase<T, Rows, Cols>::Type Array[], 
                                         std::size_t nArraySize)
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!GetLinearSize() || nRow >= GetNoRows() || nArraySize <= 0 || !Array)
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  std::size_t i = nRow * GetNoCols();    // memory block starting point
  std::size_t j = 0;
  for (;i < GetNoCols() && j < nArraySize; ++i, ++j)
    m_pDataMassive[i] = Array[j];
}

template <typename T, std::size_t Rows, std::size_t Cols>
void DS_2DBase<T, Rows, Cols>::CopyOnRow(const DS_2DBase<T, Misc::ValueIfRtDynamic<Cols, Misc::DynamicSize, 1>::Value, Misc::ValueIfRtDynamic<Cols, Misc::DynamicSize, Cols>::Value>& RowLine, std::size_t nRow) 
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!GetLinearSize() || nRow >= GetNoRows() || !RowLine.GetLinearSize())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  for (std::size_t k = 0; k < GetNoCols(); ++k)
    (*this)(nRow * GetNoCols() + k) = RowLine(k); 
}

template <typename T, std::size_t Rows, std::size_t Cols>
void DS_2DBase<T, Rows, Cols>::CopyOnCol(std::size_t nRow, typename DS_2DBase<T, Rows, Cols>::Type Array[], 
                                         std::size_t nArraySize)
{
  std::size_t nLinearSize(GetLinearSize());
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!nLinearSize || nCol >= GetNoCols() || nArraySize <= 0 || !Array)
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  std::size_t i = 0, j = 0;
  while (i < nLinearSize && j < nArraySize)
  {
    m_pDataMassive[i] = Array[j++];
    i += GetNoCols();
  }
}

template <typename T, std::size_t Rows, std::size_t Cols>
void DS_2DBase<T, Rows, Cols>::CopyOnCol(const DS_2DBase<T, Misc::ValueIfRtDynamic<Rows, Misc::DynamicSize, Rows>::Value, Misc::ValueIfRtDynamic<Rows, Misc::DynamicSize, 1>::Value>& ColLine, std::size_t nCol)
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!GetLinearSize() || nCol >= GetNoCols() || !ColLine.GetLinearSize())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  for (std::size_t k = 0; k < GetNoRows(); ++k)
    (*this)(k * GetNoCols() + nCol) = ColLine(k); 
}

template <typename T, std::size_t Rows, std::size_t Cols>
void DS_2DBase<T, Rows, Cols>::CopyFromMassive(typename DS_2DBase<T, Rows, Cols>::Type Array[], std::size_t nArraySize)
{
  std::size_t nBlockSz = GetLinearSize();
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!nBlockSz || nArraySize <= 0 || !Array)
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  std::size_t i = 0, j = 0;
  for (;i < nBlockSz && j < nArraySize; ++i, ++j)
    m_pDataMassive[i] = Array[j];
}

template <typename T, std::size_t Rows, std::size_t Cols>
void DS_2DBase<T, Rows, Cols>::SwapRows(std::size_t nRowI, std::size_t nRowJ)
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!GetLinearSize() || nRowI >= GetNoRows() || nRowJ >= GetNoRows())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidIndex);
#endif

  if (nRowI == nRowJ)
    return;

  std::size_t i(nRowI * GetNoCols()), j(nRowJ * GetNoCols());
  std::size_t nEnd = i + GetNoCols();

  while (i < nEnd)
    std::swap(m_pDataMassive[i++], m_pDataMassive[j++]);
}

template <typename T, std::size_t Rows, std::size_t Cols>
void DS_2DBase<T, Rows, Cols>::SwapCols(std::size_t nColI, std::size_t nColJ)
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!GetLinearSize() || nColI >= GetNoCols() || nColI >= GetNoCols())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidIndex);
#endif

  if (nColI == nColJ)
    return;

  std::size_t i = 0;
  std::size_t k1 = nColI;
  std::size_t k2 = nColJ;

  while (i < GetNoRows())
  {
    std::swap(m_pDataMassive[k1], m_pDataMassive[k2]);
    k1 += GetNoCols();
    k2 += GetNoCols();
    ++i;
  }
}

template <typename T, std::size_t Rows, std::size_t Cols>
inline void DS_2DBase<T, Rows, Cols>::SwapCells(const typename Math::Matrix::IdxPair& Cell_1, const Math::Matrix::IdxPair& Cell_2)
{
  std::swap(operator()(Cell_1.first, Cell_1.second), operator()(Cell_2.first, Cell_2.second));
}

template <typename T, std::size_t Rows, std::size_t Cols>
void DS_2DBase<T, Rows, Cols>::ApplyRowPermutation(const Math::Misc::Permutation<Rows>& RowPermutation)
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!GetLinearSize())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
  if (RowPermutation.GetSize() != GetNoRows())
    throw Math::Basic::GeneralException(Math::Basic::GeneralException::Type::InvalidParams);
#endif

  Misc::DynamicPermutation PermSwaps(GetNoRows());
  std::size_t i = 0;

  while (PermSwaps != RowPermutation && i < RowPermutation.GetSize())
  {
    SwapRows(i, RowPermutation(i));
    PermSwaps.Permute(i, RowPermutation(i));
    ++i;
  }
}

template <typename T, std::size_t Rows, std::size_t Cols>
void DS_2DBase<T, Rows, Cols>::ApplyColPermutation(const Math::Misc::Permutation<Cols>& ColPermutation)
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!GetLinearSize())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
  if (ColPermutation.GetSize() != GetNoCols())
    throw Math::Basic::GeneralException(Math::Basic::GeneralException::Type::InvalidParams);
#endif

  Misc::Permutation<Cols> PermSwaps;
  std::size_t i = 0;

  while (PermSwaps != ColPermutation && i < ColPermutation.GetSize())
  {
    SwapCols(i, ColPermutation(i));
    PermSwaps.Permute(i, ColPermutation(i));
    ++i;
  }
}

template <typename T, std::size_t Rows, std::size_t Cols>
template <typename F>
void DS_2DBase<T, Rows, Cols>::ProcessEachElement(F Function)
{
  std::size_t nLinearSize(GetLinearSize());
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!nLinearSize)
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  std::size_t i = 0;
  while (i < nLinearSize)
  {
    Function(m_pDataMassive[i]);
    ++i;
  }
}

template <typename T, std::size_t Rows, std::size_t Cols>
template <typename F>
void DS_2DBase<T, Rows, Cols>::ProcessRow(F Function, std::size_t nIdxRow)
{
  std::size_t nLinearSize(GetLinearSize());
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!nLinearSize || nIdxRow >= GetNoRows())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  std::size_t i = nIdxRow * GetNoCols();
  std::size_t j = i + GetNoCols();
  for (; i < j; ++i)
    Function(m_pDataMassive[i]);
}

template <typename T, std::size_t Rows, std::size_t Cols>
template <typename F>
void DS_2DBase<T, Rows, Cols>::ProcessCol(F Function, std::size_t nIdxCol)
{
  std::size_t nLinearSize(GetLinearSize());
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!nLinearSize || nIdxCol >= GetNoCols())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  std::size_t i = nIdxCol;
  while (i < nLinearSize)
  {
    Function(m_pDataMassive[i]);
    i += GetNoCols();
  }
}

template <typename T, std::size_t Rows, std::size_t Cols>
std::ostream& DS_2DBase<T, Rows, Cols>::Serialize(std::ostream& Stream) const
{
  std::size_t nSize(GetLinearSize());
  for (std::size_t k = 0; k < nSize; ++k)
    Stream << m_pDataMassive[k] << " ";

  return Stream;
}

template <typename T, std::size_t Rows, std::size_t Cols>
void DS_2DBase<T, Rows, Cols>::CopyToMassive(typename DS_2DBase<T, Rows, Cols>::Type*& Array, std::size_t& nSize) const
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (Array)
    delete [] Array;
#endif

  nSize = GetLinearSize();
  Array = new Type[nSize];

  for (std::size_t k = 0; k < nSize; ++k)
    *(Array + k) = *(m_pDataMassive + k);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename T, std::size_t Rows, std::size_t Cols>
DenseStorage<T, Rows, Cols>::DenseStorage(std::size_t nRows /* = Misc::DynamicSize */, std::size_t nCols /* = Misc::DynamicSize */) : DS_2DBase<T, Rows, Cols>()
{
  static_assert(Rows > 0 && Cols > 0, "The rows and columns must both be positive numbers");
  m_pDataMassive = new Type[Rows * Cols];
}

template <typename T, std::size_t Rows, std::size_t Cols>
DenseStorage<T, Rows, Cols>::DenseStorage(const DenseStorage<T, Rows, Cols>& RHS) : DS_2DBase<T, Rows, Cols>(RHS)
{
  if (RHS.GetLinearSize())
  {
    m_pDataMassive = new Type[RHS.GetLinearSize()];

    for (std::size_t k = 0; k < GetLinearSize(); ++k)
      m_pDataMassive[k] = RHS.m_pDataMassive[k];
  }
}

template <typename T, std::size_t Rows, std::size_t Cols>
DenseStorage<T, Rows, Cols>::DenseStorage(DenseStorage<T, Rows, Cols>&& RHS) : DS_2DBase<T, Rows, Cols>(static_cast<DS_2DBase<T, Rows, Cols>&&>(RHS))
{
}

template <typename T, std::size_t Rows, std::size_t Cols>
DenseStorage<T, Rows, Cols>::~DenseStorage()
{
  if (m_pDataMassive)
  {
    std::size_t nSize = GetLinearSize();
    for (std::size_t k = 0; k < nSize; ++k)      // pseudo destructor call in case
      m_pDataMassive[k].~T();                    // some crazy objects are inserted other than scalars
    delete[] m_pDataMassive;                     // the user is responsible for clearing up heap memory
    m_pDataMassive = nullptr;
  }
}

template <typename T, std::size_t Rows, std::size_t Cols>
inline DenseStorage<T, Rows, Cols>* DenseStorage<T, Rows, Cols>::MakeCopy() const
{
  return new DenseStorage<T, Rows, Cols>(*this);
}

template <typename T, std::size_t Rows, std::size_t Cols>
DenseStorage<T, Rows, Cols>& DenseStorage<T, Rows, Cols>::operator = (const DenseStorage<T, Rows, Cols>& RHS_Object)
{
  DS_2DBase<T, Rows, Cols>::operator = (RHS_Object);
  return *this;
}

template <typename T, std::size_t Rows, std::size_t Cols>
DenseStorage<T, Rows, Cols>& DenseStorage<T, Rows, Cols>::operator = (DenseStorage<T, Rows, Cols>&& RHS_Object)
{
  DS_2DBase<T, Rows, Cols>::operator = (static_cast<DS_2DBase<T, Rows, Cols>&&>(RHS_Object));
  return *this;
}

template <typename T, std::size_t Rows, std::size_t Cols>
std::size_t DenseStorage<T, Rows, Cols>::GetNoRows() const
{
  return Rows;
}

template <typename T, std::size_t Rows, std::size_t Cols>
std::size_t DenseStorage<T, Rows, Cols>::GetNoCols() const
{
  return Cols;
}

template <typename T, std::size_t Rows, std::size_t Cols>
std::size_t DenseStorage<T, Rows, Cols>::GetLinearSize() const
{
  return const_cast<DenseStorage<T, Rows, Cols>*>(this)->m_pDataMassive ? GetNoCols() * GetNoRows() : 0;
}

template <typename T, std::size_t Rows, std::size_t Cols>
inline DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize> DenseStorage<T, Rows, Cols>::CreateRTDynamicStorage() const
{
  return DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>(GetNoRows(), GetNoCols());
}


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>::DenseStorage(std::size_t nSize) : DS_2DBase<T, Misc::DynamicSize, Misc::DynamicSize>(), m_nNoCols(nSize), m_nNoRows(nSize)
{
  if (nSize > 0)
    m_pDataMassive = new Type[m_nNoRows * m_nNoCols];
  else
    assert(false && L"The size cannot be null !");
}

template<typename T>
DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>::DenseStorage(std::size_t nRows, std::size_t nCols) : DS_2DBase<T, Misc::DynamicSize, Misc::DynamicSize>(), m_nNoCols(nCols), m_nNoRows(nRows)
{
  if (nRows > 0 && nCols > 0)
    m_pDataMassive = new Type[m_nNoRows * m_nNoCols];
  else
    assert(false && L"The sizes cannot be null !");
}

template<typename T>
DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>::DenseStorage(const DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>& RHS) : DS_2DBase<T, Misc::DynamicSize, Misc::DynamicSize>(RHS), 
                                                                                                                                        m_nNoCols(RHS.m_nNoCols), m_nNoRows(RHS.m_nNoRows)
{
  if (RHS.GetLinearSize())
  {
    m_pDataMassive = new Type[RHS.GetLinearSize()];

    for (std::size_t k = 0; k < GetLinearSize(); ++k)
      m_pDataMassive[k] = RHS.m_pDataMassive[k];
  }
}

template<typename T>
DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>::DenseStorage(DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>&& RHS) : DS_2DBase<T, Misc::DynamicSize, Misc::DynamicSize>(static_cast<DS_2DBase<T, Misc::DynamicSize, Misc::DynamicSize>&&>(RHS)), 
                                                                                                                                   m_nNoCols(RHS.m_nNoCols), m_nNoRows(RHS.m_nNoRows)
{
  RHS.m_nNoRows = RHS.m_nNoCols = 0;
}

template<typename T>
DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>::~DenseStorage()
{
  if (m_pDataMassive)
  {
    std::size_t nSize = GetLinearSize();
    for (std::size_t k = 0; k < nSize; ++k)      // pseudo destructor call in case
      m_pDataMassive[k].~T();                    // some crazy objects are inserted other than scalars
    delete[] m_pDataMassive;                     // the user is responsible for clearing up heap memory
    m_pDataMassive = nullptr;
  }
}

template<typename T>
inline DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>* DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>::MakeCopy() const
{
  return new DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>(*this);
}

template<typename T>
DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>& DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>::operator = (const DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>& RHS_Object)
{
  DS_2DBase<T, Misc::DynamicSize, Misc::DynamicSize>::operator = (RHS_Object);
  m_nNoCols = RHS_Object.m_nNoCols;
  m_nNoRows = RHS_Object.m_nNoRows;
  return *this;
}

template<typename T>
DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>& DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>::operator = (DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>&& RHS_Object)
{
  DS_2DBase<T, Misc::DynamicSize, Misc::DynamicSize>::operator = (static_cast<DS_2DBase<T, Misc::DynamicSize, Misc::DynamicSize>&&>(RHS_Object));
  m_nNoCols = RHS_Object.m_nNoCols;
  m_nNoRows = RHS_Object.m_nNoRows;
  return *this;
}

template <typename T>
void DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>::Reshape(std::size_t nRowsNo, std::size_t nColsNo)
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!(nRowsNo && nColsNo))
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidDimsForOperation);
#endif

  if (GetNoRows() == nRowsNo && GetNoCols() == nColsNo)
    return;

  DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize> NewTable(nRowsNo, nColsNo);

  std::size_t i(0), j(0);
  std::size_t nRowsConsidered = min(GetNoRows(), nRowsNo);
  std::size_t nColsConsidered = min(GetNoCols(), nColsNo);

  for (std::size_t i = 0; i < nRowsConsidered; ++i)
  {
    for (std::size_t j = 0; j < nColsConsidered; ++j)
    {
      NewTable(i * NewTable.GetNoCols() + j) = (*this)(i * GetNoCols() + j);  // to improve if the crt type T is actually a custom scalar that holds a large amount of memory-> move semantics could be used
      (*this)(i * GetNoCols() + j).~T();                                      // for now copying and destroying the type T should do
    }
  }

  (*this) = static_cast<DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>&&>(NewTable);
}

template<typename T>
void DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>::AppendRow(const DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>& Row)
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (Row.GetNoRows() > 1 || Row.GetNoCols() != GetNoCols() || !Row.GetLinearSize())
    throw Math::Basic::GeneralException(Math::Basic::GeneralException::Type::InvalidParams);
#endif

  std::size_t nPrevSize = GetLinearSize();
  Reshape(GetNoRows() + 1, GetNoCols());
  std::size_t nNewSize = GetLinearSize(), j = 0;

  for (std::size_t k = nPrevSize; k < nNewSize; ++k)
    *GetAddress(k) = Row(j++);
}

template<typename T>
template <std::size_t Cols>
void DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>::AppendRow(const DenseStorage<T, 1, Cols>& Row)
{
  static_assert(Cols > 0, "The number of column must be > 0");
  std::size_t nPrevSize = GetLinearSize();
  Reshape(GetNoRows() + 1, GetNoCols());
  std::size_t nNewSize = GetLinearSize(), j = 0;

  for (std::size_t k = nPrevSize; k < nNewSize; ++k)
    *GetAddress(k) = Row(j++);
}

template<typename T>
void DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>::AppendCol(const DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>& Col)
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (Col.GetNoCols() > 1 || Col.GetNoRows() != GetNoRows() || !Col.GetLinearSize())
    throw Math::Basic::GeneralException(Math::Basic::GeneralException::Type::InvalidParams);
#endif

  Reshape(GetNoRows(), GetNoCols() + 1);

  std::size_t j = Col.GetLinearSize() - 1;
  std::size_t k = GetLinearSize() - 1;

  while (k > 0)
  {
    *GetAddress(k) = Col(j--);

    if (k < GetNoCols())
      break;
    else
      k -= GetNoCols();
  }
}

template<typename T>
template <std::size_t Rows>
void DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>::AppendCol(const DenseStorage<T, Rows, 1>& Col)
{
  static_assert(Rows > 0, "The number of rows must be > 0");
  Reshape(GetNoRows(), GetNoCols() + 1);

  std::size_t j = Col.GetLinearSize() - 1;
  std::size_t k = GetLinearSize() - 1;

  while (k > 0)
  {
    m_pDataMassive[k] = Col(j--);

    if (k < GetNoCols())
      break;
    else
      k -= GetNoCols();
  }
}

template<typename T>
std::size_t DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>::GetNoRows() const
{
  return m_nNoRows;
}

template<typename T>
std::size_t DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>::GetNoCols() const
{
  return m_nNoCols;
}

template<typename T>
std::size_t DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>::GetLinearSize() const
{
  return const_cast<DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>*>(this)->m_pDataMassive ? GetNoCols() * GetNoRows() : 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename T, std::size_t Cols>
DenseStorage<T, 1, Cols>::DenseStorage(std::size_t nRows /* = 1 */, std::size_t nCols /* = Cols */) : DS_2DBase<T, 1, Cols>()
{
  static_assert(Cols > 0, "The number of columns must be positive");
  m_pDataMassive = new Type[Cols];
}

template <typename T, std::size_t Cols>
DenseStorage<T, 1, Cols>::DenseStorage(const DenseStorage<T, 1, Cols>& RHS) : DS_2DBase<T, 1, Cols>(RHS)
{
  if (RHS.GetLinearSize())
  {
    m_pDataMassive = new Type[RHS.GetLinearSize()];

    for (std::size_t k = 0; k < GetLinearSize(); ++k)
      m_pDataMassive[k] = RHS.m_pDataMassive[k];
  }
}

template <typename T, std::size_t Cols>
DenseStorage<T, 1, Cols>::DenseStorage(DenseStorage<T, 1, Cols>&& RHS) : DS_2DBase<T, 1, Cols>(static_cast<DS_2DBase<T, 1, Cols>&&>(RHS))
{
}

template <typename T, std::size_t Cols>
DenseStorage<T, 1, Cols>::~DenseStorage()
{
  if (m_pDataMassive)
  {
    std::size_t nSize = GetLinearSize();
    for (std::size_t k = 0; k < nSize; ++k)      // pseudo destructor call in case
      m_pDataMassive[k].~T();                    // some crazy objects are inserted other than scalars
    delete[] m_pDataMassive;                     // the user is responsible for clearing up heap memory
    m_pDataMassive = nullptr;
  }
}

template <typename T, std::size_t Cols>
inline DenseStorage<T, 1, Cols>* DenseStorage<T, 1, Cols>::MakeCopy() const
{
  return new DenseStorage<T, 1, Cols>(*this);
}

template <typename T, std::size_t Cols>
DenseStorage<T, 1, Cols>& DenseStorage<T, 1, Cols>::operator = (const DenseStorage<T, 1, Cols>& RHS_Object)
{
  DS_2DBase<T, 1, Cols>::operator = (RHS_Object);
  return *this;
}

template <typename T, std::size_t Cols>
DenseStorage<T, 1, Cols>& DenseStorage<T, 1, Cols>::operator = (DenseStorage<T, 1, Cols>&& RHS_Object)
{
  DS_2DBase<T, 1, Cols>::operator = (static_cast<DS_2DBase<T, 1, Cols>&&>(RHS_Object));
  return *this;
}

template <typename T, std::size_t Cols>
inline typename T& DenseStorage<T, 1, Cols>::operator [] (std::size_t nIdx) 
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (m_pDataMassive == nullptr || nIdx >= GetNoCols())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  return *(m_pDataMassive + nIdx);
}

template <typename T, std::size_t Cols>
inline const typename T& DenseStorage<T, 1, Cols>::operator [] (std::size_t nIdx) const
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (m_pDataMassive == nullptr || nIdx >= GetNoCols())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  return *(m_pDataMassive + nIdx);
}

template <typename T, std::size_t Cols>
inline typename T& DenseStorage<T, 1, Cols>::operator () (std::size_t nIdxCol)
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (m_pDataMassive == nullptr || nIdxCol >= GetNoCols())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  return *(m_pDataMassive + nIdxCol);
}

template <typename T, std::size_t Cols>
inline const typename T& DenseStorage<T, 1, Cols>::operator () (std::size_t nIdxCol) const
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (m_pDataMassive == nullptr || nIdxCol >= GetNoCols())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  return *(m_pDataMassive + nIdxCol);
}

template <typename T, std::size_t Cols>
DenseStorage<T, Cols, 1> DenseStorage<T, 1, Cols>::GetColumnVector()
{
  DenseStorage<T, Cols, 1> Transposed;
  for (std::size_t k = 0; k < GetLinearSize(); ++k)
    Transposed[k] = (*this)[k];
  return Transposed;
}

template <typename T, std::size_t Cols>
std::size_t DenseStorage<T, 1, Cols>::GetNoRows() const
{
  return 1;
}

template <typename T, std::size_t Cols>
std::size_t DenseStorage<T, 1, Cols>::GetNoCols() const
{
  return Cols;
}

template<typename T, std::size_t Cols>
std::size_t DenseStorage<T, 1, Cols>::GetLinearSize() const
{
  return const_cast<DenseStorage<T, 1, Cols>*>(this)->m_pDataMassive ? GetNoCols() * GetNoRows() : 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename T, std::size_t Rows>
DenseStorage<T, Rows, 1>::DenseStorage(std::size_t nRows /* = Rows */, std::size_t nCols /* = 1 */) : DS_2DBase<T, Rows, 1>()
{
  static_assert(Rows > 0, "The number of rows must be positive");
  m_pDataMassive = new Type[Rows];
}

template <typename T, std::size_t Rows>
DenseStorage<T, Rows, 1>::DenseStorage(const DenseStorage<T, Rows, 1>& RHS) : DS_2DBase<T, Rows, 1>(RHS)
{
  if (RHS.GetLinearSize())
  {
    m_pDataMassive = new Type[RHS.GetLinearSize()];

    for (std::size_t k = 0; k < GetLinearSize(); ++k)
      m_pDataMassive[k] = RHS.m_pDataMassive[k];
  }
}

template <typename T, std::size_t Rows>
DenseStorage<T, Rows, 1>::DenseStorage(DenseStorage<T, Rows, 1>&& RHS) : DS_2DBase<T, Rows, 1>(static_cast<DS_2DBase<T, Rows, 1>&&>(RHS))
{
}

template <typename T, std::size_t Rows>
DenseStorage<T, Rows, 1>::~DenseStorage()
{
  if (m_pDataMassive)
  {
    std::size_t nSize = GetLinearSize();
    for (std::size_t k = 0; k < nSize; ++k)      // pseudo destructor call in case
      m_pDataMassive[k].~T();                    // some crazy objects are inserted other than scalars
    delete[] m_pDataMassive;                     // the user is responsible for clearing up heap memory
    m_pDataMassive = nullptr;
  }
}

template <typename T, std::size_t Rows>
inline DenseStorage<T, Rows, 1>* DenseStorage<T, Rows, 1>::MakeCopy() const
{
  return new DenseStorage<T, Rows, 1>(*this);
}

template <typename T, std::size_t Rows>
DenseStorage<T, Rows, 1>& DenseStorage<T, Rows, 1>::operator = (const DenseStorage<T, Rows, 1>& RHS_Object)
{
  DS_2DBase<T, Rows, 1>::operator = (RHS_Object);
  return *this;
}

template <typename T, std::size_t Rows>
DenseStorage<T, Rows, 1>& DenseStorage<T, Rows, 1>::operator = (DenseStorage<T, Rows, 1>&& RHS_Object)
{
  DS_2DBase<T, Rows, 1>::operator = (static_cast<DS_2DBase<T, Rows, 1>&&>(RHS_Object));
  return *this;
}

template <typename T, std::size_t Rows>
inline typename T& DenseStorage<T, Rows, 1>::operator [] (std::size_t nIdx) 
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (m_pDataMassive == nullptr || nIdx >= GetNoRows())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  return *(m_pDataMassive + nIdx);
}

template <typename T, std::size_t Rows>
inline const typename T& DenseStorage<T, Rows, 1>::operator [] (std::size_t nIdx) const
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (m_pDataMassive == nullptr || nIdx >= GetNoRows())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  return *(m_pDataMassive + nIdx);
}

template <typename T, std::size_t Rows>
inline typename T& DenseStorage<T, Rows, 1>::operator () (std::size_t nIdxRow)
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (m_pDataMassive == nullptr || nIdxRow >= GetNoRows())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  return *(m_pDataMassive + nIdxRow);
}

template <typename T, std::size_t Rows>
inline const typename T& DenseStorage<T, Rows, 1>::operator () (std::size_t nIdxRow) const
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (m_pDataMassive == nullptr || nIdxRow >= GetNoRows())
    throw Math::Basic::StorageException(Math::Basic::StorageException::Type::InvalidMemoryAcces);
#endif

  return *(m_pDataMassive + nIdxRow);
}

template <typename T, std::size_t Rows>
DenseStorage<T, 1, Rows> DenseStorage<T, Rows, 1>::GetRowVector()
{
  DenseStorage<T, 1, Rows> Transposed;
  for (std::size_t k = 0; k < GetLinearSize(); ++k)
    Transposed[k] = (*this)[k];
  return Transposed;
}

template <typename T, std::size_t Rows>
std::size_t DenseStorage<T, Rows, 1>::GetNoRows() const
{
  return Rows;
}

template <typename T, std::size_t Rows>
std::size_t DenseStorage<T, Rows, 1>::GetNoCols() const
{
  return 1;
}

template <typename T, std::size_t Rows>
std::size_t DenseStorage<T, Rows, 1>::GetLinearSize() const
{
  return const_cast<DenseStorage<T, Rows, 1>*>(this)->m_pDataMassive ? GetNoCols() * GetNoRows() : 0;
}

template <typename T, std::size_t Rows>
inline DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize> DenseStorage<T, Rows, 1>::CreateRTDynamicStorage() const
{
  return DenseStorage<T, Misc::DynamicSize, Misc::DynamicSize>(GetNoRows(), GetNoCols());
}

