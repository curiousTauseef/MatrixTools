template <std::size_t Dim>
Perm_Base<Dim>::Perm_Base(std::size_t Size /* = Dim */) : m_Array(new std::size_t [Size])
{
  for (std::size_t k = 0; k < Size; ++k)     // identical permutation by default starting from zero
    *(m_Array + k) = k;
}

template <std::size_t Dim>
Perm_Base<Dim>::Perm_Base(const std::map<std::size_t, std::size_t>& Mapping)
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (!Mapping.size() || Mapping.rbegin()->first >= Mapping.size())
    throw Basic::GeneralException(Basic::GeneralException::Type::InvalidParams);
  std::set<std::size_t> VisitedNodes;
#endif

  m_Array = new std::size_t [Mapping.size()];

  std::map<std::size_t, std::size_t>::const_iterator itCrt(Mapping.begin()), itEnd(Mapping.end());
  for (; itCrt != itEnd; ++itCrt)
  {
    m_Array[itCrt->first] = itCrt->second;
#ifndef MATH_SAFETY_CHECKS_OFF
    if (!(VisitedNodes.insert(itCrt->second)).second)
      throw Basic::GeneralException(Basic::GeneralException::Type::InvalidParams);   // duplicate elements
#endif
  }
}

template <std::size_t Dim>
Perm_Base<Dim>::Perm_Base(const Perm_Base<Dim>& RHS)
{
  m_Array = new std::size_t [RHS.GetSize()];
  for (std::size_t k = 0; k < RHS.GetSize(); ++k)
    *(m_Array + k) = *(RHS.m_Array + k);
}

template <std::size_t Dim>
Perm_Base<Dim>::Perm_Base(Perm_Base<Dim>&& RHS)
{
  m_Array = RHS.m_Array;
  RHS.m_Array = nullptr;
}

template <std::size_t Dim>
inline Perm_Base<Dim>::~Perm_Base()
{
  if (m_Array)
    delete [] m_Array;
}

template <std::size_t Dim>
Perm_Base<Dim>& Perm_Base<Dim>::operator = (const Perm_Base<Dim>& RHS)
{
  if (this != &RHS)
  {
    m_Array = new std::size_t [RHS.GetSize()];
    for (std::size_t k = 0; k < RHS.GetSize(); ++k)
      *(m_Array + k) = *(RHS.m_Array + k);
  }
  return *this;
}

template <std::size_t Dim>
Perm_Base<Dim>& Perm_Base<Dim>::operator = (Perm_Base<Dim>&& RHS)
{
  if (this != &RHS)
  {
    if (m_Array)
      delete [] m_Array;

    m_Array = RHS.m_Array;
    RHS.m_Array = nullptr;
  }
  return *this;
}

template <std::size_t Dim>
inline std::size_t Perm_Base<Dim>::operator () (std::size_t Index) const
{
  if (Index < GetSize())
    return m_Array[Index];
  return m_Array[GetSize() - 1];     // else return the last element
}

template <std::size_t Dim>
Perm_Base<Dim> Perm_Base<Dim>::operator * (const Perm_Base<Dim>& RHS) const
{
  Perm_Base<Dim> Result;  // identity for now

  for (std::size_t i = 0; i < GetSize(); ++i)
    Result.m_Array[i] = m_Array[RHS.m_Array[i]]; 

  return Result;
}

template <std::size_t Dim>
bool Perm_Base<Dim>::operator == (const Perm_Base<Dim>& RHS) const
{
  for (std::size_t i = 0; i < GetSize(); ++i)
    if (m_Array[i] != RHS.m_Array[i])
      return false;

  return true;
}

template <std::size_t Dim>
bool Perm_Base<Dim>::operator != (const Perm_Base<Dim>& RHS) const
{
  return !((*this) == RHS);
}

template <std::size_t Dim>
void Perm_Base<Dim>::GetInverse(Perm_Base<Dim>& Inverse) const
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (Inverse.GetSize() != GetSize())
    throw Basic::GeneralException(Basic::GeneralException::Type::InvalidParams);
#endif

  Inverse = *this;

  for (std::size_t i = 0; i < GetSize(); ++i)
    Inverse.m_Array[Inverse.m_Array[i]] = i; 
}

template <std::size_t Dim>
std::size_t Perm_Base<Dim>::SearchValue(std::size_t Value) const
{
#ifndef MATH_SAFETY_CHECKS_OFF
  if (Value >= GetSize())
    throw Basic::GeneralException(Basic::GeneralException::Type::InvalidParams);
#endif

  for (std::size_t i = 0; i < GetSize(); ++i)
  {
    if (m_Array[i] == Value)
      return i;
  }

  return GetSize() - 1;      // I can't actually get here but I want to suppress any compiler warnings...
}

template <std::size_t Dim>
template <typename PermType>
void Perm_Base<Dim>::DisjointCycles(std::list<PermType>& CyclesList) const
{
  bool bContinue = true, bSkip = false;
  std::vector<bool> Visited;
  Visited.reserve(GetSize());

  std::size_t k = 0, nVisited = 0;

  for (std::size_t i = 0; i < GetSize(); ++i)
    Visited.push_back(false);

  PermType CrtCycle(GetSize());

  for (std::size_t j = 0; j < GetSize(); ++j)   // Go to the next unvisited cycle
  {
    if (!Visited[j])
    {
      CrtCycle.ResetToIdentity();
      k = j;
      do 
      {
        CrtCycle.m_Array[k] = m_Array[k];
        Visited[k] = true;
        k = m_Array[k];
      } while (k != j);
      CyclesList.push_back(CrtCycle);
    }
  }
}

template <std::size_t Dim>
Perm_Base<Dim>& Perm_Base<Dim>::ResetToIdentity()
{
  for (std::size_t k = 0; k < GetSize(); ++k)
    *(m_Array + k) = k;
  return *this;
}

template <std::size_t Dim>
Perm_Base<Dim>& Perm_Base<Dim>::Permute(std::size_t i, std::size_t j)
{
  if (i < GetSize() && j < GetSize() && i != j)
  {
    *(m_Array + i) = *(m_Array + i) ^ *(m_Array + j);
    *(m_Array + j) = *(m_Array + i) ^ *(m_Array + j);
    *(m_Array + i) = *(m_Array + i) ^ *(m_Array + j);
  }
  return *this;
}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
template <std::size_t Dim>
Permutation<Dim>::Permutation(std::size_t Size /* = Dim */) : Perm_Base<Dim>(Size)
{
}

template <std::size_t Dim>
Permutation<Dim>::Permutation(const std::map<std::size_t, std::size_t>& Mapping) : Perm_Base<Dim>(Mapping)
{
}

template <std::size_t Dim>
Permutation<Dim>::Permutation(const Permutation<Dim>& RHS) : Perm_Base<Dim>(RHS)
{
}

template <std::size_t Dim>
Permutation<Dim>::Permutation(Permutation<Dim>&& RHS) : Perm_Base<Dim>(static_cast<Perm_Base<Dim>&&>(RHS))
{
}

template <std::size_t Dim>
inline Permutation<Dim>::~Permutation()
{
  // Base will handle things correctly...
}

template <std::size_t Dim>
Permutation<Dim>& Permutation<Dim>::operator = (const Permutation<Dim>& RHS)
{
  Perm_Base::operator = (RHS);
  return *this;
}

template <std::size_t Dim>
Permutation<Dim>& Permutation<Dim>::operator = (Permutation<Dim>&& RHS)
{
  Perm_Base::operator = (static_cast<Perm_Base<Dim>&&>(RHS));
  return *this;
}

template <std::size_t Dim>
std::size_t Permutation<Dim>::GetSize() const
{
  return Dim;
}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
Permutation<Misc::DynamicSize>::Permutation(std::size_t Size) : Perm_Base<Misc::DynamicSize>(Size)
{
  m_Size = Size > 0 ? Size : 0;
}

Permutation<Misc::DynamicSize>::Permutation(const std::map<std::size_t, std::size_t>& Mapping) : Perm_Base<Misc::DynamicSize>(Mapping)
{
}

Permutation<Misc::DynamicSize>::Permutation(const Permutation<Misc::DynamicSize>& RHS) : Perm_Base<Misc::DynamicSize>(RHS), m_Size(RHS.GetSize())
{
}

Permutation<Misc::DynamicSize>::Permutation(Permutation<Misc::DynamicSize>&& RHS) : Perm_Base<Misc::DynamicSize>(static_cast<Perm_Base<Misc::DynamicSize>&&>(RHS)), m_Size(RHS.GetSize())
{
  RHS.m_Size = 0;
}

Permutation<Misc::DynamicSize>::~Permutation()
{
   // Base will take care of things of course...
}

Permutation<Misc::DynamicSize>& Permutation<Misc::DynamicSize>::operator = (const Permutation<Misc::DynamicSize>& RHS) 
{
  if (this != &RHS)
  {
    Perm_Base<Misc::DynamicSize>::operator = (RHS);
    m_Size = RHS.GetSize();
  }
  return *this;
}

Permutation<Misc::DynamicSize>& Permutation<Misc::DynamicSize>::operator = (Permutation<Misc::DynamicSize>&& RHS)
{
  if (this != &RHS)
  {
    Perm_Base<Misc::DynamicSize>::operator = (static_cast<Perm_Base<Misc::DynamicSize>&&>(RHS));
    RHS.m_Size = 0;
  }
  return *this;
}

std::size_t Permutation<Misc::DynamicSize>::GetSize() const
{
  return m_Size;
}

Permutation<Misc::DynamicSize>& Permutation<Misc::DynamicSize>::Resize(std::size_t nNewSize)
{
  assert(nNewSize > 0 && L"The size must be a positive integer!");

  if (nNewSize != GetSize() && nNewSize > 0)
  {
    if (m_Array)
      delete [] m_Array;

    m_Size = nNewSize;
    m_Array = new std::size_t [GetSize()];
    ResetToIdentity();
  }

  return *this;
}