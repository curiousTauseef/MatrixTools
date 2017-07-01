template <typename Type>
RawIterator<Type>::RawIterator(Type* Ptr) : m_CrtPtr(Ptr)
{
}

template <typename Type>
RawIterator<Type>::RawIterator(const RawIterator<Type>& RHS) : m_CrtPtr(RHS.m_CrtPtr)
{
}

template <typename Type>
RawIterator<Type>::RawIterator(const RawIterator<typename AddConst<Type>::TypeWithCst>& RHS) : 
  m_CrtPtr(const_cast<RawIterator<typename AddConst<Type>::TypeWithCst>&>(RHS).GetWrappedPtr())  // Just to keep things in order
{
}

template <typename Type>
RawIterator<Type>::~RawIterator()
{
}

template <typename Type>
inline RawIterator<Type>& RawIterator<Type>::operator = (const RawIterator<Type>& RHS)
{
  m_CrtPtr = RHS.m_CrtPtr;
  return *this;
}

template <typename Type>
RawIterator<Type>& RawIterator<Type>::operator = (const RawIterator<typename AddConst<Type>::TypeWithCst>& RHS)
{
  m_CrtPtr = const_cast<RawIterator<typename AddConst<Type>::TypeWithCst>&>(RHS).GetWrappedPtr();
  return *this;
}

template <typename Type>
inline RawIterator<Type>& RawIterator<Type>::operator = (Type* TypePtr)
{
  m_CrtPtr = TypePtr;
  return *this;
}

template <typename Type>
inline RawIterator<Type>& RawIterator<Type>::operator = (typename AddConst<Type>::TypeWithCst* TypePtr)
{
  m_CrtPtr = TypePtr;
  return *this;
}

template <typename Type>
inline RawIterator<Type>::operator bool () const
{
  return m_CrtPtr != nullptr ? true : false;
}

template <typename Type>
inline bool RawIterator<Type>::operator == (const RawIterator<Type>& Iterator) const
{
  return m_CrtPtr == Iterator.m_CrtPtr;
}

template <typename Type>
inline bool RawIterator<Type>::operator != (const RawIterator<Type>& Iterator) const
{
  return m_CrtPtr != Iterator.m_CrtPtr;
}

template <typename Type>
inline const RawIterator<Type> RawIterator<Type>::operator + (const ptrdiff_t& Delta) const
{
  return RawIterator(m_CrtPtr + Delta);
}

template <typename Type>
inline const RawIterator<Type> RawIterator<Type>::operator - (const ptrdiff_t& Delta) const
{
  return RawIterator(m_CrtPtr - Delta);
}

template <typename Type>
inline ptrdiff_t RawIterator<Type>::operator - (const RawIterator<Type>& Iterator) const
{
  return m_CrtPtr - Iterator.m_CrtPtr;
}

template <typename Type>
inline RawIterator<Type>& RawIterator<Type>::operator += (const ptrdiff_t& Delta)
{
  m_CrtPtr = m_CrtPtr + Delta;
  return *this;
}

template <typename Type>
inline RawIterator<Type>& RawIterator<Type>::operator -= (const ptrdiff_t& Delta)
{
  m_CrtPtr = m_CrtPtr - Delta;
  return *this;
}

template <typename Type>
inline RawIterator<Type>& RawIterator<Type>::operator ++ ()
{
  m_CrtPtr = m_CrtPtr + 1;
  return *this;
}

template <typename Type>
inline RawIterator<Type>& RawIterator<Type>::operator -- ()
{
  m_CrtPtr = m_CrtPtr - 1;
  return *this;
}

template <typename Type>
inline RawIterator<Type> RawIterator<Type>::operator ++ (int)
{
  RawIterator<Type> Copy(*this);
  m_CrtPtr = m_CrtPtr + 1;
  return Copy;
}

template <typename Type>
inline RawIterator<Type> RawIterator<Type>::operator -- (int)
{
  RawIterator<Type> Copy(*this);
  m_CrtPtr = m_CrtPtr - 1;
  return Copy;
}

template <typename Type>
inline Type& RawIterator<Type>::operator * () const
{
  return *m_CrtPtr ;
}

template <typename Type>
inline Type* RawIterator<Type>::operator -> () const
{
  return m_CrtPtr ;
}

template <typename Type>
inline Type& RawIterator<Type>::operator [] (ptrdiff_t Delta) const
{
  return *(m_CrtPtr + Delta) ;
}

template <typename Type>
inline Type* RawIterator<Type>::GetWrappedPtr()
{
  return m_CrtPtr;
}

template <typename Type>
RevRawIterator<Type>::RevRawIterator(Type* Ptr = nullptr) : RawIterator(Ptr)
{
}

template <typename Type>
RevRawIterator<Type>::RevRawIterator(const RevRawIterator<Type>& RHS) : RawIterator(RHS)
{
}

template <typename Type>
RevRawIterator<Type>::RevRawIterator(const RevRawIterator<typename AddConst<Type>::TypeWithCst>& RHS) : RawIterator(RHS)
{
}

template <typename Type>
RevRawIterator<Type>::~RevRawIterator()
{
}

template <typename Type>
inline RevRawIterator<Type>& RevRawIterator<Type>::operator = (const RevRawIterator<Type>& RHS)
{
  RawIterator::operator = (RHS);
  return *this;
}

template <typename Type>
RevRawIterator<Type>& RevRawIterator<Type>::operator = (const RevRawIterator<typename AddConst<Type>::TypeWithCst>& RHS)
{
  RawIterator::operator = (RHS);
  return *this;
}

template <typename Type>
inline RevRawIterator<Type>& RevRawIterator<Type>::operator = (Type* TypePtr)
{
  RawIterator::operator = (TypePtr);
  return *this;
}

template <typename Type>
inline RevRawIterator<Type>& RevRawIterator<Type>::operator = (typename AddConst<Type>::TypeWithCst* TypePtr)
{
  RawIterator::operator = (TypePtr);
  return *this;
}

template <typename Type>
inline bool RevRawIterator<Type>::operator == (const RevRawIterator<Type>& Iterator) const
{
  return m_CrtPtr == Iterator.m_CrtPtr;
}

template <typename Type>
inline bool RevRawIterator<Type>::operator != (const RevRawIterator<Type>& Iterator) const
{
  return m_CrtPtr != Iterator.m_CrtPtr;
}

template <typename Type>
inline const RevRawIterator<Type> RevRawIterator<Type>::operator + (const ptrdiff_t& Delta) const
{
  return RevRawIterator(m_CrtPtr - Delta);
}

template <typename Type>
inline const RevRawIterator<Type> RevRawIterator<Type>::operator - (const ptrdiff_t& Delta) const
{
  return RevRawIterator(m_CrtPtr + Delta);
}

template <typename Type>
inline ptrdiff_t RevRawIterator<Type>::operator - (const RevRawIterator<Type>& Iterator) const
{
  return Iterator.m_CrtPtr - m_CrtPtr;
}

template <typename Type>
inline RevRawIterator<Type>& RevRawIterator<Type>::operator += (const ptrdiff_t& Delta)
{
  m_CrtPtr = m_CrtPtr - Delta;
  return *this;
}

template <typename Type>
inline RevRawIterator<Type>& RevRawIterator<Type>::operator -= (const ptrdiff_t& Delta)
{
  m_CrtPtr = m_CrtPtr + Delta;
  return *this;
}

template <typename Type>
inline RevRawIterator<Type>& RevRawIterator<Type>::operator ++ ()
{
  m_CrtPtr = m_CrtPtr - 1;
  return *this;
}

template <typename Type>
inline RevRawIterator<Type>& RevRawIterator<Type>::operator -- ()
{
  m_CrtPtr = m_CrtPtr + 1;
  return *this;
}

template <typename Type>
inline RevRawIterator<Type> RevRawIterator<Type>::operator ++ (int)
{
  RevRawIterator<Type> Copy(*this);
  m_CrtPtr = m_CrtPtr - 1;
  return Copy;
}

template <typename Type>
inline RevRawIterator<Type> RevRawIterator<Type>::operator -- (int)
{
  RevRawIterator<Type> Copy(*this);
  m_CrtPtr = m_CrtPtr + 1;
  return Copy;
}

template <typename Type>
inline Type& RevRawIterator<Type>::operator [] (ptrdiff_t Delta) const
{
  return *(m_CrtPtr - Delta) ;
}