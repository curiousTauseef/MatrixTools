template<class Scalar>
Scalar Abs(const Scalar& X)
{
  if (X > 0 || X == 0)
    return X;
  else
    return X * (-1);
}

template<class Scalar, bool ABS>
inline Scalar Min(const Scalar& A, const Scalar& B)
{
  if (ABS)
    return (Abs<Scalar>(A) < Abs<Scalar>(B)) ? A : B;
  return (A < B) ? A : B;
}

template<class Scalar, bool ABS>
inline Scalar Max(const Scalar& A, const Scalar& B)
{
  if (ABS)
    return (Abs<Scalar>(A) > Abs<Scalar>(B)) ? A : B;
  return (A > B) ? A : B;
}

template <class Scalar>
int ValueCmp(const Scalar& A, const Scalar& B, const Scalar& Tolerance)
{
  if (A > B + Tolerance)
    return 1;
  else if (B > A + Tolerance)
    return -1;
  else 
    return 0;
}

template<class Scalar>
int ValueCmpZero(const Scalar& A, const Scalar& Tolerance)
{
  if (ValueCmp<Scalar>(A, 0, Tolerance) == 1)
    return 1;
  else if (ValueCmp<Scalar>(A, 0, Tolerance) == -1)
    return -1;
  else 
    return 0;
}

template<class Scalar>
bool ValueBtw(const Scalar& ASup, const Scalar& BInf, const Scalar& X, 
              const Scalar& Tolerance, compType eCompType /* = compType::EQ */)
{
  switch (eCompType)
  {
  case compType::EQ_SUP :
    {
      if (ValueCmp(X, BInf, Tolerance) == 1 && ValueCmp(X, ASup, Tolerance) <= 0)
        return true;
      else
        return false;
    }
    break;
  case compType::EQ_INF :
    {
      if (ValueCmp(X, BInf, Tolerance) >= 0 && ValueCmp(X, ASup, Tolerance) == -1)
        return true;
      else
        return false;
    }
    break;
  case compType::EQ :
    {
      if (ValueCmp(X, BInf, Tolerance) >= 0 && ValueCmp(X, ASup, Tolerance) <= 0)
        return true;
      else
        return false;
    }
    break;
  }
  return false;
}

template<class Scalar>
inline Scalar Square(const Scalar& A)
{ 
  return A * A; 
}

template<class Scalar>
inline Scalar Cube(const Scalar& A) 
{ 
  return A * A * A ; 
}

template<class Scalar>
inline Scalar Deg_inRad(const Scalar& X) 
{ 
  static_assert(std::is_floating_point<Scalar>::value && "Only basic floating point types are usable");
  return X * (Pi / 180.); 
}

template<class Scalar>
inline Scalar Rad_inDeg(const Scalar& X) 
{ 
  static_assert(std::is_floating_point<Scalar>::value && "Only basic floating point types are usable");
  return X * (180. / Pi); 
}

template<class Scalar>
Scalar Interpoli_Linear2D(const Scalar& X1, const Scalar& Y1, 
                          const Scalar& X2, const Scalar& Y2, const Scalar& Xc, 
					                const Scalar& Tolerance)
{
  static_assert(std::is_floating_point<Scalar>::value && "Only basic floating point types are usable");
  Scalar Yc;
  if (ValueCmp(Xc - X1, Tolerance) != 0 && ValueCmp(Xc - X1, Tolerance) != 0)
    Yc = Y1 + (Y2 - Y1) * (Xc - X1) / (X2 - X1);
  else
    Yc = Y2;
  return Yc;
}
