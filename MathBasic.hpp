#ifndef MathFuncDLL_hpp__
#define MathFuncDLL_hpp__

// Math Header
#include "MathDefines.h"

#ifndef min
#define min(x, y) x < y ? x : y
#endif

#ifndef max
#define max(x, y) x > y ? x : y
#endif

MATH_START
MATH_BASIC_START

enum class compType : char { EQ_SUP = 1, EQ_INF = -1, EQ = 0 };

//////////////////////////////////////////////////////////////////////////
/// Get the absolute value of an object
/// @param[in]  X  - the actual object
/// @return     the absolute value of X
//////////////////////////////////////////////////////////////////////////
template<class Scalar>
Scalar Abs(const Scalar& X);
//////////////////////////////////////////////////////////////////////////
/// Get the maximum value
/// @param[in]  A         - the first value
/// @param[in]  B         - the second value
/// @return     the largest value from the 2
//////////////////////////////////////////////////////////////////////////
template<class Scalar, bool ABS>
Scalar Min(const Scalar& A, const Scalar& B);
//////////////////////////////////////////////////////////////////////////
/// Get the minimum value
/// @param[in]  A         - the first value
/// @param[in]  B         - the second value
/// @return     the smallest value from the 2
//////////////////////////////////////////////////////////////////////////
template<class Scalar, bool ABS>
Scalar Max(const Scalar& A, const Scalar& B);
//////////////////////////////////////////////////////////////////////////
/// Compare two types with a given tolerance.
/// Function returns: 1 if A > B ; -1 if A < B ; 0 if A = B.
/// @param[in]  A         - the first value
/// @param[in]  B         - the second value
/// @param[in]  Tolerance - the tolerance value
/// @return     the result of the comparison: 1 A > B ; -1 A < B ; 0 A = B
//////////////////////////////////////////////////////////////////////////
template <class Scalar>
int ValueCmp(const Scalar& A, const Scalar& B, const Scalar& Tolerance);
//////////////////////////////////////////////////////////////////////////
/// Compare a type with zero with a given tolerance.
/// Function returns: 1 if A > 0 ; -1 if A < 0 ; 0 if A = 0.
/// @param[in]  A         - the number to be compared
/// @param[in]  Tolerance - the tolerance value
/// @return     the result of the comparison: 1 A > 0 ; -1 A < 0 ; 0 A = 0
//////////////////////////////////////////////////////////////////////////
template<class Scalar>
int ValueCmpZero(const Scalar& A, const Scalar& Tolerance);
//////////////////////////////////////////////////////////////////////////
/// Test if a type is inside an interval.
/// Interval can be open/closed or any combination.
/// @param[in]  BInf      - the lower margin
/// @param[in]  ASup      - the upper margin
/// @param[in]  X         - the number to be tested
/// @param[in]  eCompType  - the type of margin (closed/open or variations)
/// @param[in]  Tolerance - the tolerance value
/// @return     the result of the test true/false
//////////////////////////////////////////////////////////////////////////
template<class Scalar>
bool ValueBtw(const Scalar& ASup, const Scalar& BInf, const Scalar& X, const Scalar& Tolerance, compType eCompType = compType::EQ);
//////////////////////////////////////////////////////////////////////////
/// Get the square of a type.
/// Compute X^2.
/// @param[in]	A - The object to raise.
/// @return		  The square of A.
//////////////////////////////////////////////////////////////////////////
template<class Scalar>
Scalar Square(const Scalar& A);
//////////////////////////////////////////////////////////////////////////
/// Get the cube of a type.
/// Compute X^3.
/// @param[in]	A - The object to raise.
/// @return		  The cube of A.
//////////////////////////////////////////////////////////////////////////
template<class Scalar>
Scalar Cube(const Scalar& A);
//////////////////////////////////////////////////////////////////////////
/// Convert angle from degrees to radians.
/// @param[in]	X - The angle value [degrees].
/// @return		  The angle in radians.
//////////////////////////////////////////////////////////////////////////
template<class Scalar>
Scalar Deg_inRad(const Scalar& X);
//////////////////////////////////////////////////////////////////////////
/// Convert angle from radians to degrees.
/// @param[in]	X - The angle value [radians].
/// @return		  The angle in degrees.
//////////////////////////////////////////////////////////////////////////
template<class Scalar>
Scalar Rad_inDeg(const Scalar& X);
//////////////////////////////////////////////////////////////////////////
/// Linear interpolation for a given object.
/// The function gets y = f(x) assuming a linear variation.
/// @param[in]  X1      - the lower margin x value
/// @param[in]  Y1      - the lower margin y value
/// @param[in]  X2      - the upper margin x value
/// @param[in]  Y2      - the upper margin y value
/// @param[in]  Xc      - the x value
/// @return     the f(Xc) value for a linear variation
//////////////////////////////////////////////////////////////////////////
template<class Scalar>
Scalar Interpoli_Linear2D(const Scalar& X1, const Scalar& Y1, const Scalar& X2, const Scalar& Y2, const Scalar& Xc, const Scalar& Tolerance);

#include "MathBasic.inl"

MATH_BASIC_END
MATH_END

#endif // MathFuncDLL_hpp__
