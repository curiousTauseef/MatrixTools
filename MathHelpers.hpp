#ifndef MathHelpers_hpp__
#define MathHelpers_hpp__

#include "MathDefines.h"
#include <cassert>

// Math helpers

MATH_START
MATH_MISC_START

enum { DynamicSize = 0 };

// - Trait for identifying compile time (static) vs runtime (dynamic) template parameters
// - returns ValTrue if Val is equal to Misc::DynamicSize ans ValFalse otherwise
template <std::size_t Val, std::size_t ValTrue, std::size_t ValFalse>
struct ValueIfRtDynamic
{
  static const std::size_t Value = ValFalse;
};

template <std::size_t ValTrue, std::size_t ValFalse>
struct ValueIfRtDynamic<Misc::DynamicSize, ValTrue, ValFalse>
{
  static const std::size_t Value = ValTrue;
};

template <bool CompileTimeExpr, typename TypeValIfTrue, typename TypeValIfFalse>
struct CompileIf
{
  typedef TypeValIfFalse RetValue;
};

template <typename TypeValIfTrue, typename TypeValIfFalse>
struct CompileIf<true, TypeValIfTrue, TypeValIfFalse>
{
  typedef TypeValIfTrue RetValue;
};

// - Generic base permutation class (can't be instantiated by itself since it is meant this way)
template <std::size_t Dim>
class Perm_Base
{
protected :
  /// Ctor (has a parameter in order to have a common interface)
  explicit Perm_Base(std::size_t Size = Dim);
  /// Ctor (explicitly specify the pairs of values)
  explicit Perm_Base(const std::map<std::size_t, std::size_t>& Mapping);
  /// Copy ctor
  Perm_Base(const Perm_Base<Dim>& RHS);
  /// Move ctor
  Perm_Base(Perm_Base<Dim>&& RHS);
public :
  /// Dtor
  virtual ~Perm_Base();

  // - Operators

  /// Assignment
  Perm_Base<Dim>& operator = (const Perm_Base<Dim>& RHS);
  /// Move assignment
  Perm_Base<Dim>& operator = (Perm_Base<Dim>&& RHS);
  /// Access value (as a copy)
  std::size_t operator () (std::size_t Index) const;
  /// Compose 2 permutations -> s * t (i) = s(t(i)) 
  Perm_Base<Dim> operator * (const Perm_Base<Dim>& RHS) const;
  /// Test if 2 permutations are equal
  bool operator == (const Perm_Base<Dim>& RHS) const;
  /// Test if 2 permutations are different
  bool operator != (const Perm_Base<Dim>& RHS) const;

  // - Accessors

  /// Get the size
  virtual std::size_t GetSize() const = 0;
  /// Compute the inverse of the current permutation and insert it into the specified object (the sizes must match)
  void GetInverse(Perm_Base<Dim>& Inverse) const;
  /// Search for the specified value and return the corresponding index
  std::size_t SearchValue(std::size_t Value) const;
  /// Decompose the permutation into disjoint cycles (will only compile when PermType is a valid permutation)
  template <typename PermType>
  void DisjointCycles(std::list<PermType>& CyclesList) const;

  // - Mutators

  /// Reset to identity
  Perm_Base<Dim>& ResetToIdentity();
  /// Permute 2 positions
  Perm_Base<Dim>& Permute(std::size_t i, std::size_t j);

protected :
  std::size_t* m_Array;
};

// - Generic class for handling permutations
template <std::size_t Dim = Misc::DynamicSize>
class Permutation : public Perm_Base<Dim>
{
public :
  /// Ctor (has a parameter in order to have a common interface)
  Permutation(std::size_t Size = Dim);
  /// Ctor (explicitly specify the pairs of values the values (first and second) must be < Size and unique!)
  explicit Permutation(const std::map<std::size_t, std::size_t>& Mapping);
  /// Copy ctor
  Permutation(const Permutation<Dim>& RHS);
  /// Move ctor
  Permutation(Permutation<Dim>&& RHS);
  /// Dtor
  ~Permutation();

  // - Accessors
  /// Get the size
  std::size_t GetSize() const;

  // - Operators

  /// Assignment
  Permutation<Dim>& operator = (const Permutation<Dim>& RHS);
  /// Move assignment
  Permutation<Dim>& operator = (Permutation<Dim>&& RHS);
};

// - Partial specialization for handling a runtime variable number of elements
template<>
class Permutation<Misc::DynamicSize> : public Perm_Base<Misc::DynamicSize>
{
public :
  /// Ctor
  explicit Permutation(std::size_t Size); 
  /// Ctor (explicitly specify the pairs of values the values (first and second) must be < Size and unique!)
  explicit Permutation(const std::map<std::size_t, std::size_t>& Mapping);
  /// Copy ctor
  Permutation(const Permutation<Misc::DynamicSize>& RHS);
  /// Move ctor
  Permutation(Permutation<Misc::DynamicSize>&& RHS);
  /// Dtor
  ~Permutation();

  // - Operators

  /// Assignment
  Permutation<Misc::DynamicSize>& operator = (const Permutation<Misc::DynamicSize>& RHS);
  /// Move assignment
  Permutation<Misc::DynamicSize>& operator = (Permutation<Misc::DynamicSize>&& RHS);

  // - Accessors

  /// Get the size
  std::size_t GetSize() const;

  // - Mutators

  /// Resize
  Permutation<Misc::DynamicSize>& Resize(std::size_t nNewSize);
private :
  std::size_t m_Size;
};

typedef Permutation<> DynamicPermutation;

#include "MathHelpers.inl"

MATH_MISC_END
MATH_END


#endif // MathHelpers_h__
