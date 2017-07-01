#ifndef MathDefines_h__
#define MathDefines_h__

#ifdef MATHFUNCDLL_EXPORTS
  #define MATHDLL_API __declspec(dllexport)
#else
  #define MATHDLL_API __declspec(dllexport)
#endif // MATHFUNCDLL_EXPORTS

// Namespace reserved for math methods
#define MATH_START namespace Math {
#define MATH_END                  }

// Namespace reserved for basic mathematical functions
#define MATH_BASIC_START namespace Basic {
#define MATH_BASIC_END                   }

// Namespace reserved for misc utilities
#define MATH_MISC_START  namespace Misc {
#define MATH_MISC_END                   }

// Namespace reserved for general matrices
#define MATH_MATRIX_START namespace Matrix {
#define MATH_MATRIX_END                    }

// Namespace reserved for storage
#define MATH_STORAGE_START namespace Storage {
#define MATH_STORAGE_END                     }

// Namespace reserved for test sequences
#define MATH_MATRIX_TEST_START  namespace MatrixTest {
#define MATH_MATRIX_END                              }

// Parallelize for loops with Open MP
#ifndef _DEBUG                                    
  #ifdef  _OPENMP                                 
     #define PARALLEL_BLOCK #pragma omp parallel        
  #else                                           
     #define PARALLEL_BLOCK                       
  #endif                                          
#else                                             
  #define PARALLEL_BLOCK                          
#endif                         

#endif // MathDefines_h__
