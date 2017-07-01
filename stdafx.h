// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
// Windows Header Files:
#include <windows.h>

// TODO: reference additional headers your program requires here

// C/C++
#include <cassert>
// STL
#include <vector>
#include <list>
#include <set>
#include <map>
#include <queue>
#include <deque>
#include <stack>
#include <string>
#include <bitset>
#include <algorithm>
#include <memory>
#include <iostream>
#include <type_traits>
#include <sstream>

// others
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif