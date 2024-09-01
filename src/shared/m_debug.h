/*
!! ------------------------------------------------------------------------- !!
!!                                                          -*- mode:f90 -*-
!! Header file for m_debug.F90
!!
!! ------------------------------------------------------------------------- !!
*/

#if defined(_DEBUG) || defined(DEBUG)
#if defined(__INTEL_COMPILER)
#define debug(a)   debug__macro( a, #a, __FILE__, __LINE__ )
#define DEBUG(a)   debug__macro( a, #a, __FILE__, __LINE__ )
#else              
#define debug(a)   debug__macro( a, __FILE__, __LINE__ )
#define DEBUG(a)   debug__macro( a, __FILE__, __LINE__ )
#endif             
#else              
#define debug(a)   debug__void()
#define DEBUG(a)   debug__void()
#endif

#if defined(_INFO) || defined(INFO)
#define info(a)    info__macro( a, __FILE__, __LINE__ )
#define INFO(a)    info__macro( a, __FILE__, __LINE__ )
#else              
#define info(a)    debug__void()
#define INFO(a)    debug__void()
#endif

#if defined(__INTEL_COMPILER)
#define assert(a)  assert__macro( a, #a, __FILE__, __LINE__ )
#define ASSERT(a)  assert__macro( a, #a, __FILE__, __LINE__ )
#else 
#define assert(a)  assert__macro( a, __FILE__, __LINE__ )
#define ASSERT(a)  assert__macro( a, __FILE__, __LINE__ )
#endif             
