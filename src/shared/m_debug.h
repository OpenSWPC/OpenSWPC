#if defined(_DEBUG)
#define debug(a)   debug__macro(a, __FILE__, __LINE__)
#else              
#define debug(a)   debug__macro()
#endif
#define info(a)    info__macro(a, __FILE__, __LINE__)
#define assert(a)  assert__macro(a, __FILE__, __LINE__)
