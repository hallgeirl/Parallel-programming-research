#ifndef types_h
// Include this in every modules that uses floating point
// arithmetic, and declare all floating point values as "DOUBLE"
// With a switch of a command line macro set up in the Makefile
// we can then change the arithmetic
//
#define types_h
#ifdef FLOAT
#define DOUBLE float
#else
#define DOUBLE double
#endif

#endif

