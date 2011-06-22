/**
 * Trace.h
 *
 * Manage the output of trace messages.
 *
 */

#ifndef TRACE_H_
#define TRACE_H_


#include <cstdarg>
#include <cstring>
#include <map>
#include <iostream>
#include <fstream>


using namespace std;

//the map structure needs to know how to evaluate if two keys are the same
struct cmp_str
{
  bool operator()(const char* a, const char* b)
  {
  	// strcmp returns a negative value if two strings are the same
	return strcmp(a,b) < 0;
  }
};

class Trace{
public:

	Trace();
	Trace(const char*);

	~Trace();

	FILE* traceFile;

	// char* keys, int values	
	// a tag is enabled if the value is nonzero
	map<const char*, int, cmp_str> traceTypes;

	void addTraceType(const char*, int);
	#ifdef NOTRACING 
	inline void trace(const char*, const char*, ...){0;}
	#else
	void trace(const char*, const char*, ...);
	#endif
	FILE* getTraceFile();
	FILE* setTraceFile(FILE*);

	void enableTraceType(const char*);
	void disableTraceType(const char*);
	
	
};




#endif

