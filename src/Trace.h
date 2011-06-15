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

using namespace std;

struct cmp_str
{
  bool operator()(const char* a, const char* b)
  {
  	return strcmp(a,b) < 0;
  }
};

class Trace{
public:


	Trace();
	~Trace();
	
	map<const char*, int, cmp_str> traceTypes;

	void addTraceType(const char*);
	void trace(const char*, const char*, ...);

	

};




#endif

