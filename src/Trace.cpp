
#include <cstdio>

#include "Trace.h"


Trace::Trace()
{



}


Trace::~Trace()
{



}

void Trace::addTraceType(const char* type){

	traceTypes[type] = 1;
}

void Trace::trace(const char* type, const char* format, ...){


	
	va_list args;
	va_start(args, format);
	if(traceTypes[type]){
		printf("[ %s ] \t",type);//, __FILE__, __LINE__);
		vprintf(format, args);
	}	
	va_end(args);

}
