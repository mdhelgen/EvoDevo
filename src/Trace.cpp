/**
 * Trace.cpp
 *
 * Trace implementation.
 *
 * Trace messages are called with a tag that can be optionally turned on or off with a simple call.
 *
 * This allows trace messages to be given context and turned on or off very flexibly.
 *
 */

#include <cstdio>
#include "Trace.h"
#include <ostream>
/**
 * Trace::Trace()
 *
 * Default Constructor.
 *
 */
Trace::Trace()
{
 
 printf("Tracing loaded. (location %u)\n",(unsigned int) this);
 traceFile = stdout;

}

Trace::Trace(const char* c)
{
	traceFile = fopen(c, "a");
	
}

/**
 * Trace::~Trace()
 *
 * Default Destructor.
 */
Trace::~Trace()
{


}

/**
 * void Trace::addTraceType(const char*, int)
 *
 * Adds a new trace tag and enables it.
 *
 *  e.g.
 *    Trace::addTraceType("output")
 *    //output related
 *    t.trace("output", "Output complete");
 *
 * @param tag The tag to be used for tracing
 * @param enabled Initial state of the trace type. Nonzero is enabled.
 */
void Trace::addTraceType(const char* tag, int enabled){

	if(enabled)
		Trace::enableTraceType(tag);
	else
		Trace::disableTraceType(tag);
}

/**
 * void Trace::trace(const char*, const char*, ...)
 *
 * Outputs a trace message with the given format if the trace tag is enabled.
 * Output is not automatically terminated with a newline character.
 * 
 * @param tag Trace type
 * @param format string
 * @param ... variable arguments corresponding to format string
 */
#ifndef NOTRACING
void Trace::trace(const char* tag, const char* format, ...){

	va_list args;
	// variable arguments start after format
	va_start(args, format);

	//if the value associated with the tag is nonzero it is enabled
	if(traceTypes[tag]){
		
		//prefix message with tag
		fprintf(traceFile,"[ %-6s ] \t",tag);
		
		//output formatted trace message
		vfprintf(traceFile,format, args);
	}	
	
	va_end(args);

}
#endif

/**
 * Trace::enableTraceType(const char*)
 *
 * Enable the trace type, causing future trace messages tagged with this type to be output.
 * 
 * @param tag the trace tag to enable.
 */
void Trace::enableTraceType(const char* tag){
	
	//set the value with this key to 1
	traceTypes[tag] = 1;

	Trace::trace("trce","Trace type \'%s\' enabled.\n", tag);

}
/**
 * Trace::disableTraceType(const char*)
 *
 * Disable the trace type, causing future trace messages tagged with this type to be surpressed.
 *
 * @param tag the trace tag to disable.
 */
void Trace::disableTraceType(const char* tag){

	//set the value with this key to 0
	traceTypes[tag] = 0;

	Trace::trace("trce","Trace type \'%s\' disabled.\n", tag);

}

FILE* Trace::getTraceFile(){
	return traceFile;
}

FILE* Trace::setTraceFile(FILE * tf){
	traceFile = tf;
	return traceFile;
}
