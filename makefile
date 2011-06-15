CC	= g++
CFLAGS	= -g -O0 -Wall
IFLAGS = -I"../include"
LFLAGS = -L"../lib"


OBJ_FILE	= Main.o Cell.o Experiment.o DerivGraph.o Trace.o
EXE_FILE	= EvoDevo
OUTPUT_DIR	= ./output


${EXE_FILE}: ${OBJ_FILE}
	${CC} ${LFLAGS} -o ${EXE_FILE} ${OBJ_FILE} 

Main.o: Main.cpp
	${CC} ${IFLAGS} ${CFLAGS} -c Main.cpp

Cell.o: Cell.cpp Cell.h
	${CC} ${IFLAGS} ${CFLAGS} -c Cell.cpp

DerivGraph.o: DerivGraph.cpp DerivGraph.h
	${CC} ${IFLAGS} ${CFLAGS} -c DerivGraph.cpp

Experiment.o: Experiment.cpp Experiment.h
	${CC} ${IFLAGS} ${CFLAGS} -c Experiment.cpp

Trace.o: Trace.cpp Trace.h
	${CC} ${IFLAGS} ${CFLAGS} -c Trace.cpp

clean:
	rm -rf ${OBJ_FILE} ${EXE_FILE} ${OUTPUT_DIR} 
