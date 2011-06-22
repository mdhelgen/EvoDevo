CC	= g++
CFLAGS	= -g -O2 -Wall -DNOTRACING
IFLAGS = -I"../include"
LFLAGS = -L"../lib"
VPATH =../src

OBJ_FILE	= Main.o Cell.o Experiment.o DerivGraph.o Interaction.o Molecule.o CustomMolecules.o CustomInteractions.o Trace.o
EXE_FILE	= EvoDevo
OUTPUT_DIR	= ./output


${EXE_FILE}: ${OBJ_FILE}
	${CC} ${LFLAGS} -o ${EXE_FILE} ${OBJ_FILE} 

Main.o: ${VPATH}/Main.cpp
	${CC} ${IFLAGS} ${CFLAGS} -c $< -o $@ 

Cell.o: ${VPATH}/Cell.cpp ${VPATH}/Cell.h
	${CC} ${IFLAGS} ${CFLAGS} -c ${VPATH}/Cell.cpp

DerivGraph.o: ${VPATH}/DerivGraph.cpp ${VPATH}/DerivGraph.h
	${CC} ${IFLAGS} ${CFLAGS} -c ${VPATH}/DerivGraph.cpp

Experiment.o: ${VPATH}/Experiment.cpp ${VPATH}/Experiment.h
	${CC} ${IFLAGS} ${CFLAGS} -c ${VPATH}/Experiment.cpp

Trace.o: ${VPATH}/Trace.cpp ${VPATH}/Trace.h
	${CC} ${IFLAGS} ${CFLAGS} -c ${VPATH}/Trace.cpp

Interaction.o: ${VPATH}/Interaction.cpp ${VPATH}/Interaction.h
	${CC} ${IFLAGS} ${CFLAGS} -c ${VPATH}/Interaction.cpp

Molecule.o: ${VPATH}/Molecule.cpp ${VPATH}/Molecule.h
	${CC} ${IFLAGS} ${CFLAGS} -c ${VPATH}/Molecule.cpp

CustomMolecules.o: ${VPATH}/CustomMolecules.cpp ${VPATH}/CustomMolecules.h
	${CC} ${IFLAGS} ${CFLAGS} -c ${VPATH}/CustomMolecules.cpp

CustomInteractions.o: ${VPATH}/CustomInteractions.cpp ${VPATH}/CustomInteractions.h
	${CC} ${IFLAGS} ${CFLAGS} -c ${VPATH}/CustomInteractions.cpp

run: ${EXE_FILE}
	EvoDevo

leakcheck: ${EXE_FILE}
	valgrind --leak-check=full EvoDevo

clean:
	rm -rf ${OBJ_FILE} ${EXE_FILE} ${OUTPUT_DIR} 
