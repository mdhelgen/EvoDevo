CC	= g++
CFLAGS	= -g -O0 -Wall
IFLAGS = -I"../include"
LFLAGS = -L"../lib"


OBJ_FILE	= Main.o Cell.o Experiment.o DerivGraph.o 
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

DNA.o: ./data/DNA.cpp 
	${CC} ${IFLAGS} ${CFLAGS} -c ./data/DNA.cpp

Interaction.o: ./data/Interaction.cpp
	${CC} ${IFLAGS} ${CFLAGS} -c ./data/Interaction.cpp

MRNA.o: ./data/MRNA.cpp
	${CC} ${IFLAGS} ${CFLAGS} -c ./data/MRNA.cpp

Molecule.o: ./data/Molecule.cpp
	${CC} ${IFLAGS} ${CFLAGS} -c ./data/Molecule.cpp

Protein.o: ./data/Protein.cpp
	${CC} ${IFLAGS} ${CFLAGS} -c ./data/Protein.cpp

Experiment.o: Experiment.cpp Experiment.h
	${CC} ${IFLAGS} ${CFLAGS} -c Experiment.cpp

Term.o: ./org/Term.cpp
	${CC} ${IFLAGS} ${CFLAGS} -c ./org/Term.cpp

clean:
	rm -rf ${OBJ_FILE} ${EXE_FILE} ${OUTPUT_DIR} 
