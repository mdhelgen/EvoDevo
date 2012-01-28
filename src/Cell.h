/**
 * Cell header file.
 * 
 * The Cell class mainly interfaces with the DerivGraph class which it contains.
 */

#ifndef CELL_H_
#define CELL_H_

#include "MersenneTwister.h"

#include "DerivGraph.h"

class Cell{

public:
	Cell(int, int, int, int,float,float, float, float, float);
	~Cell();
	int mutate();
	
	//output functions
	void outputDotImage(const char*, int);
	void outputDataPlot(const char*, int);
	void outputDataCsv(const char*, int);
	void outputInteractionCsv(const char*, int);
	
	// runge kutta functions
	void rk();
	void stochasticSim();
	int getScore();
	
	int getID(){ return CellID; };
private:
	
	// cell properties
	static int CellCounter;
	int CellID;
	int currentGen;

	// random generator
	MTRand r;

	// molecules/interactions class
	DerivGraph * equations;

	// runge kutta values
	float rkTimeStep;
	float rkTimeLimit;
};

#endif
