/**
 * Cell header file.
 *
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
	void outputDotImage();
	void outputDataPlot();
	int getScore();
	int getID(){ return CellID; };
	void rkTest();
	void rk();
private:
	static int CellCounter;

	int CellID;
	int currentGen;
	MTRand r;
	DerivGraph * equations;
	float rkTimeStep;
	float rkTimeLimit;
};





#endif
