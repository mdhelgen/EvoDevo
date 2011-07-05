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
	Cell();
	~Cell();
	int mutate();
	void outputDotImage();
	void outputDataPlot();
	void getScore();
private:
	static int CellCounter;

	int CellID;
	int currentGen;
	MTRand r;
	DerivGraph* equations;

};





#endif
