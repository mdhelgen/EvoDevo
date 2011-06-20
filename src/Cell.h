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


private:
	static int CellID;
	int currentGen;
	MTRand r;
	DerivGraph* equations;

};





#endif
