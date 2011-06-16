/**
 * Experiment.h
 *
 * Experiment header file.
 *
 * Experiment contains a vector of Cell* pointers, as well as some information concerning generations and scoring.
 *
 */

#ifndef EXPERIMENT_H_
#define EXPERIMENT_H_

#include <stdlib.h>
#include <vector>
#include <fstream>
#include "Cell.h"

using namespace std;


class Experiment {
public:
	Experiment(int, int);
	~Experiment();

	void start();

private:
	vector<Cell*> cells;

	int maxGenerations;
	int numCells;
	
	int scoringInterval;
	int numHighScores;
};

#endif /* EXPERIMENT_H_ */
