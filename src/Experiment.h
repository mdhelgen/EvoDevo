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

#include <sys/stat.h>
#include <sys/types.h>
#include <cstdio>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include "Cell.h"

using namespace std;


class Experiment {
public:
	Experiment(int ncells, int generations, int max_basic, int max_ptm, int max_comp, int max_prom, float min_kinetic_rate, float max_kinetic_rate, float rk_time_limit, float rk_time_step, float initial_conc);
	~Experiment();

	void start();
	
	//set commandline output options
	void setOutputOptions(int, int, int, int, int, int);
private:
	vector<Cell*> cells;

	// default cell properties
	int maxGenerations;
	int numCells;
	int maxBasic;
	int maxPTM;
	int maxComp;
	int maxProm;
	
	// default interaction properties
	float minKineticRate;
	float maxKineticRate;

	// runge kutta options
	float rkTimeLimit;
	float rkTimeStep;
	int scoringInterval;

	// default molecule properties
	float initialConc;

	// output directory prefix (prefix/pid/<outputgoeshere>)
	const char* prefix;
	int pid;

	// output flags
	int graphviz_enabled;
	int gnuplot_enabled;
        int output_each_gen;
	int output_csv_interactions;
	int output_csv_data;

	// unused ?
	int numHighScores;
};

#endif /* EXPERIMENT_H_ */
