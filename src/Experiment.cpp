/**
 * Experiment.cpp
 *
 * Create cells, loop generations, delete cells.
 *
 *
 */

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <vector>

#include "Experiment.h"

using namespace std;

//external declaration of Trace t
#include "ExternTrace.h"


/**
 * Experiment::Experiment(int, int)
 *
 * Experiment constructor.
 *
 * @param number of Cell objects to be created.
 * @param number of Generations the Experiment will run for.
 *
 */
Experiment::Experiment(int ncells, int generations) {


	
	t.trace("init","Creating new Experiment\n");
	
	t.trace("mloc","Experiment location at %u\n",(unsigned int) this);

	t.trace("args","%d Cells\n",ncells);
	t.trace("args","%d Generations\n", generations);

	maxGenerations = generations;

	//create the cell objects and add them to our cells vector
	for (int i = 0; i < ncells; i++){
		t.trace("init","Creating Cell (%d)\n",i);
		cells.push_back(new Cell());
	}
/*	
	char buf[200];
	
	const char* prefix = "../output";
	mkdir(prefix , S_IRWXU | S_IRWXG | S_IRWXO);
	
	sprintf(buf, "%s/%d", prefix, getpid());
	mkdir(buf, S_IRWXU | S_IRWXG | S_IRWXO);

	sprintf(buf, "%s/%d/cell", prefix, getpid());
	mkdir(buf, S_IRWXU | S_IRWXG | S_IRWXO);

	sprintf(buf, "%s/%d/trace.txt", prefix, getpid());
	FILE* tracef = fopen(buf, "a+");
	fprintf(tracef, "test\n");

	t.setTraceFile(tracef);

	sprintf(buf, "%s/%d/cell/genXstd.csv", prefix, getpid());
	FILE* stdfile = fopen(buf, "a+");
	fprintf(stdfile, "test\n");

	sprintf(buf, "%s/%d/cell/genXprec.csv", prefix, getpid());
	FILE* precfile = fopen(buf, "a+");
	fprintf(precfile, "test\n");
*/

	t.trace("init","New Experiment created\n");
}

/**
 * Experiment::~Experiment(int, int)
 *
 * Experiment destructor.
 *
 * Deletes the Cell objects from the cells vector, then deletes the vector itself.
 */
Experiment::~Experiment() {


	//call the destructor for each cell in the vector
	for(unsigned i = 0; i < cells.size(); i++){
		t.trace("free","Deleting Cells[%d] at location %d\n",i,&(cells[i]));	
		delete cells[i];
	}

	t.trace("free","Deleting Cell[] object at location %d\n", &cells);
	cells.clear();



}

//void Experiment::outputGeneration(){


//}


/**
 * Deprecated 
 */
void Experiment::start()
{
//getscore before anything
for(int i = 1; i <= maxGenerations; i++)
{
	for(unsigned int c = 0; c < cells.size(); c++)
	{
		t.trace("mutate","Gen %-3d Cell loc %u\n", i, (unsigned int) cells[c]);
		//mutate
		cells[c]->mutate();
		//getscore
	}
}

for(unsigned int c = 0; c < cells.size(); c++)
	cells[c]->outputDotImage();

/*
	ofstream output;


	char filename[181];
	snprintf(filename, 179, "./output/%d/score.log",getpid());

	//Comment this and document in the report.
	output.open(filename);

	// loop once for each generation
	for (int generation = 1; generation <= maxGenerations; generation++)
	{
		Cell* bestCell;
		int bestCellScore = -1;

		TRACE(4,"Beginning Generation " << generation << "...");
		vector<Cell*>::iterator cellIterator = cells.begin();
		while( cellIterator != cells.end() ) {


			(*cellIterator)->nextGeneration();
			++cellIterator;
		}

		// if this is a scoring generation
		if (generation % scoringInterval == 0){

			TRACE(5,"Beginning RK scoring...")
			// loop through all of the cells
			for(size_t i = 0; i < cells.size(); i++){

				cells[i]->rungeKuttaEvaluate(RUNGE_KUTTA_STD_STEP,RUNGE_KUTTA_STD_LIMIT);

				//get the score of the current cell
				int score = cells[i]->getScore(1);

				//if this is better than the previous best cell, save as the best cell
				if(score > bestCellScore){
					bestCell = cells[i];
					bestCellScore = score;

					TRACE(3,"Generation " << generation << " Best cell is now " << bestCell->cellID << " (score=" <<bestCellScore << ").")
				}
			}

			TRACE(4,"RK scoring complete.")

			TRACE(3,"Generation " << generation << " scoring finished, best cell is " << bestCell->cellID <<" (score=" <<bestCellScore <<").")

			output << "Generation " << generation << ", best cell is " << bestCell->cellID <<" (score=" <<bestCellScore <<")."<< endl;

			//redo Runge-Kutta with finer precision
			bestCell->rungeKuttaEvaluate(RUNGE_KUTTA_PRECISE_STEP,RUNGE_KUTTA_PRECISE_LIMIT);
			int betterScore = bestCell->getScore(1);
			TRACE(3, "Generation " << generation << " best cell score is "<<betterScore <<" after precise rk.")


		}

		TRACE(3,"Generation " << generation << " Finished.");

	}

	output.close();

	TRACE(2,"Experiment Finished.")
*/
}

