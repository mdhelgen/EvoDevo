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
 * @param ncells number of Cell objects to be created.
 * @param generations number of Generations the Experiment will run for.
 * @param max_basic maximum number of basic proteins allowed in each Cell
 * @param max_ptm maximum number of PTM proteins allowed in each Cell
 * @param max_comp maximum number of complexed proteins allowed in each Cell
 * @param max_prom maximum number of protein-promoter interactions allowed in each cell
 * @param min_kinetic_rate the lower bound on randomly generated kinetic rates
 * @param max_kinetic_rate the upper bound on randomly generated kinetic rates
 * @param rk_time_limit the stopping condition for runge-kutta iteration
 * @param rk_time_step how much time to advance each iteration
 * @param initial_conc the initial concentration for molecules
 */
Experiment::Experiment(int ncells, int generations, int max_basic, int max_ptm, int max_comp, int max_prom, float min_kinetic_rate, float max_kinetic_rate, float rk_time_limit, float rk_time_step, float initial_conc)
	   :maxBasic(max_basic), maxPTM(max_ptm), maxComp(max_comp), maxProm(max_prom), minKineticRate(min_kinetic_rate), maxKineticRate(max_kinetic_rate), rkTimeLimit(rk_time_limit), rkTimeStep(rk_time_step), initialConc(initial_conc){

	
	t.trace("init","Creating new Experiment\n");
	
	t.trace("mloc","Experiment location at %p\n", this);

	t.trace("args","%d Cells\n",ncells);
	t.trace("args","%d Generations\n", generations);

	maxGenerations = generations;

	char buf[200];
	pid = getpid();
	
	prefix = "../output";
	mkdir(prefix , S_IRWXU | S_IRWXG | S_IRWXO);
	
	sprintf(buf, "%s/%d", prefix, pid);
	mkdir(buf, S_IRWXU | S_IRWXG | S_IRWXO);

	//create the cell objects and add them to our cells vector
	for (int i = 0; i < ncells; i++){
		t.trace("init","Creating Cell (%d)\n",i);
		cells.push_back(new Cell(maxBasic, maxPTM, maxComp, maxProm,minKineticRate,maxKineticRate, rkTimeStep, rkTimeLimit, initialConc));
		
		//set up directory for the new cell in output folder		
		sprintf(buf, "%s/%d/cell%d", prefix, pid, cells.back()->getID());
		mkdir(buf, S_IRWXU | S_IRWXG | S_IRWXO); 
		sprintf(buf, "%s/%d/cell%d/csv", prefix, pid, cells.back()->getID());
		mkdir(buf, S_IRWXU | S_IRWXG | S_IRWXO); 
	}

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
		t.trace("free","Deleting Cells[%d] at location %p\n",i,&(cells[i]));	
		delete cells[i];
	}

	t.trace("free","Deleting Cell[] object at location %d\n", &cells);
	cells.clear();



}

/**
 * Experiment::setOutputOptions(int, int, int, int, int, int)
 * 
 * Set options received from the command line parameters which deal with output files or formats
 * A value of 1 enables the option, and 0 disables it. 
 *
 * @param gv_flag Triggers GraphViz output of png files which display the directed graph representation of the cell 
 * @param gp_flag Triggers Gnuplot output of png files which display the concentration of molecules over time at a given generation
 * @param eachgen_flag Triggers output of files for every generation, not only generations on the scoring interval
 * @param csv_cell Triggers output of csv files which contain the interactions and rates within a cell
 * @param csv_data Triggers output of csv files which contain concentration data of molecules over time at a given generation
 * @param scoring_interval How many generations should occur between scoring and output
 *
 */
void Experiment::setOutputOptions(int gv_flag, int gp_flag, int eachgen_flag, int csv_cell, int csv_data, int scoring_interval){

	graphviz_enabled = gv_flag;
	gnuplot_enabled = gp_flag;
	output_each_gen = eachgen_flag;
        output_csv_interactions = csv_cell;
	output_csv_data = csv_data;	
	scoringInterval = scoring_interval;
	
}
/**
 * Experiment::start()
 *
 * The main experiment loop for the simulation. For each generation, every cell is mutated. Depending on the current generation number, and the value
 * of scoringInterval, The cells may also be evaluated by runge-kutta, and output the files related to the best cell.
 *
 * The experiment terminates once the generations reach the generation limit.
*/
void Experiment::start()
{

	t.trace("args","Graphviz: %d\n",graphviz_enabled);
	t.trace("args","Gnuplot: %d\n", gnuplot_enabled);

	int bestScore = -1;
	Cell* bestCell = 0;

	//generational loop
	for(int i = 1; i <= maxGenerations; i++)
	{
		bestScore = -1;
		bestCell = 0;
		
		t.trace("gens","Generation %d started (max %d)\n",i, maxGenerations);
		for(unsigned int c = 0; c < cells.size(); c++)
		{
			t.trace("mutate","Gen %-3d Cell loc %p\n", i, cells[c]);
			//mutate
			cells[c]->mutate();
			
			//find the best cell
			//if scoring interval is 5, this runs every 5 generations
			if(i % scoringInterval == 0){
				cells[c]->rk();
				if(cells[c]->getScore() < -1  ){
					cells[c]->outputDataPlot(prefix, pid);
					cells[c]->outputDotImage(prefix, pid);
				}	

				//keep track of the cell with the highest score so far
				if(cells[c]->getScore() > bestScore){
					bestCell = cells[c];
					bestScore = bestCell->getScore();
					t.trace("score","Best cell is cell %d with score %d\n",bestCell->getID(), bestCell->getScore());
				}		
			}	

			//if the flag is set, generate output for every cell during every generation
			//this will significantly increase the runtime of the simulation
			if(output_each_gen){
				//if the flag is set, generate output every generation
				cells[c]->rk();

				cells[c]->stochasticSim();

				if(graphviz_enabled)
					cells[c]->outputDotImage(prefix, pid);
				if(gnuplot_enabled)
					cells[c]->outputDataPlot(prefix, pid);
				if(output_csv_data)
					cells[c]->outputDataCsv(prefix, pid);
				if(output_csv_interactions)
					cells[c]->outputInteractionCsv(prefix, pid);
			}
		}
		//if the scoring interval is 5, this runs every 5 generations
		if(i % scoringInterval == 0){
			
			//all cells have been checked, so the bestCell variable holds the cell with the highest score
			t.trace("score","Best cell at end of Generation %d is cell %d with score %d\n", i, bestCell->getID(), bestCell->getScore());
			
			//output the best cell
			if(graphviz_enabled)
				bestCell->outputDotImage(prefix, pid);
			if(gnuplot_enabled)
				bestCell->outputDataPlot(prefix, pid);
			if(output_csv_data)	
				bestCell->outputDataCsv(prefix, pid);
			if(output_csv_interactions)
				bestCell->outputInteractionCsv(prefix, pid);
		}
		t.trace("gens","Generation %d finished (max %d)\n",i, maxGenerations);
	}
	
	return;

	//this can probably be removed

	//generate output at the end of the experiment
	for(unsigned int c = 0; c < cells.size(); c++){
		if(graphviz_enabled)
			cells[c]->outputDotImage(prefix, pid);
		if(gnuplot_enabled)
			cells[c]->outputDataPlot(prefix, pid);
	}

}

