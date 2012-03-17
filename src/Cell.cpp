/**
 * Cell implementation file.
 *
 */
#include <iostream>
#include "Cell.h"
using namespace std;


//external declaration of Trace t
#include "ExternTrace.h"

int Cell::CellCounter = 0;

/**
 * Cell::Cell()
 *
 * Cell default constructor.
 *
 * Allocates:
 * 	1 derivGraph object
 */
Cell::Cell(int max_basic, int max_ptm, int max_comp, int max_promoter,float min_kinetic_rate, float max_kinetic_rate, float rk_time_step, float rk_time_limit, float initial_conc){

    t.trace("init", "Creating new Cell\n");
    t.trace("mloc", "Cell location at %p\n", this);
    equations = new DerivGraph();
    
    //set the limits of mutation occurrences
    equations->setLimits(max_basic, max_ptm, max_comp, max_promoter);
    
    //kinetic rate lower and upper bounds
    equations->setKineticRateLimits(min_kinetic_rate, max_kinetic_rate);

    //runge kutta timestep and time limit
    equations->setRungeKuttaEval(rk_time_step, rk_time_limit);

    //molecule initial concentration
    equations->setDefaultInitialConc(initial_conc);

    rkTimeStep = rk_time_step;
    rkTimeLimit = rk_time_limit;

    currentGen = 0;
  
    //assign the cell a unique number 
    CellID = CellCounter++;


    //initialize the cell with a single basic protein
    equations->newBasic();

    t.trace("init", "New Cell created\n");

}

/**
 * Cell::~Cell()
 *
 * Cell default destructor.
 *
 * Frees:
 * 	1 derivGraph object
 *
 */
Cell::~Cell(){

t.trace("free","Deleting DerivGraph object at %p\n", equations);
delete equations;

}

/**
 * Cell::mutate()
 * 
 * Choose a random mutation to apply to the cell.
 *  A mutation category is chosen first, and a random mutation is then selected from within that category.
 *    Mutations and their effects can be added fairly easily by modifying the probabilities, and coding
 *    the effect in the DerivGraph class.
 *
 * (Small)
 *   Forward Rate Change
 *   Reverse Rate Change
 *   Degradation Rate Change
 *   New PTM
 *   Histone Modification
 * (Large)
 *   New Complex
 *   New Basic Protein
 *   New Protein-Promoter
 * (Null Mutation) -- no op
 *
 */
int Cell::mutate(){
	//advance generation
	currentGen++;
	
	//random values for mutation selection
	double mutationCategory = r.rand(1);
	double mutationType = r.rand(1);

	//small mutation category
	if(mutationCategory < .4)
	{
		t.trace("mutate","Mutation Category: Small\n");
	
		//forward rate change	
		if(mutationType < .2)
		{
			t.trace("mutate","Mutation Type: Forward Rate Change\n");	
			equations->forwardRateChange();
		}
		//reverse rate change
		else if(mutationType < .4)
		{
			t.trace("mutate","Mutation Type: Reverse Rate Change\n");	
			equations->reverseRateChange();
		}
		//degradation rate change
		else if(mutationType < .6)
		{
			t.trace("mutate","Mutation Type: Degradation Rate Change\n");	
			equations->degradationRateChange();
		}
		//new Post Translational Modification
		else if(mutationType < .8)
		{
			t.trace("mutate","Mutation Type: New PTM\n");	
			equations->newPTM();
		}
		//histone modification
		else
		{
			t.trace("mutate","Mutation Type: Histone Modification\n");	
			equations->histoneMod();
		}
	}
	//large mutation category
	else if(mutationCategory < .7)
	{
		t.trace("mutate","Mutation Category: Large\n");
	 	//new complex	
		if(mutationType < .33)
		{
			t.trace("mutate","Mutation Type: New Protein-Protein Complex\n");	
			equations->newComplex();
		}	
		//new basic protein
		else if(mutationType < .67)
		{
			t.trace("mutate","Mutation Type: New Basic Protein\n");	
			equations->newBasic();
		}
		//new protein-promoter
		else
		{
			t.trace("mutate","Mutation Type: New Protein-Promoter Interaction\n");
			equations->newPromoter();
		}
	}
	//null mutation
	else
		t.trace("mutate","Mutation Category: Null\n"); 
	
	
	return -1;
}
/**
 * int Cell::getScore()
 * 
 * Get the "score" associated with the cell, for ranking purposes. Each scoring generation can optionally only display the "best" cell to cut
 * down on generated data.
 *
 * The default behavior or getScore is to assign the score of the highest scored molecule within the cell to the cell itself.
 *
 * @return the score of the cell
 */
int Cell::getScore(){

	//get highest scored molecule within the cell
	Molecule* m = equations->getBestMolecule(CellID);
	
	// return the score of the best molecule
	return m->getScore();
		
}
/**
 * void Cell::outputDotImage(const char*, int)
 *
 * Output a png image using graphviz of the current cell's directed graph representation
 *
 * @param prefix the prefix of the output folder, relative to the execution directory (typically "../output")
 * @param pid a relatively unique value for the output directory (the pid of the current process)
 */
void Cell::outputDotImage(const char* prefix, int pid){
	equations->outputDotImage(prefix, pid, CellID, currentGen);
}

/**
 * void Cell::outputDotImage(const char*, int)
 *
 * Output png images of the concentration vs time for each molecule in the cell at current generation
 * 
 *   Note: Cell::rk() should be called during the current generation prior to this call to re-solve the equations
 *
 * @param prefix the prefix of the output folder, relative to the execution directory (typically "../output")
 * @param pid a relatively unique value for the output directory (the pid of the current process)
 */
void Cell::outputDataPlot(const char* prefix, int pid){
	equations->outputDataPlot(prefix, pid, CellID, currentGen, rkTimeStep);
}

/**
 * void Cell::outputDataCsv(const char*, int)
 *
 * Output csv files of the concentration vs time for each molecule in the cell at current generation
 *
 *   Note: Cell::rk() should be called during the current generation prior to this call to re-solve the equations
 *
 * @param prefix the prefix of the output folder, relative to the execution directory (typically "../output")
 * @param pid a relatively unique value for the output directory (the pid of the current process)
 */
void Cell::outputDataCsv(const char* prefix, int pid){
	equations->outputDataCsv(prefix, pid, CellID, currentGen, rkTimeStep);
}

/**
 * void Cell::outputDotImage(const char*, int)
 *
 * Output a csv file containing the interactions and corresponding rates within the current cell at the current generation
 *
 * @param prefix the prefix of the output folder, relative to the execution directory (typically "../output")
 * @param pid a relatively unique value for the output directory (the pid of the current process)
 */
void Cell::outputInteractionCsv(const char* prefix, int pid){
	equations->outputInteractionCsv(prefix, pid, CellID, currentGen);
}

/**
 * void Cell::rk()
 *
 * Solve the equations based on the interactions and molecules within the cell by numerical approximation, using
 * the runge-kutta 4th order method.
 *
 * This method is computationally intensive.
 */
void Cell::rk(){
	equations->rungeKuttaEvaluate(rkTimeStep, rkTimeLimit);
}

/**
 * void Cell::stochasticSim()
 *
 * Simulate the cell behavior based on a stochastic model (the gillespie algorithm), as opposed to the deterministic model offered by rk().
 *
 * 
 *
 *
 * This method is computationally intensive
 */
void Cell::stochasticSim(){

	equations->gillespieEvaluate();

}
