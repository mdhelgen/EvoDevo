/**
 * Cell implementation file.
 *
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
    t.trace("mloc", "Cell location at %u\n", (unsigned int) this);
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
   
    CellID = CellCounter++;

    //equations->>test();

	equations->newBasic();
//	equations->newPromoter();
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

t.trace("free","Deleting DerivGraph object at %u\n",(unsigned int) equations);
delete equations;

}

int Cell::mutate(){
	
	currentGen++;
	double mutationCategory = r.rand(1);
	double mutationType = r.rand(1);

	if(mutationCategory < .4)
	{
		t.trace("mutate","Mutation Category: Small\n");
		
		if(mutationType < .2)
		{
			t.trace("mutate","Mutation Type: Forward Rate Change\n");	
			equations->forwardRateChange();

		}//end fwd rate change
		else if(mutationType < .4)
		{
			t.trace("mutate","Mutation Type: Reverse Rate Change\n");	
			equations->reverseRateChange();

		}//end rev rate change
		else if(mutationType < .6)
		{
			t.trace("mutate","Mutation Type: Degradation Rate Change\n");	
			equations->degradationRateChange();
		
		}//end deg rate change
		else if(mutationType < .8)
		{
			t.trace("mutate","Mutation Type: New PTM\n");	
			equations->newPTM();
		}//end new ptm
		else
		{
			t.trace("mutate","Mutation Type: Histone Modification\n");	
			equations->histoneMod();
		}//end histone mod

	}//end small category
	else if(mutationCategory < .7)
	{
		t.trace("mutate","Mutation Category: Large\n");
		
		if(mutationType < .33)
		{
			t.trace("mutate","Mutation Type: New Protein-Protein Complex\n");	
			equations->newComplex();
		}//end new complex
		else if(mutationType < .67)
		{
			t.trace("mutate","Mutation Type: New Basic Protein\n");	
			equations->newBasic();
				
		}//end new basic
		else
		{
			t.trace("mutate","Mutation Type: New Protein-Promoter Interaction\n");
			equations->newPromoter();
		}//end new promoter

	}//end large category
	else
	{
		t.trace("mutate","Mutation Category: Null\n"); 



	}//end null mutation
	
	//equations->rungeKuttaEvaluate(rkTimeStep, rkTimeLimit);
//	equations->outputDotImage(CellID,currentGen );
	
	return -1;
}

int Cell::getScore(){

	Molecule* m = equations->getBestMolecule(CellID);
	
	return m->getScore();
	
		
}

void Cell::outputDotImage(){
//	printf("dot\n");
	equations->outputDotImage(CellID, currentGen);
}

void Cell::outputDataPlot(){
	equations->outputDataPlot(CellID, currentGen, rkTimeStep);
}


void Cell::rkTest(){

	float steps[] = {1, .75, .5, .25, .1, .05, .025, .01, .001};
	for(int i = 0; i < 9 ; i++){
		equations->setRungeKuttaEval(steps[i], 20);
		equations->rungeKuttaEvaluate(steps[i], 20);
		equations->outputDataPlot(i, currentGen, steps[i]);
	}


}

void Cell::rk(){
	equations->rungeKuttaEvaluate(rkTimeStep, rkTimeLimit);
}

