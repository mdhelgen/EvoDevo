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

int Cell::CellID = 0;

/**
 * Cell::Cell()
 *
 * Cell default constructor.
 *
 * Allocates:
 * 	1 derivGraph object
 */
Cell::Cell(){

    t.trace("init", "Creating new Cell\n");
    t.trace("mloc", "Cell location at %u\n", (unsigned int) this);
    equations = new DerivGraph();
   
    currentGen = 0;
   
    CellID++;

    equations->test();

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

    t.trace("free", "Deleting DerivGraph object at location %d\n", &equations);    
    delete equations;

}

int Cell::mutate(){

	double mutationCategory = r.rand(1);
	if(mutationCategory < .4)
		t.trace("mutate","Mutation Type: Small\n");
	else if(mutationCategory < .7)
		t.trace("mutate","Mutation Type: Large\n");
	else
		t.trace("mutate","Mutation Type: Null\n"); 
	
	
	equations->rungeKuttaEvaluate(1.0);
	return -1;
}


