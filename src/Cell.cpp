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
