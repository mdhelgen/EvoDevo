#include <iostream>
#include "Cell.h"

using namespace std;

#include "ExternTrace.h"

 
Cell::Cell(){

    t.trace("init", "Creating new Cell\n");
    t.trace("mloc", "Cell location at %u\n", (unsigned int) this);
    equations = new DerivGraph();
    equations->test();

    t.trace("init", "New Cell created\n");


}


Cell::~Cell(){

    t.trace("free", "Deleting DerivGraph object at location %d\n", &equations);    
    delete equations;



}
