#include <iostream>
#include "Cell.h"

using namespace std;

#include "ExternTrace.h"

Cell::Cell(){

    t->trace("init", "New Cell Created\n");
    equations = new DerivGraph();

}


Cell::~Cell(){
  
    t->trace("free", "Deleting DerivGraph object at location %d\n", &equations);    
    delete equations;

}
