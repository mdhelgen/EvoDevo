#include "Molecule.h"

#include "ExternTrace.h"


Molecule::Molecule(){

	t.trace("init", "Creating new Molecule\n");
	t.trace("mloc", "Molecule location at %u\n", (unsigned int) this);

	currentConcentration = 4;

	t.trace("init", "New Molecule created\n");

}

Molecule::~Molecule(){

}


float Molecule::getValue(){

	return currentConcentration;

}
