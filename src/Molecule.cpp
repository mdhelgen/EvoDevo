#include "Molecule.h"

#include "ExternTrace.h"


Molecule::Molecule(){

	t.trace("init", "Creating new Molecule\n");
	t.trace("mloc", "Molecule location at %u\n", (unsigned int) this);

	currentConcentration = 4;
	rkVal[0] = 0;
	rkVal[1] = 0;
	rkVal[2] = 0;
	rkVal[3] = 0;
	
	t.trace("init", "New Molecule created\n");

}

Molecule::~Molecule(){

}

float Molecule::getrkVal(int k){
	return rkVal[k];
}
float Molecule::getValue(){

	return currentConcentration;

}

void Molecule::setValue(float v){

	initialConcentration = v;
	currentConcentration = v;
	rungeKuttaSolution.push_back(v);
}

void Molecule::updateRkVal(int index, float amount){

	rkVal[index] += amount;

}

void Molecule::nextPoint(float step){

	t.trace("rk-4","rkvals: %f %f %f %f\n",rkVal[0], rkVal[1], rkVal[2], rkVal[3]);
	currentConcentration += ((step/6) * (rkVal[0] + 2*rkVal[1] + 2*rkVal[2] + rkVal[3]));

	rungeKuttaSolution.push_back(currentConcentration);
	rkVal[0] = 0;
	rkVal[1] = 0;
	rkVal[2] = 0;
	rkVal[3] = 0;

}

void Molecule::outputRK(){
	for(unsigned int i = 0; i< rungeKuttaSolution.size(); i++)
		t.trace("rk-4","%f - %u - %f\n", initialConcentration, i, rungeKuttaSolution[i]);
}
