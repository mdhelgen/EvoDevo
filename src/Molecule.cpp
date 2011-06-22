#include "Molecule.h"

#include "ExternTrace.h"

/**
 * Molecule::Molecule()
 *
 * Default Constructor
 *
 */
Molecule::Molecule(){

	t.trace("init", "Creating new Molecule\n");
	t.trace("mloc", "Molecule location at %u\n", (unsigned int) this);

	//set default concentration
	currentConcentration = 4;
	initialConcentration = 4;

	longName = "Molecule";
	shortName = "a";
	
	moleculeID = -1;

	//initialize Runge-Kutta intermediate values
	rkVal[0] = 0;
	rkVal[1] = 0;
	rkVal[2] = 0;
	rkVal[3] = 0;
	
	t.trace("init", "New Molecule created\n");

}

/**
 * Molecule::~Molecule()
 *
 * Default Destructor
 */
Molecule::~Molecule(){

}

/**
 * float Molecule::getrkVal(int)
 *
 * Get the value of the intermediate Runge-Kutta value for a particular iteration.
 *
 * @param k Which iterations rkVal to return
 *
 * @return The value of this molecules rkVal[k]
 *
 */
float Molecule::getrkVal(int k){
	if(k >= 0 && k <= 3)
		return rkVal[k];
	return 0;
}

/**
 * float Molecule::getValue()
 *
 */
float Molecule::getValue(){

	return currentConcentration;

}

/**
 * void Molecule::setValue(float)
 *
 */
void Molecule::setValue(float v){

	initialConcentration = v;
	currentConcentration = v;
	rungeKuttaSolution.push_back(v);
}

/**
 * void Molecule::updateRkVal(int, float)
 *
 */
void Molecule::updateRkVal(int index, float amount){
	
	rkVal[index] += amount;
}

/**
 * float Molecule::rkApprox(int, float)
 *
 */
float Molecule::rkApprox(int rkIteration, float rkStepSize){
	switch(rkIteration){
	case 0:
		return getValue();
	case 1:
		return (getValue() + ( rkVal[0] * (rkStepSize/2)));
	case 2:
		return (getValue() + ( rkVal[1] * (rkStepSize/2)));
	case 3:
		return (getValue() + ( rkVal[2] * rkStepSize ));
	}

	return currentConcentration;

}

/**
 * void Molecule::nextPoint(float)
 *
 */
void Molecule::nextPoint(float step){

	t.trace("rk-4","rkvals: %f %f %f %f\n",rkVal[0], rkVal[1], rkVal[2], rkVal[3]);
	currentConcentration += ((step/6) * (rkVal[0] + 2*rkVal[1] + 2*rkVal[2] + rkVal[3]));

	rungeKuttaSolution.push_back(currentConcentration);
	rkVal[0] = 0;
	rkVal[1] = 0;
	rkVal[2] = 0;
	rkVal[3] = 0;

}

char* Molecule::getShortName(){
	
	memset(buf, '\0', 80);
	sprintf(buf, "%s%d", shortName, moleculeID);
	return buf;
}

char* Molecule::getLongName(){

	memset(buf, '\0', 80);
	sprintf(buf, "%s %d", longName, moleculeID);
	return buf;
}
void Molecule::setID(int i){

	moleculeID = i;
}

int Molecule::getID(){
	return moleculeID;
}
void Molecule::reset(){
	
	rungeKuttaSolution.erase(rungeKuttaSolution.begin(), rungeKuttaSolution.end());
	rungeKuttaSolution.push_back(initialConcentration);
	currentConcentration = initialConcentration;
	rkVal[0] = 0;
	rkVal[1] = 0;
	rkVal[2] = 0;
	rkVal[3] = 0;
}


/**
 * void Molecule::outputRK()
 *
 * TEST METHOD
 */
void Molecule::outputRK(){
	for(unsigned int i = 0; i< rungeKuttaSolution.size(); i++)
		t.trace("rk-4","%s - %u - %f\n", getShortName(), i, rungeKuttaSolution[i]);
}
