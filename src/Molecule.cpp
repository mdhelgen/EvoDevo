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
	
	wasPTM = 0;

	//initialize Runge-Kutta intermediate values
	rkVal[0] = 0;
	rkVal[1] = 0;
	rkVal[2] = 0;
	rkVal[3] = 0;
	
	t.trace("init", "New Molecule created\n");
	
	PTMArray[0] = 0;
	PTMArray[1] = 0;
	PTMArray[2] = 0;
	PTMArray[3] = 0;
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
	t.trace("rk-val","%s rkval[%d] update %f + %f = %f\n",getShortName(), index, rkVal[index], amount, rkVal[index]+amount);	
	rkVal[index] += amount;
	t.trace("rk-val","%s new rkval[%d] = %f\n", getShortName(), index, rkVal[index]);
}

/**
 * float Molecule::rkApprox(int, float)
 *
 */
float Molecule::rkApprox(int rkIteration, float rkStepSize){
	
	float approxVal = 0;
	
	switch(rkIteration){
	case 0:
		approxVal = getValue();
		break;
	case 1:
		approxVal = (getValue() + ( rkVal[0] * (rkStepSize/2)));
		break;
	case 2:
		approxVal = (getValue() + ( rkVal[1] * (rkStepSize/2)));
		break;
	case 3:
		approxVal = (getValue() + ( rkVal[2] * rkStepSize ));
		break;
	}
	return approxVal < 0 ? 0 : approxVal;
}

/**
 * void Molecule::nextPoint(float)
 *
 */
void Molecule::nextPoint(float step){

	t.trace("rk-val","%s rkvals: %f %f %f %f\n",getShortName(), rkVal[0], rkVal[1], rkVal[2], rkVal[3]);
	
	float delta = ((step/6) * (rkVal[0] + 2*rkVal[1] + 2*rkVal[2] + rkVal[3]));
	t.trace("rk-new","%s(%d) conc: %f delta: %f\n",getShortName(),rungeKuttaSolution.size(), currentConcentration, delta);
	currentConcentration+=delta;	
	
	if(currentConcentration < 0){
		t.trace("rk-new", "%s %f being set to 0\n",getShortName(),currentConcentration);
		currentConcentration = 0;
	}
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

vector<float>* Molecule::getRungeKuttaSolution(){
	return &rungeKuttaSolution;

}

int Molecule::getPTMCount(int index){
	return PTMArray[index];
}

/*
int Molecule::getScore(){



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



