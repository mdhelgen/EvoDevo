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
	t.trace("mloc", "Molecule location at %p\n", this);

	//set default concentration
	currentConcentration = 5;
	initialConcentration = 5;

	//set the long and short name prefix for this molecule
	longName = "Molecule";
	shortName = "a";
	
	//set the molecule id
	moleculeID = -1;

	numChanges = 0;
	prevDir = 0;
	currentDir = 0;


	wasPTM = 0;

	//initialize Runge-Kutta intermediate values
	rkVal[0] = 0;
	rkVal[1] = 0;
	rkVal[2] = 0;
	rkVal[3] = 0;
	
	t.trace("init", "New Molecule created\n");
	
	//initialize the PTMArray
	PTMArray[0] = 0;
	PTMArray[1] = 0;
	PTMArray[2] = 0;
	PTMArray[3] = 0;


	//gillespie_numMolecules = 50;
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
 * (Virtual Function) 
 *
 * Get the current value of this molecule. 
 *
 * @returns the current value of the concentration
 */
float Molecule::getValue(){

	return currentConcentration;

}

/**
 * void Molecule::setValue(float)
 *
 * @param v the new value to set as the concentration
 */
void Molecule::setValue(float v){

	initialConcentration = v;
	currentConcentration = v;
	rungeKuttaSolution.push_back(v);
}

/**
 * void Molecule::updateRkVal(int, float)
 *
 * Adds some amount to the specified intermediate values used by Runge-Kutta.
 *
 * @param index The index of the rkValue array to update
 * @param amount The amount to add to the rkValue array
 */
void Molecule::updateRkVal(int index, float amount){
	t.trace("rk-val","%s rkval[%d] update %f + %f = %f\n",getShortName(), index, rkVal[index], amount, rkVal[index]+amount);	
	rkVal[index] += amount;
	t.trace("rk-val","%s new rkval[%d] = %f\n", getShortName(), index, rkVal[index]);
}

/**
 * float Molecule::rkApprox(int, float)
 * (Virtual Function) 
 *
 * Returns the next approximate value of this molecule for the next timestep for the specified
 * stage of Runge-Kutta. Runge-Kutta uses successive iterations to make more accurate approximations
 * of a solution. 
 *
 * rkApprox should be used in Interaction::getEffect, to provide the Runge-Kutta corrected concentrations
 * of molecules during runge-kutta calculation instead of the base value for all iterations.
 *
 * @param rkIteration the current iteration of Runge-Kutta
 * @param rkStepSize the timestep being used by Runge-Kutta
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
	return (approxVal <= 0 ? 0 : approxVal);
}

/**
 * void Molecule::nextPoint(float)
 * 
 * Adds a data point to the rungeKuttaSolution based on the rkVals calculated by Runge-Kutta
 *
 * @param step The stepsize used to calculate the rkVals
 */
void Molecule::nextPoint(float step){

	float oldConc = currentConcentration;

	t.trace("rk-val","%s rkvals: %f %f %f %f\n",getShortName(), rkVal[0], rkVal[1], rkVal[2], rkVal[3]);
	
	//runge-kutta calculation of change in value during the current timestep
	float delta = ((step/6) * (rkVal[0] + 2*rkVal[1] + 2*rkVal[2] + rkVal[3]));
	
	t.trace("rk-new","%s(%d) conc: %f delta: %f\n",getShortName(),rungeKuttaSolution.size(), currentConcentration, delta);
	
	//add the change in value to the current value
	currentConcentration+=delta;	

	//ensure non-negative concentration
	if(currentConcentration < 0){
		t.trace("rk-new", "%s %f being set to 0\n",getShortName(),currentConcentration);
		currentConcentration = 0;
	}

	//add the new value to the solution vector
	rungeKuttaSolution.push_back(currentConcentration);
	
	//reset the rkVals to zero
	rkVal[0] = 0;
	rkVal[1] = 0;
	rkVal[2] = 0;
	rkVal[3] = 0;
	
	/*
	 * Scoring -- Oscillation counting
	 * 
	 * as points are added, check to see if the direction has changed from the previous direction
	 */
	
	float actualChange = currentConcentration-oldConc;
	
	if(actualChange <= .0001 && actualChange >= -0.0001)
		return;

	//the value just increased
	if(actualChange > 0)
		currentDir = 1;
	//the value just decreased
	if(actualChange < 0)
		currentDir = -1;
	
	t.trace("score", "%s%d - (%f , %f), dir = %d, prev = %d\n",shortName, moleculeID,currentConcentration,delta  , currentDir, prevDir);

	//if the value previously decreased and just increased, the last point was a minimum
	if(prevDir == -1 && currentDir == 1){ 
		minima.push_back(oldConc);
		numChanges++;
	}
	//if the value previously increased and just decreased, the last point was a maximum
	if(prevDir == 1 && currentDir == -1){
		maxima.push_back(oldConc);
		numChanges++;
        }

	//save for the next point
	prevDir = currentDir;

	return;
}
/**
 * char* Molecule::getShortName()
 * (Virtual function)
 *
 * Return the "short" name of a molecule.
 *
 * The short name consists of the short prefix set in the constructor appended to the moleculeID
 * with no space in between.
 * Ex. g1, p4, ptm2
 *
 * @return the short name of the current molecule
 *
 */
char* Molecule::getShortName(){
	
	memset(buf, '\0', 80);
	sprintf(buf, "%s%d", shortName, moleculeID);
	return buf;
}

/**
 * char* Molecule::getLongName()
 * (Virtual function)
 * 
 * Return the "long" name of a molecule.
 *
 * The long name consists of the long prefix set in the constructor appended to the moleculeID with a space
 * in between.
 * Ex. DNA 1, Protein 3, Complex 8
 *
 * @return the long name of the current molecule
 */
char* Molecule::getLongName(){

	memset(buf, '\0', 80);
	sprintf(buf, "%s %d", longName, moleculeID);
	return buf;
}

/**
 * void Molecule::setID(int)
 *
 * Set the ID of a molecule, used when displaying molecule names. The ID is a number which is not necessarily unique, but should only be shared
 * between strongly related molecules.
 *
 */
void Molecule::setID(int i){

	moleculeID = i;
}

/**
 * int Molecule::getID()
 *
 * Get the ID of the current molecule.
 *
 * @return the current Molecule's ID
 */
int Molecule::getID(){
	return moleculeID;
}

/**
 * void Molecule::reset()
 *
 * Reset the molecule between Runge-Kutta runs. The vector containing runge-kutta data points is erased, the initial
 * concentration is added as the first element, and the rkVals are all reset to 0.
 *
 */
void Molecule::reset(){
	
	rungeKuttaSolution.erase(rungeKuttaSolution.begin(), rungeKuttaSolution.end());
	rungeKuttaSolution.push_back(initialConcentration);
	
	minima.erase(minima.begin(), minima.end());

	maxima.erase(maxima.begin(), maxima.end());

	currentConcentration = initialConcentration;
	rkVal[0] = 0;
	rkVal[1] = 0;
	rkVal[2] = 0;
	rkVal[3] = 0;
	prevDir = 0;
	currentDir = 0;
	numChanges = 0;
}

/**
 *
 *
 *
 *
 */
vector<float>* Molecule::getRungeKuttaSolution(){
	return &rungeKuttaSolution;

}

/**
 *
 *
 *
 *
 */
int Molecule::getPTMCount(int index){
	return PTMArray[index];
}

/**
 *
 *
 *
 *
 */
int Molecule::getScore(){

	t.trace("test","%s numChanges: %d, minima: %d, maxima: %d\n",getShortName(),numChanges, minima.size(), maxima.size());
	return numChanges;

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



