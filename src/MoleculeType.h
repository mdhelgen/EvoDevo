/**
 * MoleculeType.h
 *
 * MoleculeType holds default information about a type of molecule.
 */

#ifndef MOLECULE_TYPE_H_
#define MOLECULE_TYPE_H_

class MoleculeType{

public:
	MoleculeType();
	~MoleculeType();



private:
	
	
	//Name associated with the type of molecule
	const char* name;

	//Short name associated with the type of molecule
	const char* shortName;

	//Default factor applied to incoming interaction effects
	float defaultIncomingFactor;

	//Default factor applied to outgoing interactino effects
	float defaultOutgoingFactor;

	//Maximum allowed instances of this molecule type allowed
	int maxAllowed;
	
	//Display this type of molecule in the output files
	bool showInOutput;

	//Calculate value updates during integration for this molecule type
	bool calcDuringRK;

	

}





#endif
